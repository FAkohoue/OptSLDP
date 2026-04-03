# ==============================================================================
# pruning.R
# Core SLDP pruning: important SNP expansion and background pruning.
# Both functions are scale-aware and accept an optional GDS context.
# ==============================================================================


#' Expand candidate SNPs into the full set of important (protected) SNPs
#'
#' Builds the protected marker set in two layers:
#' 1. **Positional window** -- every SNP within `window_kb` kb of each candidate
#'    (same chromosome) is added.
#' 2. **LD neighborhood** (optional) -- SNPs in the window with r^2 >= `r2_flag`
#'    against the focal candidate are also added.
#'
#' The LD computation dispatches to the GDS backend automatically when a GDS
#' context (`ctx`) is supplied.
#'
#' @param candidate_snps       Character vector of candidate SNP IDs.
#' @param snp_info             Marker metadata `data.table`.
#' @param geno_mat             Numeric genotype matrix **or** `NULL` (GDS path).
#' @param window_kb            Positional expansion window in kb (each side).
#'   Default `50`.
#' @param include_ld_neighbors Logical; include LD neighbors. Default `TRUE`.
#' @param r2_flag              r^2 threshold for neighbor inclusion. Default `0.90`.
#' @param ctx                  Optional GDS context list. Default `NULL`.
#' @param verbose              Logical. Default `TRUE`.
#'
#' @return Character vector of unique important SNP IDs (superset of
#'   `candidate_snps`).
#' @export
expand_important_snps <- function(candidate_snps,
                                  snp_info,
                                  geno_mat,
                                  window_kb            = 50,
                                  include_ld_neighbors = TRUE,
                                  r2_flag              = 0.90,
                                  ctx                  = NULL,
                                  verbose              = TRUE) {

  data.table::setDT(snp_info)
  important <- unique(candidate_snps)
  window_bp <- as.integer(window_kb * 1000L)
  n_cands   <- length(candidate_snps)
  use_gds   <- !is.null(ctx) && identical(ctx$strategy, "gds")

  .report_progress(
    "Expanding ", n_cands,
    " candidate SNPs (window = +/-", window_kb, " kb",
    if (use_gds) "; batched GDS" else "",
    ")",
    verbose = verbose
  )

  # --------------------------------------------------------------------------
  # BATCHED CHROMOSOME STRATEGY (GDS path)
  # --------------------------------------------------------------------------
  # Instead of calling snpgdsLDMat() once per candidate (~5000+ GDS reads),
  # we group candidates by chromosome and:
  #   1. Find the union of all candidate windows on that chromosome.
  #   2. Extract that union as a genotype matrix ONCE via .extract_geno_gds().
  #   3. Compute the full LD matrix ONCE in RAM with cor().
  #   4. Slice it per candidate -- no further GDS reads needed.
  # This reduces GDS reads from n_candidates to n_chromosomes (~11 for rice).
  # --------------------------------------------------------------------------
  if (use_gds && isTRUE(include_ld_neighbors)) {
    # -- Build a keyed data.table for fast positional queries ------------------
    si_key <- data.table::copy(snp_info[, c("SNP","CHR","POS"), with = FALSE])
    data.table::setkey(si_key, CHR, POS)

    cands_dt <- merge(
      data.table::data.table(SNP = candidate_snps),
      si_key, by = "SNP", sort = FALSE
    )

    # Accumulate important SNPs as a set (integer index into si_key)
    important_set <- data.table::data.table(SNP = important)

    for (chr in unique(cands_dt$CHR)) {
      chr_cands <- cands_dt[cands_dt$CHR == chr, ]
      chr_info  <- si_key[si_key$CHR == chr, ]
      n_chr_cands <- nrow(chr_cands)

      .report_progress(
        "  LD expansion chromosome ", chr,
        " (", n_chr_cands, " candidates, ",
        nrow(chr_info), " chr SNPs) ...",
        verbose = verbose
      )

      # ---- Layer 1: all positional windows vectorised ----------------------
      # For each candidate, find window SNPs using a range join.
      # data.table foverlaps() does this in O(n log n) instead of O(n * m).
      windows_dt <- data.table::data.table(
        CHR   = chr,
        start = chr_cands$POS - window_bp,
        end   = chr_cands$POS + window_bp
      )
      chr_pos_dt <- data.table::data.table(
        CHR   = chr,
        start = chr_info$POS,
        end   = chr_info$POS,
        SNP   = chr_info$SNP
      )
      data.table::setkey(windows_dt,  CHR, start, end)
      data.table::setkey(chr_pos_dt,  CHR, start, end)
      overlaps   <- data.table::foverlaps(chr_pos_dt, windows_dt,
                                          type = "within", nomatch = NULL)
      window_snp_ids <- unique(overlaps$SNP)
      important_set  <- unique(rbind(important_set,
                                     data.table::data.table(SNP = window_snp_ids)))
      rm(overlaps, windows_dt, chr_pos_dt); gc(FALSE)

      if (length(window_snp_ids) <= 1L) next

      # ---- Layer 2: extract genotypes once per chromosome ------------------
      geno_chr <- tryCatch(
        .extract_geno_gds(ctx$genofile,
                          snp_ids    = window_snp_ids,
                          sample_ids = NULL),
        error = function(e) NULL
      )
      if (is.null(geno_chr)) next

      # ---- Layer 3: candidate-subset LD via C++ r2_subset_cpp() --------------
      # Only compute r^2 for candidate rows vs all window SNPs.
      # Cost: O(n_cands * n_window * n_samp) vs O(n_window^2 * n_samp).
      # Uses BLAS DGEMM on the candidate sub-rows only.
      cand_in_window <- intersect(chr_cands$SNP, rownames(geno_chr))
      if (length(cand_in_window) == 0L) { rm(geno_chr); gc(FALSE); next }

      r2_chr <- tryCatch(
        .compute_ld_subset_cpp(geno_chr, cand_in_window),
        error = function(e) {
          # Pure-R fallback: full matrix then slice candidate rows
          .r2_tcrossprod(geno_chr)[cand_in_window, , drop = FALSE]
        }
      )
      rm(geno_chr); gc(FALSE)

      # ---- Layer 4: C++ sparse threshold filter on candidate-subset matrix ---
      # r2_chr is already (cands x all_window_snps) from Layer 3.
      # Use above_threshold_subset_cpp() which handles the asymmetric shape.
      above <- tryCatch(
        .above_threshold_subset(r2_chr, threshold = r2_flag),
        error = function(e) {
          idx <- which(r2_chr >= r2_flag, arr.ind = TRUE)
          idx[rownames(r2_chr)[idx[,1L]] != colnames(r2_chr)[idx[,2L]], ,
              drop = FALSE]
        }
      )

      if (nrow(above) > 0L) {
        rn <- rownames(r2_chr)
        cn <- colnames(r2_chr)
        ld_pairs <- data.table::data.table(
          focal    = rn[above[, 1L]],
          neighbor = cn[above[, 2L]]
        )
        if (nrow(ld_pairs) > 0L)
          important_set <- unique(rbind(important_set,
                                        data.table::data.table(SNP = ld_pairs$neighbor)))
        rm(ld_pairs)
      }
      rm(above, r2_chr); gc(FALSE)
    }

    return(unique(important_set$SNP))
  }

  # --------------------------------------------------------------------------
  # ORIGINAL PER-SNP LOOP (in-memory path or LD neighbors disabled)
  # --------------------------------------------------------------------------
  for (i in seq_along(candidate_snps)) {
    if (i %% 100L == 0L) {
      .report_progress(
        "  Important SNP expansion: ", i, " / ", n_cands,
        verbose = verbose
      )
      gc(FALSE)
    }

    snp <- candidate_snps[i]
    row <- snp_info[snp_info$SNP == snp, , drop = FALSE]
    if (nrow(row) == 0L) next

    chr <- row$CHR[1L]
    pos <- row$POS[1L]

    region <- snp_info$SNP[
      snp_info$CHR == chr &
        snp_info$POS >= (pos - window_bp) &
        snp_info$POS <= (pos + window_bp)
    ]
    important <- union(important, region)

    if (isTRUE(include_ld_neighbors) && length(region) > 1L) {
      ld_nbrs <- find_ld_neighbors(
        focal_snp     = snp,
        geno_mat      = geno_mat,
        candidate_ids = unique(region),
        r2_threshold  = r2_flag,
        ctx           = ctx
      )
      important <- union(important, ld_nbrs)
      rm(ld_nbrs)
    }
    rm(region)
  }

  unique(important)
}


#' LD-prune background (non-important) SNPs chromosome-wise
#'
#' Applies greedy LD pruning to all SNPs outside the protected set.
#'
#' * **In-memory** (`strategy = "in_memory"` / `"chunked"`): loads the
#'   chromosome sub-matrix and applies a forward-selection greedy loop.
#' * **GDS** (`strategy = "gds"`): delegates to `snpgdsLDpruning()` which
#'   processes genotypes in streaming blocks -- the chromosome sub-matrix is
#'   never fully materialised in RAM.
#'
#' @param remaining_snps Character vector of non-protected SNP IDs to prune.
#' @param snp_info       Marker metadata `data.table`.
#' @param geno_mat       Numeric genotype matrix **or** `NULL` (GDS path).
#' @param r2_genome      Pruning r^2 threshold. Default `0.80`.
#' @param slide_max_bp   Sliding window in bp for the GDS path. Default
#'   `1 000 000`.
#' @param ctx            Optional GDS context list. Default `NULL`.
#' @param verbose        Logical. Default `TRUE`.
#'
#' @return Character vector of retained SNP IDs.
#' @export
prune_background_snps <- function(remaining_snps,
                                  snp_info,
                                  geno_mat,
                                  r2_genome    = 0.80,
                                  slide_max_bp = 1000000L,
                                  ctx          = NULL,
                                  verbose      = TRUE) {

  data.table::setDT(snp_info)
  rem_info <- snp_info[snp_info$SNP %in% remaining_snps, , drop = FALSE]
  rem_info <- rem_info[order(rem_info$CHR, rem_info$POS), , drop = FALSE]
  kept_all <- character(0L)

  .report_progress(
    "Pruning background SNPs chromosome-wise at r^2 >= ", r2_genome,
    if (!is.null(ctx)) paste0(" [strategy=", ctx$strategy, "]") else "",
    verbose = verbose
  )

  # -- GDS path -----------------------------------------------------------------
  if (!is.null(ctx) && identical(ctx$strategy, "gds")) {
    for (chr in unique(rem_info$CHR)) {
      chr_snps <- rem_info$SNP[rem_info$CHR == chr]
      .report_progress(
        "  Background pruning chromosome ", chr,
        " (", length(chr_snps), " SNPs) [GDS]",
        verbose = verbose
      )
      if (length(chr_snps) <= 1L) {
        kept_all <- c(kept_all, chr_snps)
        next
      }
      kept <- .prune_background_chr_gds(
        ctx$genofile, chr_snps,
        r2_genome    = r2_genome,
        slide_max_bp = slide_max_bp,
        n_cores      = ctx$n_cores
      )
      kept_all <- c(kept_all, kept)
    }
    return(unique(kept_all))
  }

  # -- In-memory path ------------------------------------------------------------
  for (chr in unique(rem_info$CHR)) {
    chr_snps <- rem_info$SNP[rem_info$CHR == chr]
    .report_progress(
      "  Background pruning chromosome ", chr,
      " (", length(chr_snps), " SNPs)",
      verbose = verbose
    )

    if (length(chr_snps) <= 1L) {
      kept_all <- c(kept_all, chr_snps)
      next
    }

    r2_mat <- compute_r2_subset(geno_mat, chr_snps)
    # C++ greedy pruning -- same algorithm as the R loop but compiled
    keep_lgl <- .greedy_prune_r2(r2_mat, threshold = r2_genome)
    kept_all <- c(kept_all, chr_snps[keep_lgl])
    rm(r2_mat, keep_lgl); gc(FALSE)
  }

  unique(kept_all)
}
