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

  .report_progress(
    "Expanding ", n_cands,
    " candidate SNPs (window = +/-", window_kb, " kb",
    if (!is.null(ctx)) paste0("; strategy=", ctx$strategy) else "",
    ")",
    verbose = verbose
  )

  for (i in seq_along(candidate_snps)) {
    if (i %% 100L == 0L) {
      .report_progress(
        "  Important SNP expansion: ", i, " / ", n_cands,
        verbose = verbose
      )
    }

    snp <- candidate_snps[i]
    row <- snp_info[snp_info$SNP == snp, , drop = FALSE]
    if (nrow(row) == 0L) next

    chr <- row$CHR[1L]
    pos <- row$POS[1L]

    # Layer 1: positional window
    region <- snp_info$SNP[
      snp_info$CHR == chr &
      snp_info$POS >= (pos - window_bp) &
      snp_info$POS <= (pos + window_bp)
    ]
    important <- union(important, region)

    # Layer 2: LD neighbors within window (scale-aware)
    if (isTRUE(include_ld_neighbors) && length(region) > 1L) {
      ld_nbrs <- find_ld_neighbors(
        focal_snp     = snp,
        geno_mat      = geno_mat,
        candidate_ids = unique(region),
        r2_threshold  = r2_flag,
        ctx           = ctx
      )
      important <- union(important, ld_nbrs)
    }
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
    kept   <- rep(TRUE, length(chr_snps))
    names(kept) <- chr_snps

    for (i in seq_along(chr_snps)) {
      if (!kept[chr_snps[i]]) next
      hi <- which(r2_mat[i, ] >= r2_genome)
      hi <- hi[hi > i]
      if (length(hi) > 0L) kept[chr_snps[hi]] <- FALSE
    }

    kept_all <- c(kept_all, names(kept)[kept])
    rm(r2_mat, kept)   # Inspiration 5: free r^2 matrix before next chromosome
    gc(FALSE)
  }

  unique(kept_all)
}
