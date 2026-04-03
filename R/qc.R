# ==============================================================================
# qc.R
# Quality-control functions: MAF filtering and high-LD pre-pruning.
# All functions are scale-aware: they accept an optional GDS context and
# dispatch to the GDS backend for large datasets.
# ==============================================================================


#' Compute allele frequency and minor allele frequency for every SNP
#'
#' For the in-memory path, estimates AF from the dosage matrix as
#' `rowMeans / 2`. For the GDS path, delegates to `snpgdsSNPRateFreq()` which
#' streams genotypes from disk without loading the full matrix.
#'
#' @param geno_mat Numeric genotype matrix (SNPs x samples, coded 0/1/2/NA)
#'   **or** `NULL` when `ctx$strategy == "gds"`.
#' @param ctx      Optional GDS context list. Default `NULL`.
#'
#' @return A `data.table` with columns `SNP`, `AF`, and `MAF`.
#'   When using the GDS path, `SNP` values are the SNP IDs stored in the file.
#' @export
#' @examples
#' \dontrun{
#' maf_dt <- compute_maf(geno_mat)
#' }
compute_maf <- function(geno_mat, ctx = NULL) {

  # -- GDS path -----------------------------------------------------------------
  if (!is.null(ctx) && identical(ctx$strategy, "gds")) {
    freq <- SNPRelate::snpgdsSNPRateFreq(
      ctx$genofile,
      with.id    = TRUE
    )
    af  <- freq$AlleleFreq
    maf <- pmin(af, 1 - af)
    return(data.table::data.table(SNP = freq$snp.id, AF = af, MAF = maf))
  }

  # -- In-memory path ------------------------------------------------------------
  if (is.null(geno_mat)) {
    stop("compute_maf: geno_mat is NULL and no GDS context was provided.",
         call. = FALSE)
  }
  af  <- rowMeans(geno_mat, na.rm = TRUE) / 2
  maf <- pmin(af, 1 - af)
  data.table::data.table(SNP = rownames(geno_mat), AF = af, MAF = maf)
}


#' Filter SNPs by minor allele frequency
#'
#' Removes SNPs whose MAF is below `maf_threshold` or is `NA`. When a GDS
#' context is active the filtering is performed entirely in the GDS layer
#' without loading genotypes into RAM; when in-memory, the matrix is subset
#' in place.
#'
#' @param snp_info      Marker metadata `data.table`.
#' @param geno_mat      Numeric genotype matrix **or** `NULL` (GDS path).
#' @param maf_threshold Minimum MAF. Default `0.05`.
#' @param ctx           Optional GDS context list. Default `NULL`.
#'
#' @return A named list:
#'   \describe{
#'     \item{`snp_info`}{Filtered metadata, preserving original order.}
#'     \item{`geno_mat`}{Filtered matrix (`NULL` for GDS path).}
#'     \item{`maf_table`}{`data.table` of SNP / AF / MAF for retained markers.}
#'   }
#' @export
filter_snps_by_maf <- function(snp_info, geno_mat, maf_threshold = 0.05,
                                ctx = NULL) {

  # -- GDS path -----------------------------------------------------------------
  if (!is.null(ctx) && identical(ctx$strategy, "gds")) {
    keep    <- .filter_maf_gds(ctx$genofile,
                                maf_min  = maf_threshold,
                                n_cores  = ctx$n_cores)
    maf_all <- compute_maf(NULL, ctx = ctx)
    snp_info_f <- data.table::copy(snp_info[snp_info$SNP %in% keep, , drop = FALSE])
    data.table::setorder(snp_info_f, CHR, POS)
    return(list(
      snp_info  = snp_info_f,
      geno_mat  = NULL,
      maf_table = maf_all[SNP %in% keep]
    ))
  }

  # -- In-memory path ------------------------------------------------------------
  maf_dt <- compute_maf(geno_mat)
  keep   <- maf_dt$SNP[!is.na(maf_dt$MAF) & maf_dt$MAF >= maf_threshold]

  snp_info_f <- data.table::copy(
    snp_info[snp_info$SNP %in% keep, , drop = FALSE]
  )
  snp_info_f <- snp_info_f[match(keep, snp_info_f$SNP), , drop = FALSE]
  geno_mat_f <- geno_mat[keep, , drop = FALSE]

  list(
    snp_info  = snp_info_f,
    geno_mat  = geno_mat_f,
    maf_table = maf_dt[maf_dt$SNP %in% keep, ]
  )
}


#' Chromosome-wise high-LD pre-pruning
#'
#' Performs a fast pre-pruning pass at a very high LD threshold (default
#' `r^2 >= 0.99`) to remove near-duplicate SNPs before the main SLDP pipeline.
#'
#' * **In-memory** (`strategy = "in_memory"` / `"chunked"`): greedy
#'   position-order loop using `compute_r2_subset()`.
#' * **GDS** (`strategy = "gds"`): uses `snpgdsLDpruning()` which processes
#'   genotypes in streaming blocks -- fast enough for 10M+ SNP panels.
#'
#' @param snp_info    Marker metadata `data.table`.
#' @param geno_mat    Numeric genotype matrix **or** `NULL` (GDS path).
#' @param r2_pre      High-LD removal threshold. Default `0.99`.
#' @param slide_max_bp Sliding window in bp for the GDS path. Default `500 000`.
#' @param ctx         Optional GDS context list. Default `NULL`.
#' @param verbose     Logical. Default `TRUE`.
#'
#' @return A named list with `snp_info`, `geno_mat` (or `NULL`), and
#'   `kept_snps`.
#' @export
preprune_high_ld <- function(snp_info, geno_mat,
                              r2_pre       = 0.99,
                              slide_max_bp = 500000L,
                              ctx          = NULL,
                              verbose      = TRUE) {

  data.table::setDT(snp_info)
  data.table::setorder(snp_info, CHR, POS)

  .report_progress(
    "Starting chromosome-wise pre-pruning at r^2 >= ", r2_pre,
    verbose = verbose
  )

  # -- GDS path -----------------------------------------------------------------
  if (!is.null(ctx) && identical(ctx$strategy, "gds")) {
    keep_all <- character(0L)

    for (chr in unique(snp_info$CHR)) {
      chr_snps <- snp_info$SNP[snp_info$CHR == chr]
      .report_progress(
        "  Pre-pruning chromosome ", chr,
        " (", length(chr_snps), " SNPs) [GDS]",
        verbose = verbose
      )
      kept <- .preprune_chr_gds(
        ctx$genofile, chr_snps,
        r2_pre       = r2_pre,
        slide_max_bp = slide_max_bp,
        n_cores      = ctx$n_cores
      )
      keep_all <- c(keep_all, kept)
    }

    snp_info_f <- snp_info[snp_info$SNP %in% keep_all, , drop = FALSE]
    snp_info_f <- snp_info_f[order(snp_info_f$CHR, snp_info_f$POS), , drop = FALSE]
    .report_progress(
      "Finished pre-pruning: retained ", nrow(snp_info_f), " SNPs",
      verbose = verbose
    )
    return(list(snp_info = snp_info_f, geno_mat = NULL,
                kept_snps = snp_info_f$SNP))
  }

  # -- In-memory path ------------------------------------------------------------
  keep_all <- character(0L)

  for (chr in unique(snp_info$CHR)) {
    chr_snps <- snp_info$SNP[snp_info$CHR == chr]
    .report_progress(
      "  Pre-pruning chromosome ", chr,
      " (", length(chr_snps), " SNPs)",
      verbose = verbose
    )

    if (length(chr_snps) <= 1L) {
      keep_all <- c(keep_all, chr_snps)
      next
    }

    r2_mat <- compute_r2_subset(geno_mat, chr_snps)
    kept   <- rep(TRUE, length(chr_snps))
    names(kept) <- chr_snps

    for (i in seq_along(chr_snps)) {
      if (!kept[chr_snps[i]]) next
      hi <- which(r2_mat[i, ] >= r2_pre)
      hi <- hi[hi > i]
      if (length(hi) > 0L) kept[chr_snps[hi]] <- FALSE
    }

    keep_all <- c(keep_all, names(kept)[kept])
    rm(r2_mat, kept)   # Inspiration 5: free r^2 matrix before next chromosome
    gc(FALSE)
  }

  snp_info_f <- snp_info[snp_info$SNP %in% keep_all, , drop = FALSE]
  snp_info_f <- snp_info_f[order(snp_info_f$CHR, snp_info_f$POS), , drop = FALSE]
  geno_mat_f <- geno_mat[snp_info_f$SNP, , drop = FALSE]

  .report_progress(
    "Finished pre-pruning: retained ", nrow(snp_info_f), " SNPs",
    verbose = verbose
  )

  list(
    snp_info  = snp_info_f,
    geno_mat  = geno_mat_f,
    kept_snps = snp_info_f$SNP
  )
}
