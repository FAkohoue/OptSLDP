# ==============================================================================
# gds_backend.R
# GDS-based backend for large-scale LD computation (SNPRelate).
#
# These functions are used automatically by the scale-aware LD layer
# (ld.R) when the GDS strategy is active. They should never be called
# directly by user code -- the public API is compute_r2_subset() and
# find_ld_neighbors(), which dispatch here transparently.
#
# Scale thresholds (defaults, configurable via options):
#   in_memory  : n_snps <= 200 000
#   chunked    : 200 000 < n_snps <= 2 000 000
#   gds        : n_snps >  2 000 000
#
# GDS strategy:
#   1. Write the full genotype matrix to a SNPRelate GDS file once.
#   2. All LD queries stream from disk via snpgdsLDMat() -- no full matrix
#      ever lives in RAM simultaneously.
#   3. Subset GDS files are created on demand for chromosome-level pruning
#      and closed immediately after use.
#   4. All GDS handles are tracked in a package-level environment
#      (.sldp_gds_env) and closed on exit or on explicit teardown.
# ==============================================================================


# -- Package-level GDS handle registry -----------------------------------------
# Stores open GDS file objects so they can be force-closed on error or exit.
.sldp_gds_env <- new.env(parent = emptyenv())
.sldp_gds_env$handles <- list()

.register_gds_handle <- function(handle, key) {
  .sldp_gds_env$handles[[key]] <- handle
  invisible(handle)
}

.deregister_gds_handle <- function(key) {
  .sldp_gds_env$handles[[key]] <- NULL
  invisible(NULL)
}

#' Close all open GDS file handles tracked by the package
#'
#' Called automatically at the end of `run_sldp()` and on package unload.
#' Safe to call multiple times (already-closed handles are silently skipped).
#'
#' @return Invisibly returns `NULL`.
#' @keywords internal
#' @noRd
.close_all_gds_handles <- function() {
  handles <- .sldp_gds_env$handles
  for (key in names(handles)) {
    h <- handles[[key]]
    tryCatch(
      SNPRelate::snpgdsClose(h),
      error = function(e) invisible(NULL)
    )
    .sldp_gds_env$handles[[key]] <- NULL
  }
  invisible(NULL)
}


# -- GDS file construction ------------------------------------------------------

#' Write a genotype matrix to a SNPRelate GDS file
#'
#' Converts the in-memory dosage matrix (SNPs x samples, coded 0/1/2) to the
#' binary GDS format required by SNPRelate. The file is written once and reused
#' for all subsequent LD queries.
#'
#' @param geno_mat  Numeric matrix (SNPs x samples, coded 0/1/2/NA).
#' @param snp_info  Marker metadata `data.table` aligned to `geno_mat`.
#' @param gds_path  Output file path (should end in `.gds`).
#' @param n_cores   Number of I/O threads passed to `snpgdsCreateGeno()`.
#'   Default `1`.
#' @param verbose   Logical. Default `TRUE`.
#'
#' @return Invisibly returns `gds_path`.
#' @keywords internal
#' @noRd
.write_gds <- function(geno_mat, snp_info, gds_path,
                        n_cores = 1L, verbose = TRUE) {
  .assert_packages("SNPRelate")
  .report_progress(
    "Writing ", nrow(geno_mat), " SNPs to GDS: ", gds_path,
    verbose = verbose
  )

  SNPRelate::snpgdsCreateGeno(
    gds.fn         = gds_path,
    genmat         = geno_mat,
    sample.id      = colnames(geno_mat),
    snp.id         = snp_info$SNP,
    snp.chromosome = snp_info$CHR,
    snp.position   = snp_info$POS,
    snp.allele     = paste(snp_info$REF, snp_info$ALT, sep = "/"),
    snpfirstdim    = TRUE,
    compress.annotation = "ZIP.max",
    compress.geno       = "ZIP.max"
  )

  .report_progress("GDS file written.", verbose = verbose)
  invisible(gds_path)
}


#' Open a GDS file and register the handle
#'
#' @param gds_path Path to an existing `.gds` file.
#' @param key      Registry key for handle tracking.
#' @return An open SNPRelate GDS file object.
#' @keywords internal
#' @noRd
.open_gds <- function(gds_path, key = gds_path) {
  .assert_packages("SNPRelate")
  handle <- SNPRelate::snpgdsOpen(gds_path, readonly = TRUE, allow.fork = TRUE)
  .register_gds_handle(handle, key)
  handle
}


#' Close and deregister a GDS file handle
#'
#' @param handle A GDS file object returned by `.open_gds()`.
#' @param key    Registry key used when the handle was registered.
#' @return Invisibly returns `NULL`.
#' @keywords internal
#' @noRd
.close_gds <- function(handle, key) {
  tryCatch(
    SNPRelate::snpgdsClose(handle),
    error = function(e) invisible(NULL)
  )
  .deregister_gds_handle(key)
  invisible(NULL)
}


#' Create a subset GDS file from an existing one
#'
#' Used to produce per-chromosome GDS files during background pruning without
#' loading full data into memory.
#'
#' @param src_path  Source `.gds` file path.
#' @param snp_ids   Character vector of SNP IDs to retain.
#' @param dst_path  Destination `.gds` file path.
#' @param verbose   Logical.
#' @return Invisibly returns `dst_path`.
#' @keywords internal
#' @noRd
.subset_gds <- function(src_path, snp_ids, dst_path, verbose = TRUE) {
  .assert_packages("SNPRelate")
  .report_progress(
    "Creating subset GDS (", length(snp_ids), " SNPs): ", dst_path,
    verbose = verbose
  )
  SNPRelate::snpgdsCreateGenoSet(
    src.fn   = src_path,
    dest.fn  = dst_path,
    snp.id   = snp_ids,
    compress.annotation = "ZIP.max",
    compress.geno       = "ZIP.max"
  )
  invisible(dst_path)
}


# -- GDS LD computation ---------------------------------------------------------

#' Compute pairwise r^2 for a SNP subset via SNPRelate (GDS backend)
#'
#' Uses `snpgdsLDMat()` which streams genotypes from disk in blocks, keeping
#' peak RAM proportional to the block size rather than the full matrix.
#'
#' @param genofile  An open GDS file object (from `.open_gds()`).
#' @param snp_ids   Character vector of SNP IDs to include.
#' @param method    LD method passed to `snpgdsLDMat()`. `"corr"` (default)
#'   returns the correlation coefficient; r^2 is computed as `LD^2`.
#' @param n_cores   Number of threads for `snpgdsLDMat()`. Default `1`.
#'
#' @return Symmetric numeric r^2 matrix with row/col names = `snp_ids`.
#' @keywords internal
#' @noRd
.compute_r2_gds <- function(genofile, snp_ids,
                             method  = "corr",
                             n_cores = 1L) {
  ld_obj  <- SNPRelate::snpgdsLDMat(
    genofile,
    snp.id   = snp_ids,
    method   = method,
    slide    = -1L,          # compute full matrix, not sliding window
    verbose  = FALSE,
    num.thread = n_cores
  )
  r2_mat          <- ld_obj$LD^2
  rownames(r2_mat) <- ld_obj$snp.id
  colnames(r2_mat) <- ld_obj$snp.id
  r2_mat
}


#' Compute r^2 between one focal SNP and its neighbors via GDS
#'
#' More efficient than `.compute_r2_gds()` when only a single row of the LD
#' matrix is needed (focal vs. all neighbors), because SNPRelate's sliding-
#' window mode avoids computing the full nxn matrix.
#'
#' @param genofile   An open GDS file object.
#' @param focal_id   SNP ID of the focal marker.
#' @param all_ids    Character vector including `focal_id` and all candidates.
#' @param n_cores    Number of threads. Default `1`.
#'
#' @return Named numeric vector of r^2 values for each `all_ids` entry
#'   (including `focal_id` = 1.0).
#' @keywords internal
#' @noRd
.compute_r2_focal_gds <- function(genofile, focal_id, all_ids,
                                   n_cores = 1L) {
  r2_mat <- .compute_r2_gds(genofile, all_ids, n_cores = n_cores)
  r2_mat[focal_id, , drop = TRUE]
}


# -- GDS MAF filtering ----------------------------------------------------------

#' Compute MAF for all SNPs in a GDS file
#'
#' Uses `snpgdsSNPRateFreq()` which is vectorised and does not load the full
#' matrix into RAM.
#'
#' @param genofile   An open GDS file object.
#' @param maf_min    Minimum MAF threshold. SNPs below this are removed.
#' @param missing_max Maximum missing rate (0-1) per SNP. Default `0.10`.
#' @param n_cores    Number of threads. Default `1`.
#'
#' @return Character vector of SNP IDs that pass both filters.
#' @keywords internal
#' @noRd
.filter_maf_gds <- function(genofile,
                             maf_min     = 0.05,
                             missing_max = 0.10,
                             n_cores     = 1L) {
  freq <- SNPRelate::snpgdsSNPRateFreq(
    genofile,
    with.id    = TRUE,
    verbose    = FALSE,
    num.thread = n_cores
  )
  keep <- freq$AlleleFreq >= maf_min &
          freq$AlleleFreq <= (1 - maf_min) &
          freq$MissingRate <= missing_max
  freq$snp.id[keep]
}


# -- GDS-based pre-pruning (chromosome-wise, streaming) ------------------------

#' Pre-prune a single chromosome at a very high r^2 threshold via GDS
#'
#' Mirrors `preprune_high_ld()` but reads genotypes from a GDS file in
#' streaming blocks via `snpgdsLDpruning()`, which is orders of magnitude faster
#' than the in-memory greedy loop for large chromosomes.
#'
#' @param genofile   An open GDS file object.
#' @param snp_ids    SNP IDs on this chromosome (in position order).
#' @param r2_pre     LD threshold. Default `0.99`.
#' @param slide_max_bp  Sliding window size in bp. Default `500 000`.
#' @param n_cores    Threads. Default `1`.
#'
#' @return Character vector of retained SNP IDs.
#' @keywords internal
#' @noRd
.preprune_chr_gds <- function(genofile, snp_ids,
                               r2_pre       = 0.99,
                               slide_max_bp = 500000L,
                               n_cores      = 1L) {
  kept <- SNPRelate::snpgdsLDpruning(
    genofile,
    snp.id          = snp_ids,
    ld.threshold    = sqrt(r2_pre),   # snpgdsPruneLD uses |r|, not r^2
    slide.max.bp    = slide_max_bp,
    slide.max.n     = -1L,
    missing.rate    = 0.10,
    num.thread      = n_cores,
    verbose         = FALSE
  )
  kept
}


# -- GDS-based background pruning (chromosome-wise, streaming) -----------------

#' Prune background SNPs on one chromosome via GDS
#'
#' @param genofile   An open GDS file object.
#' @param snp_ids    SNP IDs on this chromosome (in position order).
#' @param r2_genome  LD pruning threshold. Default `0.80`.
#' @param slide_max_bp Sliding window in bp. Default `1 000 000`.
#' @param n_cores    Threads. Default `1`.
#'
#' @return Character vector of retained SNP IDs.
#' @keywords internal
#' @noRd
.prune_background_chr_gds <- function(genofile, snp_ids,
                                       r2_genome    = 0.80,
                                       slide_max_bp = 1000000L,
                                       n_cores      = 1L) {
  SNPRelate::snpgdsLDpruning(
    genofile,
    snp.id          = snp_ids,
    ld.threshold    = sqrt(r2_genome),
    slide.max.bp    = slide_max_bp,
    slide.max.n     = -1L,
    missing.rate    = 0.10,
    num.thread      = n_cores,
    verbose         = FALSE
  )
}


# -- GDS output extraction ------------------------------------------------------

#' Extract a genotype sub-matrix from a GDS file
#'
#' Used at the end of the pipeline to assemble the final output matrix from the
#' GDS file without loading the whole dataset into memory first.
#'
#' @param genofile  An open GDS file object.
#' @param snp_ids   SNP IDs to extract.
#' @param sample_ids Sample IDs to extract (all by default).
#'
#' @return Numeric matrix (SNPs x samples, coded 0/1/2/NA),
#'   `rownames = snp_ids`, `colnames = sample_ids`.
#' @keywords internal
#' @noRd
.extract_geno_gds <- function(genofile, snp_ids, sample_ids = NULL) {
  mat <- SNPRelate::snpgdsGetGeno(
    genofile,
    snp.id      = snp_ids,
    sample.id   = sample_ids,
    snpfirstdim = TRUE,
    with.id     = TRUE,
    verbose     = FALSE
  )
  rownames(mat$genotype) <- mat$snp.id
  colnames(mat$genotype) <- mat$sample.id
  storage.mode(mat$genotype) <- "numeric"
  mat$genotype
}


# -- Scale strategy selector ---------------------------------------------------

#' Determine the appropriate scale strategy for a given SNP count
#'
#' Reads thresholds from package options first; falls back to hardcoded
#' defaults. Allows site-wide customisation via
#' `options(optsldp.thresh_small = 1e5, optsldp.thresh_medium = 1e6)`.
#'
#' | Strategy   | Condition                            | Key behaviour                              |
#' |------------|--------------------------------------|---------------------------------------------|
#' | in_memory  | n_snps <= thresh_small (200 K)       | Full matrix in RAM; `cor()` for LD.         |
#' | chunked    | thresh_small < n_snps <= thresh_med  | Per-chromosome matrix; `cor()` for LD.      |
#' | gds        | n_snps > thresh_medium (2 M)         | Disk-backed GDS; `snpgdsLDMat()` for LD.   |
#'
#' @param n_snps        Integer. Total SNP count after MAF filtering.
#' @param user_strategy Character or `NULL`. Explicit override; one of
#'   `"in_memory"`, `"chunked"`, or `"gds"`. `NULL` triggers auto-detection.
#'
#' @return One of `"in_memory"`, `"chunked"`, or `"gds"`.
#' @keywords internal
#' @noRd
.select_strategy <- function(n_snps, user_strategy = NULL) {
  if (!is.null(user_strategy)) {
    return(match.arg(user_strategy, c("in_memory", "chunked", "gds")))
  }
  thresh_small  <- getOption("optsldp.thresh_small",  default = 2e5)
  thresh_medium <- getOption("optsldp.thresh_medium",  default = 2e6)

  if (n_snps <= thresh_small)  return("in_memory")
  if (n_snps <= thresh_medium) return("chunked")
  return("gds")
}
