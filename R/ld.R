# ==============================================================================
# ld.R
# Scale-aware LD computation layer.
#
# Public functions (compute_r2_subset, find_ld_neighbors) are unchanged from
# the user's perspective. Internally they dispatch to either the fast in-memory
# tcrossprod() BLAS path or the GDS/SNPRelate streaming path depending on the
# active scale strategy stored in the GDS context object.
#
# In-memory LD computation uses tcrossprod() -> DGEMM (BLAS level-3) rather
# than stats::cor(use="pairwise.complete.obs"). This avoids cor()'s per-pair
# NA handling overhead and maps directly to the fastest multi-threaded BLAS
# routine. NAs are mean-imputed per SNP before the matrix multiply; zero-
# variance SNPs are flagged and returned as NA rows/columns.
#
# A "GDS context" is a list produced by .build_gds_context() in OptSLDP (sldp.R):
#   $strategy   : "in_memory" | "chunked" | "gds"
#   $geno_mat   : numeric matrix (NULL when strategy == "gds")
#   $genofile   : open GDS handle (NULL when strategy != "gds")
#   $gds_path   : path to main GDS file (NULL when strategy != "gds")
#   $n_cores    : integer
# ==============================================================================


# -- Internal: fast BLAS r^2 via tcrossprod() ----------------------------------

#' Compute r^2 matrix via mean-imputation + tcrossprod() (DGEMM)
#'
#' Faster than stats::cor(use="pairwise.complete.obs") because it avoids
#' per-pair NA handling. NAs are mean-imputed per SNP so that a single
#' tcrossprod() call covers all pairs simultaneously.
#' @importFrom stats sd
#' @keywords internal
#' @noRd
.r2_tcrossprod <- function(mat) {
  # Pure-R full r^2 matrix via mean-imputation + tcrossprod() (BLAS DGEMM).
  # The C++ kernel (r2_subset_cpp) computes candidate rows only and is called
  # from .compute_ld_subset_cpp() in utils_cpp.R. This function handles the
  # full square matrix needed by the in-memory and chunked strategies.
  n_samp <- ncol(mat)

  row_means <- rowMeans(mat, na.rm = TRUE)
  na_idx    <- which(is.na(mat), arr.ind = TRUE)
  if (nrow(na_idx) > 0L)
    mat[na_idx] <- row_means[na_idx[, 1L]]

  g_sd           <- apply(mat, 1L, stats::sd)
  zero_var       <- g_sd == 0 | is.na(g_sd)
  g_sd[zero_var] <- 1
  scale_factor   <- g_sd * sqrt(n_samp - 1L)
  mat_scaled     <- (mat - row_means) / scale_factor

  r2             <- tcrossprod(mat_scaled)^2
  r2[zero_var, ] <- NA_real_
  r2[, zero_var] <- NA_real_
  diag(r2)       <- 1
  r2
}


#' Compute pairwise LD (r^2) for a subset of SNPs
#'
#' Dispatches to an in-memory `tcrossprod()` BLAS computation or a disk-backed
#' `snpgdsLDMat()` call depending on the scale strategy embedded in `ctx`.
#' When no context is supplied the function falls back to the pure in-memory
#' path, preserving backward compatibility for direct calls.
#'
#' @param geno_mat  Numeric genotype matrix (SNPs x samples, coded 0/1/2/NA)
#'   **or** `NULL` when `ctx$strategy == "gds"`.  When provided directly
#'   (no `ctx`), this is the only data source used.
#' @param snp_ids   Character vector of SNP IDs. All IDs must be present in
#'   either `rownames(geno_mat)` or the GDS file.
#' @param ctx       Optional GDS context list produced by `.build_gds_context()`.
#'   When `NULL` (default) the function uses `geno_mat` directly.
#'
#' @return Symmetric numeric r^2 matrix with row and column names equal to
#'   `snp_ids`.
#' @export
#' @examples
#' \dontrun{
#' # Direct in-memory call (no context needed for small datasets)
#' r2 <- compute_r2_subset(geno_mat, snp_ids = c("rs1", "rs2", "rs3"))
#'
#' # Context-aware call (large datasets)
#' r2 <- compute_r2_subset(NULL, snp_ids = chr_snps, ctx = gds_ctx)
#' }
compute_r2_subset <- function(geno_mat, snp_ids, ctx = NULL) {

  # -- GDS path -----------------------------------------------------------------
  if (!is.null(ctx) && identical(ctx$strategy, "gds")) {
    return(.compute_r2_gds(ctx$genofile, snp_ids, n_cores = ctx$n_cores))
  }

  # -- In-memory path (in_memory or chunked) ------------------------------------
  mat  <- if (!is.null(ctx)) ctx$geno_mat else geno_mat
  if (is.null(mat)) {
    stop(
      "compute_r2_subset: geno_mat is NULL and no GDS context was provided.",
      call. = FALSE
    )
  }
  sub <- mat[snp_ids, , drop = FALSE]
  .r2_tcrossprod(sub)
}


#' Find LD neighbors for a focal SNP
#'
#' Returns all SNPs from `candidate_ids` whose pairwise r^2 with `focal_snp`
#' meets or exceeds `r2_threshold`. Scale-aware: uses GDS streaming when a
#' context with `strategy == "gds"` is supplied.
#'
#' @param focal_snp     SNP ID of the focal marker.
#' @param geno_mat      Numeric genotype matrix **or** `NULL` (GDS path).
#' @param candidate_ids Character vector of SNP IDs to evaluate (may include
#'   `focal_snp`; it is always excluded from the result).
#' @param r2_threshold  Minimum r^2 for neighbor classification. Default `0.90`.
#' @param ctx           Optional GDS context list. Default `NULL`.
#'
#' @return Character vector of neighbor SNP IDs (excluding `focal_snp`).
#' @export
#' @examples
#' \dontrun{
#' nbrs <- find_ld_neighbors("rs100", geno_mat, window_snps, r2_threshold = 0.90)
#' nbrs <- find_ld_neighbors("rs100", NULL,     window_snps, r2_threshold = 0.90,
#'                           ctx = gds_ctx)
#' }
find_ld_neighbors <- function(focal_snp, geno_mat, candidate_ids,
                              r2_threshold = 0.90, ctx = NULL) {
  all_ids <- unique(c(focal_snp, candidate_ids))

  # -- GDS path -----------------------------------------------------------------
  if (!is.null(ctx) && identical(ctx$strategy, "gds")) {
    r2_row <- .compute_r2_focal_gds(
      ctx$genofile, focal_snp, all_ids, n_cores = ctx$n_cores
    )
    nbrs <- names(r2_row[!is.na(r2_row) & r2_row >= r2_threshold])
    return(setdiff(nbrs, focal_snp))
  }

  # -- In-memory path ------------------------------------------------------------
  r2   <- compute_r2_subset(geno_mat, all_ids, ctx = ctx)
  nbrs <- names(which(r2[focal_snp, , drop = TRUE] >= r2_threshold))
  setdiff(nbrs, focal_snp)
}
