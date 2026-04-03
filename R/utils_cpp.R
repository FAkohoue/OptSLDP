# utils_cpp.R
# R-level wrappers around the compiled C++ LD kernel functions.
# Each wrapper validates inputs, handles the name-matching bookkeeping,
# and falls back gracefully to pure-R if the shared library is not loaded.

#' Candidate-row LD matrix via C++
#'
#' Computes r^2 only for selected candidate SNP rows against all SNPs in
#' `geno_mat`. Far cheaper than the full n x n matrix when n_cands << n_snps.
#'
#' @param geno_mat     Numeric matrix (SNPs x samples, 0/1/2/NA). Must have
#'   rownames set to SNP IDs.
#' @param candidate_ids Character vector of candidate SNP IDs. Must be a
#'   subset of `rownames(geno_mat)`.
#' @return Numeric matrix (n_cands x n_snps) with rownames = matched candidate
#'   IDs, colnames = all SNP IDs. Returns a 0-row matrix if none match.
#' @keywords internal
#' @noRd
.compute_ld_subset_cpp <- function(geno_mat, candidate_ids) {
  if (is.null(rownames(geno_mat)))
    stop("geno_mat must have rownames.", call. = FALSE)
  idx  <- match(candidate_ids, rownames(geno_mat))
  keep <- !is.na(idx)
  if (!any(keep))
    return(matrix(numeric(0L), nrow = 0L, ncol = nrow(geno_mat),
                  dimnames = list(character(0L), rownames(geno_mat))))
  r2_subset_cpp(geno_mat, as.integer(idx[keep]))
}


#' Sparse threshold scan on a candidate x SNPs r^2 matrix via C++
#'
#' @param r2_mat   Numeric matrix from `.compute_ld_subset_cpp()`.
#' @param threshold Minimum r^2 to include.
#' @return Integer matrix with columns `row` and `col` (1-based indices).
#' @keywords internal
#' @noRd
.above_threshold_subset <- function(r2_mat, threshold) {
  above_threshold_subset_cpp(r2_mat, threshold, exclude_diag = TRUE)
}


#' Greedy LD pruning on a square r^2 matrix via C++
#'
#' Falls back to a pure-R loop if the C++ function is unavailable.
#'
#' @param r2_mat   Square numeric matrix.
#' @param threshold Pruning r^2 threshold.
#' @return Logical vector: `TRUE` = retained, `FALSE` = pruned.
#' @keywords internal
#' @noRd
.greedy_prune_r2 <- function(r2_mat, threshold) {
  if (!is.matrix(r2_mat))
    stop("r2_mat must be a matrix.", call. = FALSE)
  storage.mode(r2_mat) <- "double"

  # C++ path
  if (is.function(tryCatch(greedy_prune_r2_cpp, error = function(e) NULL)))
    return(greedy_prune_r2_cpp(r2_mat, threshold))

  # Pure-R fallback
  n    <- nrow(r2_mat)
  keep <- rep(TRUE, n)
  for (i in seq_len(n)) {
    if (!keep[i]) next
    hi <- which(r2_mat[i, ] >= threshold)
    hi <- hi[hi > i]
    if (length(hi)) keep[hi] <- FALSE
  }
  keep
}
