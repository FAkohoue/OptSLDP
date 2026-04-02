# ==============================================================================
# utils.R
# Internal helper utilities for the sldp package.
# ==============================================================================

# Suppress R CMD check notes about data.table column names used as symbols
# in := and setorder() calls, and about utils::head used internally.
#' @importFrom utils head
utils::globalVariables(c("SNP", "CHR", "POS", "REF", "ALT",
                         "neg_log_p", "is_qtn", "P.value"))

#' Null-coalescing operator
#'
#' Returns `x` if non-NULL, otherwise `y`.
#'
#' @param x Primary object.
#' @param y Fallback object.
#' @return `x` if non-NULL, otherwise `y`.
#' @keywords internal
#' @noRd
`%||%` <- function(x, y) if (is.null(x)) y else x


#' Validate that required packages are installed
#'
#' Throws an informative error listing every missing package in a single call,
#' so the user can install them all at once rather than hitting sequential
#' errors.
#'
#' @param pkgs Character vector of package names to check.
#' @return Invisibly returns `TRUE` when all packages are present.
#' @keywords internal
#' @noRd
.assert_packages <- function(pkgs) {
  missing_pkgs <- pkgs[
    !vapply(pkgs, requireNamespace, logical(1L), quietly = TRUE)
  ]
  if (length(missing_pkgs) > 0L) {
    stop(
      sprintf(
        "Missing required package(s): %s. Please install them first.",
        paste(missing_pkgs, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}


#' Validate the internal genotype object contract
#'
#' Checks that `snp_info` has the five required metadata columns, that
#' `geno_mat` is a numeric matrix with matching dimensions and aligned
#' rownames, and that SNP IDs are unique.
#'
#' @param snp_info Marker metadata `data.table` or `data.frame`.
#' @param geno_mat Numeric genotype matrix (SNPs x samples, coded 0/1/2/NA).
#' @return Invisibly returns `TRUE` on success.
#' @keywords internal
#' @noRd
.validate_genotype_object <- function(snp_info, geno_mat) {
  required_cols <- c("SNP", "CHR", "POS", "REF", "ALT")
  missing_cols  <- setdiff(required_cols, names(snp_info))
  if (length(missing_cols) > 0L) {
    stop(
      "snp_info is missing required columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }
  if (!is.matrix(geno_mat)) {
    stop("geno_mat must be a matrix.", call. = FALSE)
  }
  if (!is.numeric(geno_mat)) {
    stop("geno_mat must be numeric and coded as 0/1/2/NA.", call. = FALSE)
  }
  if (nrow(geno_mat) != nrow(snp_info)) {
    stop("nrow(geno_mat) must equal nrow(snp_info).", call. = FALSE)
  }
  if (is.null(rownames(geno_mat))) {
    stop("geno_mat must have rownames corresponding to SNP IDs.", call. = FALSE)
  }
  if (anyDuplicated(snp_info$SNP)) {
    stop("SNP IDs must be unique.", call. = FALSE)
  }
  if (!all(snp_info$SNP == rownames(geno_mat))) {
    stop(
      "snp_info$SNP must align exactly with rownames(geno_mat).",
      call. = FALSE
    )
  }
  invisible(TRUE)
}


#' Create an empty pruning statistics table
#'
#' Returns a zero-row `data.table` with the schema used by
#' `.append_pruning_stat()` and `write_pruning_report()`.
#'
#' @return An empty `data.table` with columns `step`, `n_before`, `n_after`,
#'   `n_removed`, and `details`.
#' @keywords internal
#' @noRd
.new_pruning_stats <- function() {
  data.table::data.table(
    step      = character(),
    n_before  = integer(),
    n_after   = integer(),
    n_removed = integer(),
    details   = character()
  )
}


#' Append one record to the pruning statistics table
#'
#' @param stats_dt Existing statistics `data.table`.
#' @param step     Character label for this pipeline step.
#' @param n_before SNP count entering the step.
#' @param n_after  SNP count leaving the step.
#' @param details  Optional free-text annotation.
#' @return Updated `data.table` with the new row appended.
#' @keywords internal
#' @noRd
.append_pruning_stat <- function(stats_dt, step, n_before, n_after,
                                 details = "") {
  data.table::rbindlist(
    list(
      stats_dt,
      data.table::data.table(
        step      = step,
        n_before  = as.integer(n_before),
        n_after   = as.integer(n_after),
        n_removed = as.integer(n_before - n_after),
        details   = as.character(details)
      )
    ),
    use.names = TRUE,
    fill      = TRUE
  )
}


#' Print a timestamped progress message
#'
#' Writes a single line to the console only when `verbose = TRUE`. The
#' timestamp format is `YYYY-MM-DD HH:MM:SS`.
#'
#' @param ... Message fragments passed to `paste0()`.
#' @param verbose Logical; set to `FALSE` to suppress all output.
#' @return Invisibly returns `NULL`.
#' @keywords internal
#' @noRd
.report_progress <- function(..., verbose = TRUE) {
  if (isTRUE(verbose)) {
    cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                paste0(...)))
  }
  invisible(NULL)
}
