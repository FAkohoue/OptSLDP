# utils_snprelate.R
# Version-safe wrapper for SNPRelate functions that may or may not
# accept num.thread depending on the installed version.

#' Call a SNPRelate function with optional threading
#'
#' Tries calling `fn` with `num.thread = n_cores`. If that argument is
#' rejected (older SNPRelate versions), retries without it transparently.
#'
#' @param fn      A SNPRelate function object.
#' @param ...     Arguments forwarded to `fn`.
#' @param n_cores Number of threads to request. Default `1L`.
#' @keywords internal
#' @noRd
.snprelate_call <- function(fn, ..., n_cores = 1L) {
  tryCatch(
    fn(..., num.thread = n_cores),
    error = function(e) {
      msg <- conditionMessage(e)
      if (grepl("num[.]thread", msg, fixed = FALSE) ||
          grepl("unused argument", msg, fixed = FALSE))
        fn(...)
      else
        stop(e)
    }
  )
}
