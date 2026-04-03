# utils_snprelate.R
# Version-safe wrapper for SNPRelate functions.
# Checks formal arguments before passing num.thread, rather than using
# tryCatch, which is unreliable for C-level errors.

#' Call a SNPRelate function, passing num.thread only if supported
#'
#' Inspects the function's formal argument list at call time. If num.thread
#' is a declared parameter it is passed; otherwise the call proceeds without
#' it. This is more reliable than tryCatch because some SNPRelate errors
#' originate at the C level and cannot be caught by R error handlers.
#'
#' @param fn      A SNPRelate function object.
#' @param ...     Arguments forwarded to fn.
#' @param n_cores Number of threads to request. Default 1L.
#' @keywords internal
#' @noRd
.snprelate_call <- function(fn, ..., n_cores = 1L) {
  if (!is.function(fn))
    stop("'fn' must be a function.", call. = FALSE)
  fmls <- formals(fn)
  if (!is.null(fmls) && "num.thread" %in% names(fmls)) {
    fn(..., num.thread = n_cores)
  } else {
    fn(...)
  }
}
