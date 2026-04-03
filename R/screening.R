# ==============================================================================
# screening.R
# Marginal SNP screening and candidate selection.
#
# Changes in v4
# -------------
# compute_screening_stats():
#   Accepts y as either a numeric vector (single trait, backward compatible)
#   OR a numeric matrix (n_samples x n_traits).  Returns a named list of
#   per-trait data.tables when y is a matrix; returns a single data.table
#   when y is a vector, exactly as before.
#
# select_candidate_snps():
#   Unchanged public signature.  Accepts a single stats data.table.
#
# .screen_all_traits() [internal]:
#   Loops compute_screening_stats() + select_candidate_snps() over every
#   column of a phenotype matrix, then returns both the per-trait stats list
#   and the union of all candidate SNP sets.
# ==============================================================================


#' Compute marginal screening statistics
#'
#' Fits a marginal linear regression of the phenotype on each SNP.  When
#' covariates are supplied, `y` is first residualised on them so the scan is
#' covariate-adjusted.
#'
#' @section Multi-trait input:
#' When `y` is a numeric matrix (n_samples x n_traits) the function is called
#' once per column and returns a **named list** of per-trait `data.table`s.
#' Column names of `y` become the names of the list.  When `y` is a plain
#' numeric vector the function returns a single `data.table` as before
#' (backward compatible).
#'
#' @param geno_mat Numeric genotype matrix (SNPs x samples, coded 0/1/2/NA).
#' @param y        Numeric phenotype vector **or** numeric matrix
#'   (n_samples x n_traits) aligned to `colnames(geno_mat)`.
#' @param covar    Optional covariate matrix / data.frame (nrow == nrow(y)).
#'   Applied to every trait independently. Default `NULL`.
#' @param verbose  Print progress messages. Default `TRUE`.
#'
#' @return Single `data.table` (vector `y`) **or** named list of `data.table`s
#'   (matrix `y`).  Each table has columns: SNP, beta, SE, z_score, P.value,
#'   PVE, AF, MAF.
#' @export
#' @examples
#' \dontrun{
#' # Single trait
#' stats <- compute_screening_stats(geno_mat, y = pheno_vec)
#'
#' # Multiple traits -- returns named list
#' stats_list <- compute_screening_stats(geno_mat, y = pheno_matrix)
#' stats_list[["Yield"]]
#' }
compute_screening_stats <- function(geno_mat, y, covar = NULL,
                                    verbose = TRUE) {
  # -- Multi-trait dispatch ----------------------------------------------------
  if (is.matrix(y)) {
    if (nrow(y) != ncol(geno_mat))
      stop("nrow(y) must equal ncol(geno_mat) for matrix phenotype input.",
           call. = FALSE)
    trait_names <- colnames(y) %||% paste0("Trait", seq_len(ncol(y)))
    out_list <- vector("list", ncol(y))
    names(out_list) <- trait_names
    for (j in seq_len(ncol(y))) {
      .report_progress("  Screening trait '", trait_names[j], "' ...",
                       verbose = verbose)
      out_list[[j]] <- .screen_single_trait(geno_mat,
                                            y    = y[, j],
                                            covar  = covar,
                                            verbose = FALSE)
    }
    return(out_list)
  }

  # -- Single-trait path (backward compatible) ---------------------------------
  .screen_single_trait(geno_mat, y = y, covar = covar, verbose = verbose)
}


# -- Internal: single-trait scan ------------------------------------------------

#' Run the marginal regression scan for one trait
#'
#' Core engine of `compute_screening_stats()`.  Always receives a plain
#' numeric vector `y`.
#'
#' @keywords internal
#' @noRd
.screen_single_trait <- function(geno_mat, y, covar = NULL, verbose = TRUE) {
  if (length(y) != ncol(geno_mat))
    stop("Length of y must equal ncol(geno_mat).", call. = FALSE)
  y <- as.numeric(y)

  # Residualise on covariates once, before the SNP loop
  if (!is.null(covar)) {
    covar <- as.data.frame(covar)
    if (nrow(covar) != length(y))
      stop("nrow(covar) must equal length(y).", call. = FALSE)
    y <- stats::residuals(stats::lm(y ~ ., data = covar))
  }

  maf_dt <- compute_maf(geno_mat)
  n_snps <- nrow(geno_mat)
  n_samp <- ncol(geno_mat)
  snp_ids <- rownames(geno_mat)

  .report_progress("Computing screening statistics for ", n_snps, " SNPs",
                   verbose = verbose)

  # -- Vectorised OLS via matrix algebra --------------------------------------
  # For each SNP i: beta_i = cov(g_i, y) / var(g_i)
  # Computed simultaneously for all SNPs using matrix ops.
  # Missing values handled by mean-imputing each SNP before algebra,
  # then flagging SNPs with too much missingness as NA.

  # Count non-missing per SNP
  not_na <- !is.na(geno_mat)              # SNPs x samples logical
  n_obs  <- rowSums(not_na)              # number of complete obs per SNP

  # Mean-impute missing values for vectorised computation
  # (results for low-obs or zero-var SNPs will be overwritten with NA below)
  row_means <- rowMeans(geno_mat, na.rm = TRUE)
  g_imp     <- geno_mat
  na_idx    <- which(is.na(geno_mat), arr.ind = TRUE)
  if (nrow(na_idx) > 0L)
    g_imp[na_idx] <- row_means[na_idx[, 1L]]

  # y mean and deviations (y has no NA after residualisation)
  y_bar   <- mean(y)
  y_dev   <- y - y_bar                   # length n_samp

  # Per-SNP: g_bar, g_dev, var(g), cov(g,y)
  g_bar   <- row_means                   # n_snps
  # g_dev matrix: SNPs x samples (sweep row means)
  g_dev   <- sweep(g_imp, 1L, g_bar, `-`)   # SNPs x samples
  g_var   <- rowSums(g_dev^2) / (n_samp - 1L)          # n_snps
  g_cov_y <- (g_dev %*% y_dev)[, 1L] / (n_samp - 1L)  # n_snps

  beta  <- g_cov_y / g_var              # n_snps
  # RSS from beta directly (no need to materialise residual matrix)
  # RSS_i = sum((y - ybar)^2) - beta_i^2 * sum((g_i - gbar_i)^2)
  ss_tot  <- sum(y_dev^2)
  ss_reg  <- beta^2 * g_var * (n_samp - 1L)  # = beta * g_cov_y * (n-1)
  rss     <- ss_tot - ss_reg
  df_res  <- n_obs - 2L
  sigma2  <- rss / df_res
  se      <- sqrt(sigma2 / (g_var * (n_samp - 1L)))

  z_score <- beta / se
  p_value <- 2 * stats::pt(-abs(z_score), df = df_res)
  pve     <- ss_reg / ss_tot
  pve     <- pmax(0, pmin(1, pve))

  # Flag invalid SNPs (too few obs or zero variance)
  invalid <- n_obs < 5L | g_var <= 0 | !is.finite(g_var)
  beta[invalid]    <- NA_real_
  se[invalid]      <- NA_real_
  z_score[invalid] <- NA_real_
  p_value[invalid] <- NA_real_
  pve[invalid]     <- NA_real_

  data.table::data.table(
    SNP     = snp_ids,
    beta    = beta,
    SE      = se,
    z_score = z_score,
    P.value = p_value,
    PVE     = pve,
    AF      = maf_dt$AF,
    MAF     = maf_dt$MAF
  )
}


# -- Internal: multi-trait loop (used by run_sldp) -----------------------------

#' Run screening and candidate selection across all traits, returning the union
#'
#' Called by `run_sldp()` when more than one trait is requested.  For each
#' trait it:
#'   1. Calls `.screen_single_trait()` independently (separate residualisation,
#'      separate p-values).
#'   2. Calls `select_candidate_snps()` with the shared threshold settings.
#'   3. Accumulates candidate SNP sets.
#' Returns the per-trait stats list and the **union** of all candidate sets so
#' that important SNPs from every trait are protected in the final panel.
#'
#' @param geno_mat       Genotype matrix (SNPs x samples).
#' @param pheno_matrix   Numeric matrix (n_samples x n_traits), column names =
#'   trait names.
#' @param covar          Shared covariate data.frame or `NULL`.
#' @param mode           Selection mode passed to `select_candidate_snps()`.
#' @param pval_threshold Threshold passed to `select_candidate_snps()`.
#' @param z_threshold    Threshold passed to `select_candidate_snps()`.
#' @param pve_threshold  Threshold passed to `select_candidate_snps()`.
#' @param threshold_logic `"AND"` or `"OR"`.
#' @param verbose        Logical.
#'
#' @return List with:
#'   \describe{
#'     \item{`screening_stats`}{Named list of per-trait stats data.tables.}
#'     \item{`candidate_snps_per_trait`}{Named list of per-trait candidate
#'       SNP ID vectors.}
#'     \item{`candidate_snps_union`}{Character vector -- union across all
#'       traits.}
#'   }
#' @keywords internal
#' @noRd
.screen_all_traits <- function(geno_mat, pheno_matrix,
                               covar           = NULL,
                               mode            = "A",
                               pval_threshold  = NULL,
                               z_threshold     = NULL,
                               pve_threshold   = NULL,
                               threshold_logic = "AND",
                               n_cores         = 1L,
                               verbose         = TRUE) {
  trait_names <- colnames(pheno_matrix)
  n_traits    <- length(trait_names)

  # Run screening in parallel across traits when multiple cores available
  # Each trait gets its own vectorised matrix algebra pass
  n_workers <- min(n_traits, max(1L, n_cores))

  if (n_workers > 1L && .Platform$OS.type == "unix") {
    if (verbose)
      message("  Running ", n_traits, " trait screens in parallel (",
              n_workers, " workers) ...")
    cl <- parallel::makeCluster(n_workers, type = "FORK")
    on.exit(parallel::stopCluster(cl), add = TRUE)

    results <- parallel::parLapply(cl, seq_len(n_traits), function(j) {
      st <- .screen_single_trait(geno_mat, y = pheno_matrix[, j],
                                 covar = covar, verbose = FALSE)
      ca <- select_candidate_snps(st, mode = mode,
                                  pval_threshold  = pval_threshold,
                                  z_threshold     = z_threshold,
                                  pve_threshold   = pve_threshold,
                                  logic           = threshold_logic)
      list(stats = st, cands = ca)
    })
    stats_list <- stats::setNames(lapply(results, `[[`, "stats"), trait_names)
    cands_list <- stats::setNames(lapply(results, `[[`, "cands"), trait_names)
  } else {
    stats_list <- vector("list", n_traits)
    cands_list <- vector("list", n_traits)
    names(stats_list) <- names(cands_list) <- trait_names

    for (j in seq_len(n_traits)) {
      tn <- trait_names[j]
      .report_progress("  [", j, "/", n_traits, "] Screening trait '", tn,
                       "' ...", verbose = verbose)
      stats_list[[tn]] <- .screen_single_trait(
        geno_mat, y = pheno_matrix[, j], covar = covar, verbose = FALSE)
      cands_list[[tn]] <- select_candidate_snps(
        stats_list[[tn]], mode = mode, pval_threshold = pval_threshold,
        z_threshold = z_threshold, pve_threshold = pve_threshold,
        logic = threshold_logic)
      .report_progress("    -> ", length(cands_list[[tn]]),
                       " candidate SNP(s) for '", tn, "'", verbose = verbose)
    }
  }

  # Union across traits: every SNP important for any trait is protected
  union_cands <- unique(unlist(cands_list, use.names = FALSE))
  .report_progress("  Combined candidate SNPs (union across ",
                   n_traits, " traits): ", length(union_cands),
                   verbose = verbose)

  list(
    screening_stats        = stats_list,
    candidate_snps_per_trait = cands_list,
    candidate_snps_union   = union_cands
  )
}


#' Select candidate SNPs by screening threshold
#'
#' Applies one of three selection modes to a screening statistics table.
#'
#' \describe{
#'   \item{Mode A}{P-value only. Requires `pval_threshold`.}
#'   \item{Mode B}{Effect-size criteria. Requires `z_threshold` and/or
#'     `pve_threshold`.}
#'   \item{Mode C}{Combination of A and B criteria.}
#' }
#'
#' @param stats_dt       `data.table` from `compute_screening_stats()` for a
#'   **single trait**.
#' @param mode           `"A"` (default), `"B"`, or `"C"`.
#' @param pval_threshold P-value upper bound (modes A / C).
#' @param z_threshold    |z-score| lower bound (modes B / C).
#' @param pve_threshold  PVE lower bound (modes B / C).
#' @param logic          `"AND"` (default) or `"OR"`.
#'
#' @return Character vector of selected SNP IDs.
#' @export
select_candidate_snps <- function(stats_dt,
                                  mode            = c("A", "B", "C"),
                                  pval_threshold  = NULL,
                                  z_threshold     = NULL,
                                  pve_threshold   = NULL,
                                  logic           = c("AND", "OR")) {
  mode    <- match.arg(mode)
  logic   <- match.arg(logic)
  use_and <- identical(logic, "AND")
  keep    <- rep(FALSE, nrow(stats_dt))

  if (mode == "A") {
    if (is.null(pval_threshold))
      stop("mode = 'A' requires pval_threshold.", call. = FALSE)
    keep <- !is.na(stats_dt$P.value) & stats_dt$P.value <= pval_threshold
  }

  if (mode %in% c("B", "C")) {
    conds <- list()
    if (mode %in% c("A", "C") && !is.null(pval_threshold))
      conds <- c(conds, list(!is.na(stats_dt$P.value) &
                               stats_dt$P.value <= pval_threshold))
    if (!is.null(z_threshold))
      conds <- c(conds, list(!is.na(stats_dt$z_score) &
                               abs(stats_dt$z_score) >= z_threshold))
    if (!is.null(pve_threshold))
      conds <- c(conds, list(!is.na(stats_dt$PVE) &
                               stats_dt$PVE >= pve_threshold))
    if (length(conds) == 0L)
      stop("mode = '", mode, "' requires at least one threshold.", call. = FALSE)
    keep <- if (use_and) Reduce(`&`, conds) else Reduce(`|`, conds)
  }

  stats_dt$SNP[keep]
}
