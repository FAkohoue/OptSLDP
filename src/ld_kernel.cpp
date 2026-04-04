// ld_kernel.cpp
// Candidate-aware LD kernel + threshold scan + greedy pruning for OptSLDP.
//
// Exported functions:
//   r2_subset_cpp()              -- r^2 for candidate rows vs all SNPs (cheaper
//                                   than full n x n matrix when n_cands << n)
//   above_threshold_subset_cpp() -- sparse (row,col) pairs >= threshold
//   greedy_prune_r2_cpp()        -- greedy forward-selection LD pruning in C++
//
// Why r2_subset_cpp instead of the full r2_matrix_cpp:
//   The full matrix is O(n^2 * m). When n_candidates << n_snps (e.g. 5000
//   candidates vs 50000 window SNPs), computing only the candidate rows is
//   O(n_cands * n_snps * m) -- much cheaper. For the expansion step we only
//   need the rows corresponding to candidate SNPs, not all pairs.

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;


// ----------------------------------------------------------------------------
// Internal: mean-impute, centre, and unit-norm scale rows of G.
// S satisfies: S * S^T = Pearson correlation matrix.
// Sets zero_var[i] = 1 for monomorphic (zero-SD) SNPs.
// ----------------------------------------------------------------------------
static arma::mat impute_scale_rows(const arma::mat& G, arma::uvec& zero_var) {
  arma::uword n_snps = G.n_rows;
  arma::uword n_samp = G.n_cols;

  arma::mat S(n_snps, n_samp);
  zero_var.set_size(n_snps);
  zero_var.zeros();

  for (arma::uword i = 0; i < n_snps; ++i) {
    arma::rowvec row = G.row(i);

    arma::uvec finite_idx = arma::find_finite(row);
    arma::uvec na_idx     = arma::find_nonfinite(row);

    double mu = finite_idx.n_elem > 0
    ? arma::mean(row.elem(finite_idx))
      : 0.0;

    for (arma::uword j : na_idx) row(j) = mu;
    row -= mu;

    double ss = arma::dot(row, row);
    if (ss < 1e-14) {
      zero_var(i) = 1;
      S.row(i).zeros();
    } else {
      S.row(i) = row / std::sqrt(ss);
    }
  }
  return S;
}


// ----------------------------------------------------------------------------
// r2_subset_cpp
//
// Compute r^2 only for selected candidate rows against ALL SNP rows.
// Much cheaper than the full n x n matrix when n_candidates << n_snps.
// Cost: O(n_cands * n_snps * n_samples) vs O(n_snps^2 * n_samples).
//
// @param G              NumericMatrix SNPs x samples (0/1/2, NA allowed).
//                       Must have rownames set to SNP IDs.
// @param candidate_idx  IntegerVector of 1-based row indices of candidate SNPs.
// @return               NumericMatrix (n_cands x n_snps):
//                         rownames = candidate SNP IDs
//                         colnames = all SNP IDs
//                         diagonal entries for self-pairs set to 1.0
//                         NA rows/cols for monomorphic SNPs
// ----------------------------------------------------------------------------
// [[Rcpp::export]]
NumericMatrix r2_subset_cpp(NumericMatrix G,
                            IntegerVector candidate_idx) {
  if (candidate_idx.size() == 0)
    return NumericMatrix(0, G.nrow());

  arma::mat  geno = as<arma::mat>(G);
  arma::uvec zero_var;
  arma::mat  S = impute_scale_rows(geno, zero_var);

  // Convert 1-based R indices to 0-based Armadillo indices
  arma::uvec cand(candidate_idx.size());
  for (int i = 0; i < candidate_idx.size(); ++i) {
    int idx = candidate_idx[i];
    if (idx < 1 || idx > G.nrow())
      stop("candidate_idx contains out-of-range values.");
    cand[i] = static_cast<arma::uword>(idx - 1);
  }

  // BLAS DGEMM: S_candidates * S_all^T  ->  r matrix (n_cands x n_snps)
  arma::mat C = S.rows(cand) * S.t();
  C = arma::square(C);   // r -> r^2

  // Mark monomorphic candidate rows as NA
  for (arma::uword i = 0; i < cand.n_elem; ++i) {
    if (zero_var[cand[i]]) C.row(i).fill(arma::datum::nan);
  }

  // Mark monomorphic SNP columns as NA
  for (arma::uword j = 0; j < zero_var.n_elem; ++j) {
    if (zero_var[j]) C.col(j).fill(arma::datum::nan);
  }

  // Set self-pairs to exactly 1.0
  for (arma::uword i = 0; i < cand.n_elem; ++i) {
    if (!zero_var[cand[i]]) C(i, cand[i]) = 1.0;
  }

  NumericMatrix out = wrap(C);
  CharacterVector all_rn = rownames(G);
  CharacterVector cand_rn(cand.n_elem);
  for (arma::uword i = 0; i < cand.n_elem; ++i)
    cand_rn[i] = all_rn[cand[i]];
  rownames(out) = cand_rn;
  colnames(out) = all_rn;
  return out;
}


// ----------------------------------------------------------------------------
// above_threshold_subset_cpp
//
// Scan a candidate x SNPs r^2 matrix and return all (row, col) pairs meeting
// the threshold. Excludes self-pairs (where row SNP ID == col SNP ID).
//
// @param R2           NumericMatrix (n_cands x n_snps) from r2_subset_cpp().
// @param threshold    double -- minimum r^2 to include (e.g. 0.90).
// @param exclude_diag bool   -- skip self-pairs by name match. Default true.
// @return             IntegerMatrix with columns "row" and "col" (1-based).
// ----------------------------------------------------------------------------
// [[Rcpp::export]]
IntegerMatrix above_threshold_subset_cpp(NumericMatrix R2,
                                         double threshold,
                                         bool   exclude_diag = true) {
  int nr = R2.nrow();
  int nc = R2.ncol();
  CharacterVector rn = rownames(R2);
  CharacterVector cn = colnames(R2);
  bool have_names = (rn.size() == nr && cn.size() == nc);

  std::vector<int> rows, cols;
  rows.reserve(static_cast<size_t>(nr) * 20);
  cols.reserve(static_cast<size_t>(nr) * 20);

  for (int i = 0; i < nr; ++i) {
    for (int j = 0; j < nc; ++j) {
      if (exclude_diag && have_names) {
        if (Rcpp::as<std::string>(rn[i]) == Rcpp::as<std::string>(cn[j]))
          continue;
      }
      double v = R2(i, j);
      if (!ISNAN(v) && v >= threshold) {
        rows.push_back(i + 1);   // 1-based for R
        cols.push_back(j + 1);
      }
    }
  }

  int np = static_cast<int>(rows.size());
  IntegerMatrix out(np, 2);
  for (int k = 0; k < np; ++k) {
    out(k, 0) = rows[k];
    out(k, 1) = cols[k];
  }
  colnames(out) = CharacterVector::create("row", "col");
  return out;
}


// ----------------------------------------------------------------------------
// greedy_prune_r2_cpp
//
// Forward-selection greedy LD pruning on a square r^2 matrix.
// For each SNP i (in order): if retained, mark all j > i with r^2 >= threshold
// as removed. Equivalent to PLINK's greedy pruning but runs entirely in C++.
//
// @param r2_mat    NumericMatrix -- square r^2 matrix (SNPs x SNPs).
// @param threshold double        -- pruning threshold (e.g. 0.80).
// @return          LogicalVector -- TRUE = retained, FALSE = pruned.
// ----------------------------------------------------------------------------
// [[Rcpp::export]]
LogicalVector greedy_prune_r2_cpp(const NumericMatrix& r2_mat,
                                  const double threshold) {
  const int n = r2_mat.nrow();
  if (r2_mat.ncol() != n)
    stop("r2_mat must be a square matrix.");

  LogicalVector keep(n, true);

  for (int i = 0; i < n; ++i) {
    if (!keep[i]) continue;
    for (int j = i + 1; j < n; ++j) {
      if (!keep[j]) continue;
      double val = r2_mat(i, j);
      if (!ISNAN(val) && val >= threshold)
        keep[j] = false;
    }
  }
  return keep;
}


// ============================================================================
// Screening kernel: vectorised OLS for all SNPs in a chunk in one C++ pass.
//
// screen_chunk_cpp()
//
// Computes beta, SE, z, p, PVE for every SNP row in G against a single
// pre-residualised phenotype vector y. Genotype variance (g_var) can be
// passed in if already computed (multi-trait reuse); pass an empty vector
// to have it computed here.
//
// @param G        NumericMatrix SNPs x samples (0/1/2, NA allowed).
//                 NAs are mean-imputed per SNP internally.
// @param y        NumericVector length n_samples (already residualised).
// @param g_var_in NumericVector length n_snps, or length-0 to compute here.
// @return         NumericMatrix n_snps x 6: beta, SE, z, p, PVE, g_var
//                 (g_var returned so caller can cache it for other traits)
// ============================================================================
// [[Rcpp::export]]
NumericMatrix screen_chunk_cpp(NumericMatrix G,
                               NumericVector y,
                               NumericVector g_var_in) {

  int n_snps = G.nrow();
  int n_samp = G.ncol();
  bool has_gvar = (g_var_in.size() == n_snps);

  // Output: beta, SE, z, p, PVE, g_var
  NumericMatrix out(n_snps, 6);
  colnames(out) = CharacterVector::create("beta","SE","z","p","PVE","g_var");

  // Precompute y statistics once
  double y_sum = 0.0;
  for (int j = 0; j < n_samp; ++j) y_sum += y[j];
  double y_bar = y_sum / n_samp;

  double ss_tot = 0.0;
  std::vector<double> y_dev(n_samp);
  for (int j = 0; j < n_samp; ++j) {
    y_dev[j] = y[j] - y_bar;
    ss_tot  += y_dev[j] * y_dev[j];
  }

  for (int i = 0; i < n_snps; ++i) {

    // -- Mean-impute NAs and compute g_bar ----------------------------------
    double g_sum   = 0.0;
    int    n_obs   = 0;
    for (int j = 0; j < n_samp; ++j) {
      double v = G(i, j);
      if (!ISNAN(v)) { g_sum += v; ++n_obs; }
    }
    if (n_obs < 5) {
      out(i, 0) = NA_REAL; out(i, 1) = NA_REAL;
      out(i, 2) = NA_REAL; out(i, 3) = NA_REAL;
      out(i, 4) = NA_REAL; out(i, 5) = NA_REAL;
      continue;
    }
    double g_bar_i = g_sum / n_obs;

    // Fill imputed row into local buffer
    std::vector<double> g(n_samp);
    for (int j = 0; j < n_samp; ++j) {
      double v = G(i, j);
      g[j] = ISNAN(v) ? g_bar_i : v;
    }

    // -- Variance and covariance --------------------------------------------
    double gvar_i;
    if (has_gvar) {
      gvar_i = g_var_in[i];
    } else {
      gvar_i = 0.0;
      for (int j = 0; j < n_samp; ++j) {
        double d = g[j] - g_bar_i;
        gvar_i  += d * d;
      }
      gvar_i /= (n_samp - 1.0);
    }

    if (gvar_i <= 0 || !R_FINITE(gvar_i)) {
      out(i, 0) = NA_REAL; out(i, 1) = NA_REAL;
      out(i, 2) = NA_REAL; out(i, 3) = NA_REAL;
      out(i, 4) = NA_REAL; out(i, 5) = gvar_i;
      continue;
    }

    double gcov_y = 0.0;
    for (int j = 0; j < n_samp; ++j)
      gcov_y += (g[j] - g_bar_i) * y_dev[j];
    gcov_y /= (n_samp - 1.0);

    // -- OLS statistics -----------------------------------------------------
    double beta_i  = gcov_y / gvar_i;
    double ss_reg  = beta_i * beta_i * gvar_i * (n_samp - 1.0);
    double rss_i   = ss_tot - ss_reg;
    int    df_res  = n_obs - 2;
    double sigma2  = rss_i / df_res;
    double se_i    = std::sqrt(sigma2 / (gvar_i * (n_samp - 1.0)));
    double z_i     = beta_i / se_i;

    // Two-tailed p from t distribution via R API
    double p_i     = 2.0 * ::Rf_pt(-std::fabs(z_i), df_res, 1, 0);
    double pve_i   = std::max(0.0, std::min(1.0, ss_reg / ss_tot));

    out(i, 0) = beta_i;
    out(i, 1) = se_i;
    out(i, 2) = z_i;
    out(i, 3) = p_i;
    out(i, 4) = pve_i;
    out(i, 5) = gvar_i;
  }
  return out;
}
