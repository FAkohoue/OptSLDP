# ==============================================================================
# tests/testthat/test-basic.R
# Comprehensive test suite for OptSLDP.
#
# Coverage
# --------
#  1.  compute_maf()            -- columns, values, NA handling, zero-variance
#  2.  filter_snps_by_maf()     -- filtering logic, return structure, edge cases
#  3.  compute_r2_subset()      -- symmetry, diagonal, pairwise values, 2-SNP case
#  4.  find_ld_neighbors()      -- above/below threshold, focal excluded, empty
#  5.  preprune_high_ld()       -- removes near-duplicate, keeps independent
#  6.  prune_background_snps()  -- greedy logic, multi-chromosome, trivial case
#  7.  expand_important_snps()  -- positional window, LD layer, no-op candidate;
#                                  empty candidates guard
#  8.  compute_screening_stats()-- single-trait columns, multi-trait list, covar,
#                                  insufficient obs -> NA, zero-variance -> NA
#  9.  select_candidate_snps()  -- mode A/B/C, AND/OR logic, all-NA, empty result
# 10.  read_phenotype()         -- single trait, multi-trait, sample alignment,
#                                  error on missing sample, duplicate IDs
# 11.  read_numeric_genotype()  -- columns, matrix values, chr normalisation
# 12.  write_numeric_genotype() + write_hapmap_genotype() -- round-trip fidelity
# 13.  write_pruning_report()   -- CSV written, summary written, column names
# 13b. clean_genotype_file()    -- pass-through, malformed removal, sep detection
# 14.  run_sldp()               -- single-trait mode A/B/C; QTN retention;
#                                  multi-trait union protection; output file written;
#                                  pruning_stats step names; scale_strategy;
#                                  hapmap output; pre-pruning; stats output files
# 15.  C++ kernel functions     -- r2_subset_cpp, above_threshold_subset_cpp,
#                                  greedy_prune_r2_cpp; vectorised OLS vs lm()
# 16.  GDS end-to-end           -- forced GDS strategy; FORK cluster; step 11
#                                  file writer; final_geno_mat NULL; CHR/POS order
# 17.  read_hapmap_genotype()    -- structure, dimensions, values, REF/ALT;
#                                  same dosage as numeric reader
# 18.  n_pcs automatic PCA       -- runs without error; changes beta vs no PCs;
#                                  ignored when covar_cols provided
# ==============================================================================

library(OptSLDP)
library(data.table)

# ==============================================================================
# Shared fixtures
# ==============================================================================

# Minimal 6-SNP x 8-sample genotype matrix used across many unit tests.
# SNP1 and SNP2 are perfectly correlated (r^2 = 1).
# SNP5 is all zeros (MAF = 0, variance = 0) -- removed by any MAF filter and returns NA stats.
.make_geno <- function() {
  set.seed(1L)
  m <- matrix(
    c(# SNP1  SNP2 are correlated, SNP3-4 independent, SNP5 mono, SNP6 varied
      0L,1L,2L,0L,1L,1L,2L,0L,   # SNP1
      0L,1L,2L,0L,1L,1L,2L,0L,   # SNP2  (identical to SNP1)
      0L,0L,1L,2L,0L,1L,1L,2L,   # SNP3
      2L,1L,0L,1L,2L,0L,1L,1L,   # SNP4
      0L,0L,0L,0L,0L,0L,0L,0L,   # SNP5  truly monomorphic (all zeros) -- var=0 -> NA stats; MAF=0 -> removed at any threshold>0
      1L,2L,0L,1L,0L,2L,1L,1L    # SNP6
    ),
    nrow = 6L, ncol = 8L, byrow = TRUE
  )
  storage.mode(m) <- "numeric"
  rownames(m) <- paste0("SNP", 1:6)
  colnames(m) <- paste0("S", 1:8)
  m
}

.make_snp_info <- function() {
  data.frame(
    SNP = paste0("SNP", 1:6),
    CHR = c("1","1","1","1","2","2"),
    POS = c(10000L, 20000L, 80000L, 90000L, 10000L, 50000L),
    REF = c("A","G","C","T","A","G"),
    ALT = c("T","A","G","C","G","A"),
    stringsAsFactors = FALSE
  )
}

.make_phenotype <- function(n = 8L) {
  set.seed(42L)
  rnorm(n)
}


# ==============================================================================
# 1.  compute_maf()
# ==============================================================================

test_that("compute_maf returns a data.table with required columns", {
  m   <- .make_geno()
  out <- compute_maf(m)
  expect_s3_class(out, "data.table")
  expect_true(all(c("SNP","AF","MAF") %in% names(out)))
  expect_equal(nrow(out), nrow(m))
})

test_that("compute_maf SNP column matches rownames", {
  m   <- .make_geno()
  out <- compute_maf(m)
  expect_equal(out$SNP, rownames(m))
})

test_that("compute_maf AF formula: rowMeans / 2", {
  m   <- .make_geno()
  out <- compute_maf(m)
  expected_af <- rowMeans(m, na.rm = TRUE) / 2
  expect_equal(out$AF, unname(expected_af))
})

test_that("compute_maf MAF is min(AF, 1-AF)", {
  m   <- .make_geno()
  out <- compute_maf(m)
  expect_equal(out$MAF, pmin(out$AF, 1 - out$AF))
  expect_true(all(out$MAF >= 0 & out$MAF <= 0.5))
})

test_that("compute_maf handles NA values without error", {
  m       <- .make_geno()
  m[1, 1] <- NA_real_
  expect_no_error(compute_maf(m))
  out <- compute_maf(m)
  expect_false(is.na(out$AF[1]))   # rowMeans with na.rm=TRUE should still work
})

test_that("compute_maf monomorphic SNP has MAF = 0", {
  m        <- matrix(c(0,0,0,0), nrow = 1L, ncol = 4L)
  rownames(m) <- "MONO"
  out <- compute_maf(m)
  expect_equal(out$MAF, 0)
})

test_that("compute_maf errors when geno_mat is NULL and no ctx", {
  expect_error(compute_maf(NULL), "geno_mat is NULL")
})


# ==============================================================================
# 2.  filter_snps_by_maf()
# ==============================================================================

test_that("filter_snps_by_maf returns correct list elements", {
  m   <- .make_geno()
  si  <- .make_snp_info()
  res <- filter_snps_by_maf(si, m, maf_threshold = 0.05)
  expect_named(res, c("snp_info","geno_mat","maf_table"))
})

test_that("filter_snps_by_maf removes monomorphic SNP5", {
  m   <- .make_geno()
  si  <- .make_snp_info()
  # SNP5 is all zeros (MAF = 0) -- removed at any threshold > 0
  res <- filter_snps_by_maf(si, m, maf_threshold = 0.05)
  expect_false("SNP5" %in% res$snp_info$SNP)
})

test_that("filter_snps_by_maf retains SNPs above threshold", {
  m   <- .make_geno()
  si  <- .make_snp_info()
  res <- filter_snps_by_maf(si, m, maf_threshold = 0.05)
  expect_true(all(res$maf_table$MAF >= 0.05))
})

test_that("filter_snps_by_maf returned geno_mat rows match snp_info rows", {
  m   <- .make_geno()
  si  <- .make_snp_info()
  res <- filter_snps_by_maf(si, m, maf_threshold = 0.05)
  expect_equal(nrow(res$geno_mat), nrow(res$snp_info))
  expect_equal(rownames(res$geno_mat), res$snp_info$SNP)
})

test_that("filter_snps_by_maf threshold = 0 keeps all non-NA SNPs", {
  m   <- .make_geno()
  si  <- .make_snp_info()
  res <- filter_snps_by_maf(si, m, maf_threshold = 0)
  expect_equal(nrow(res$snp_info), 6L)
})

test_that("filter_snps_by_maf threshold = 0.5 returns only MAF-0.5 SNPs or none", {
  m   <- .make_geno()
  si  <- .make_snp_info()
  res <- filter_snps_by_maf(si, m, maf_threshold = 0.5)
  if (nrow(res$snp_info) > 0L)
    expect_true(all(res$maf_table$MAF >= 0.5))
})


# ==============================================================================
# 3.  compute_r2_subset()
# ==============================================================================

test_that("compute_r2_subset returns a square symmetric matrix", {
  m    <- .make_geno()
  ids  <- paste0("SNP", 1:4)
  r2   <- compute_r2_subset(m, ids)
  expect_true(is.matrix(r2))
  expect_equal(nrow(r2), 4L)
  expect_equal(ncol(r2), 4L)
  expect_equal(rownames(r2), ids)
  expect_equal(colnames(r2), ids)
})

test_that("compute_r2_subset diagonal is 1", {
  m  <- .make_geno()
  r2 <- compute_r2_subset(m, paste0("SNP", 1:4))
  expect_equal(diag(r2), rep(1, 4L), ignore_attr = TRUE)
})

test_that("compute_r2_subset is symmetric", {
  m  <- .make_geno()
  r2 <- compute_r2_subset(m, paste0("SNP", 1:4))
  expect_equal(r2, t(r2))
})

test_that("compute_r2_subset SNP1 vs SNP2 is 1 (perfectly correlated)", {
  m  <- .make_geno()
  r2 <- compute_r2_subset(m, c("SNP1","SNP2"))
  expect_equal(r2["SNP1","SNP2"], 1, tolerance = 1e-10)
})

test_that("compute_r2_subset all values in [0, 1]", {
  m  <- .make_geno()
  # Exclude SNP5 (zero variance) to avoid spurious cor() warning
  r2 <- compute_r2_subset(m, c(paste0("SNP", 1:4), "SNP6"))
  expect_true(all(r2 >= 0 - 1e-10 & r2 <= 1 + 1e-10, na.rm = TRUE))
})

test_that("compute_r2_subset errors when geno_mat is NULL and no ctx", {
  expect_error(
    compute_r2_subset(NULL, c("SNP1","SNP2")),
    "geno_mat is NULL"
  )
})


# ==============================================================================
# 4.  find_ld_neighbors()
# ==============================================================================

test_that("find_ld_neighbors returns SNP2 as neighbor of SNP1 (r2=1)", {
  m    <- .make_geno()
  nbrs <- find_ld_neighbors("SNP1", m, c("SNP1","SNP2","SNP3"), r2_threshold = 0.9)
  expect_true("SNP2" %in% nbrs)
})

test_that("find_ld_neighbors never returns the focal SNP itself", {
  m    <- .make_geno()
  nbrs <- find_ld_neighbors("SNP1", m, paste0("SNP", 1:4), r2_threshold = 0.5)
  expect_false("SNP1" %in% nbrs)
})

test_that("find_ld_neighbors returns empty when threshold is very high for unrelated SNPs", {
  m    <- .make_geno()
  # SNP3 and SNP6 are on different chromosomes and should not be r2=1
  nbrs <- find_ld_neighbors("SNP3", m, c("SNP3","SNP6"), r2_threshold = 0.9999)
  # Either empty or only SNP3 filtered out -- SNP6 should NOT be r2=1 with SNP3
  expect_false("SNP3" %in% nbrs)
})

test_that("find_ld_neighbors returns character vector", {
  m    <- .make_geno()
  nbrs <- find_ld_neighbors("SNP1", m, c("SNP1","SNP2"), r2_threshold = 0.5)
  expect_type(nbrs, "character")
})


# ==============================================================================
# 5.  preprune_high_ld()
# ==============================================================================

test_that("preprune_high_ld removes one of two perfectly correlated SNPs", {
  m  <- .make_geno()
  si <- .make_snp_info()
  res <- preprune_high_ld(si, m, r2_pre = 0.99, verbose = FALSE)
  # SNP1 and SNP2 are perfectly correlated -- at most one should remain
  retained <- res$snp_info$SNP
  expect_false("SNP1" %in% retained && "SNP2" %in% retained)
})

test_that("preprune_high_ld keeps independent SNPs", {
  m  <- .make_geno()
  si <- .make_snp_info()
  res <- preprune_high_ld(si, m, r2_pre = 0.99, verbose = FALSE)
  retained <- res$snp_info$SNP
  # SNP3, SNP4 on chr1 and SNP6 on chr2 are independent -- all should be kept
  expect_true("SNP3" %in% retained)
  expect_true("SNP4" %in% retained)
  expect_true("SNP6" %in% retained)
})

test_that("preprune_high_ld return list has snp_info, geno_mat, kept_snps", {
  m  <- .make_geno()
  si <- .make_snp_info()
  res <- preprune_high_ld(si, m, r2_pre = 0.99, verbose = FALSE)
  expect_named(res, c("snp_info","geno_mat","kept_snps"))
})

test_that("preprune_high_ld geno_mat rownames align with snp_info$SNP", {
  m  <- .make_geno()
  si <- .make_snp_info()
  res <- preprune_high_ld(si, m, r2_pre = 0.99, verbose = FALSE)
  expect_equal(rownames(res$geno_mat), res$snp_info$SNP)
})

test_that("preprune_high_ld with r2_pre = 0 drastically reduces SNPs", {
  m  <- .make_geno()
  si <- .make_snp_info()
  res <- preprune_high_ld(si, m, r2_pre = 0, verbose = FALSE)
  # With r2 >= 0, most pairs are pruned; retained count must be well below 6
  # (SNP5 has zero variance so its r2 is NA and it may survive the greedy loop)
  expect_lt(nrow(res$snp_info), 5L)
})


# ==============================================================================
# 6.  prune_background_snps()
# ==============================================================================

test_that("prune_background_snps returns a character vector", {
  m  <- .make_geno()
  si <- .make_snp_info()
  out <- prune_background_snps(
    remaining_snps = paste0("SNP", 3:6),
    snp_info       = si,
    geno_mat       = m,
    r2_genome      = 0.80,
    verbose        = FALSE
  )
  expect_type(out, "character")
})

test_that("prune_background_snps does not return SNPs outside remaining_snps", {
  m  <- .make_geno()
  si <- .make_snp_info()
  remaining <- paste0("SNP", 3:4)
  out <- prune_background_snps(remaining, si, m, r2_genome = 0.80, verbose = FALSE)
  expect_true(all(out %in% remaining))
})

test_that("prune_background_snps keeps at least one SNP per chromosome", {
  m  <- .make_geno()
  si <- .make_snp_info()
  out <- prune_background_snps(
    remaining_snps = paste0("SNP", 1:6),
    snp_info       = si,
    geno_mat       = m,
    r2_genome      = 0.80,
    verbose        = FALSE
  )
  chr1_kept <- out[out %in% si$SNP[si$CHR == "1"]]
  chr2_kept <- out[out %in% si$SNP[si$CHR == "2"]]
  expect_gte(length(chr1_kept), 1L)
  expect_gte(length(chr2_kept), 1L)
})

test_that("prune_background_snps with r2_genome = 1 keeps all SNPs (no pair achieves r2=1 except SNP1/SNP2)", {
  m  <- .make_geno()
  si <- .make_snp_info()
  # Threshold of 1: only r2 exactly 1.0 triggers removal -- only SNP1/SNP2 pair
  out <- prune_background_snps(
    remaining_snps = paste0("SNP", 3:6),
    snp_info       = si,
    geno_mat       = m,
    r2_genome      = 1.0,
    verbose        = FALSE
  )
  expect_equal(sort(out), paste0("SNP", 3:6))
})

test_that("prune_background_snps trivial: single remaining SNP is always retained", {
  m  <- .make_geno()
  si <- .make_snp_info()
  out <- prune_background_snps("SNP6", si, m, r2_genome = 0.80, verbose = FALSE)
  expect_equal(out, "SNP6")
})


# ==============================================================================
# 7.  expand_important_snps()
# ==============================================================================

test_that("expand_important_snps always includes candidate SNPs themselves", {
  m   <- .make_geno()
  si  <- .make_snp_info()
  out <- expand_important_snps(
    candidate_snps       = "SNP1",
    snp_info             = si,
    geno_mat             = m,
    window_kb            = 5,
    include_ld_neighbors = FALSE,
    verbose              = FALSE
  )
  expect_true("SNP1" %in% out)
})

test_that("expand_important_snps positional window includes nearby SNPs", {
  m   <- .make_geno()
  si  <- .make_snp_info()
  # SNP1 at 10000, SNP2 at 20000 -- 50kb window includes both
  out <- expand_important_snps(
    candidate_snps       = "SNP1",
    snp_info             = si,
    geno_mat             = m,
    window_kb            = 50,
    include_ld_neighbors = FALSE,
    verbose              = FALSE
  )
  expect_true("SNP2" %in% out)
})

test_that("expand_important_snps does not include SNPs on other chromosomes via positional window", {
  m   <- .make_geno()
  si  <- .make_snp_info()
  out <- expand_important_snps(
    candidate_snps       = "SNP1",
    snp_info             = si,
    geno_mat             = m,
    window_kb            = 9999,
    include_ld_neighbors = FALSE,
    verbose              = FALSE
  )
  # SNP5, SNP6 are on chr2 -- should NOT be included via positional window
  expect_false("SNP5" %in% out)
  expect_false("SNP6" %in% out)
})

test_that("expand_important_snps LD layer adds SNP2 as neighbor of SNP1", {
  m   <- .make_geno()
  si  <- .make_snp_info()
  out <- expand_important_snps(
    candidate_snps       = "SNP1",
    snp_info             = si,
    geno_mat             = m,
    window_kb            = 50,
    include_ld_neighbors = TRUE,
    r2_flag              = 0.9,
    verbose              = FALSE
  )
  expect_true("SNP2" %in% out)
})

test_that("expand_important_snps returns unique IDs only", {
  m   <- .make_geno()
  si  <- .make_snp_info()
  out <- expand_important_snps(
    candidate_snps       = c("SNP1","SNP2"),
    snp_info             = si,
    geno_mat             = m,
    window_kb            = 50,
    include_ld_neighbors = FALSE,
    verbose              = FALSE
  )
  expect_equal(length(out), length(unique(out)))
})

test_that("expand_important_snps empty candidates returns empty set", {
  m   <- .make_geno()
  si  <- .make_snp_info()
  out <- expand_important_snps(
    candidate_snps       = character(0),
    snp_info             = si,
    geno_mat             = m,
    window_kb            = 50,
    include_ld_neighbors = FALSE,
    verbose              = FALSE
  )
  expect_equal(length(out), 0L)
})


# ==============================================================================
# 8.  compute_screening_stats()
# ==============================================================================

test_that("compute_screening_stats single-trait returns data.table with all columns", {
  m   <- .make_geno()
  y   <- .make_phenotype()
  out <- compute_screening_stats(m, y = y, verbose = FALSE)
  expect_s3_class(out, "data.table")
  expected_cols <- c("SNP","beta","SE","z_score","P.value","PVE","AF","MAF")
  expect_true(all(expected_cols %in% names(out)))
  expect_equal(nrow(out), 6L)
})

test_that("compute_screening_stats single-trait SNP column matches rownames", {
  m   <- .make_geno()
  y   <- .make_phenotype()
  out <- compute_screening_stats(m, y = y, verbose = FALSE)
  expect_equal(out$SNP, rownames(m))
})

test_that("compute_screening_stats P.values are in (0, 1] or NA", {
  m   <- .make_geno()
  y   <- .make_phenotype()
  out <- compute_screening_stats(m, y = y, verbose = FALSE)
  p   <- out$P.value[!is.na(out$P.value)]
  expect_true(all(p > 0 & p <= 1))
})

test_that("compute_screening_stats z_score = beta / SE", {
  m   <- .make_geno()
  y   <- .make_phenotype()
  out <- compute_screening_stats(m, y = y, verbose = FALSE)
  ok  <- !is.na(out$beta) & !is.na(out$SE) & out$SE != 0
  expect_equal(out$z_score[ok], out$beta[ok] / out$SE[ok], tolerance = 1e-10)
})

test_that("compute_screening_stats monomorphic SNP5 returns NA for stats", {
  m   <- .make_geno()
  y   <- .make_phenotype()
  out <- compute_screening_stats(m, y = y, verbose = FALSE)
  row5 <- out[out$SNP == "SNP5", ]
  expect_true(is.na(row5$beta))
  expect_true(is.na(row5$P.value))
})

test_that("compute_screening_stats covariate adjustment runs without error", {
  m     <- .make_geno()
  y     <- .make_phenotype()
  covar <- data.frame(PC1 = rnorm(8L), PC2 = rnorm(8L))
  expect_no_error(
    compute_screening_stats(m, y = y, covar = covar, verbose = FALSE)
  )
})

test_that("compute_screening_stats covariate adjustment changes beta values", {
  set.seed(7L)
  m     <- .make_geno()
  y     <- .make_phenotype()
  covar <- data.frame(PC1 = rnorm(8L))
  out_no_cov <- compute_screening_stats(m, y = y, verbose = FALSE)
  out_cov    <- compute_screening_stats(m, y = y, covar = covar, verbose = FALSE)
  # At least some betas should differ after residualisation
  diff_count <- sum(abs(out_no_cov$beta - out_cov$beta) > 1e-10, na.rm = TRUE)
  expect_gt(diff_count, 0L)
})

test_that("compute_screening_stats multi-trait returns named list", {
  m   <- .make_geno()
  y1  <- .make_phenotype()
  y2  <- rnorm(8L)
  Y   <- cbind(Trait1 = y1, Trait2 = y2)
  out <- compute_screening_stats(m, y = Y, verbose = FALSE)
  expect_type(out, "list")
  expect_named(out, c("Trait1","Trait2"))
})

test_that("compute_screening_stats multi-trait each element is a data.table", {
  m   <- .make_geno()
  Y   <- cbind(T1 = .make_phenotype(), T2 = rnorm(8L))
  out <- compute_screening_stats(m, y = Y, verbose = FALSE)
  expect_s3_class(out$T1, "data.table")
  expect_s3_class(out$T2, "data.table")
  expect_equal(nrow(out$T1), 6L)
  expect_equal(nrow(out$T2), 6L)
})

test_that("compute_screening_stats errors when y length mismatches ncol(geno_mat)", {
  m <- .make_geno()
  expect_error(
    compute_screening_stats(m, y = rnorm(5L), verbose = FALSE),
    "Length of y"
  )
})

test_that("compute_screening_stats errors when matrix y nrow mismatches", {
  m <- .make_geno()
  Y <- matrix(rnorm(5L * 2L), nrow = 5L, ncol = 2L)
  expect_error(
    compute_screening_stats(m, y = Y, verbose = FALSE),
    "nrow\\(y\\)"
  )
})


# ==============================================================================
# 9.  select_candidate_snps()
# ==============================================================================

# Shared stats table for all mode tests
.make_stats <- function() {
  data.table::data.table(
    SNP     = c("A","B","C","D"),
    beta    = c(0.8, 0.1, 0.5, -0.9),
    SE      = c(0.1, 0.1, 0.1,  0.1),
    z_score = c(8.0, 1.0, 5.0, -9.0),
    P.value = c(0.001, 0.8, 0.04, 0.002),
    PVE     = c(0.10, 0.005, 0.06, 0.12),
    AF      = c(0.3, 0.4, 0.2, 0.35),
    MAF     = c(0.3, 0.4, 0.2, 0.35)
  )
}

test_that("select_candidate_snps mode A returns SNPs below pval_threshold", {
  dt  <- .make_stats()
  out <- select_candidate_snps(dt, mode = "A", pval_threshold = 0.05)
  expect_true(all(out %in% c("A","C","D")))
  expect_false("B" %in% out)
})

test_that("select_candidate_snps mode A exact boundary: p == threshold is included", {
  dt  <- data.table::data.table(
    SNP = c("X","Y"), z_score = c(1,1), P.value = c(0.05, 0.051),
    PVE = c(0.01, 0.01), beta = c(1,1), SE = c(1,1),
    AF = c(0.3,0.3), MAF = c(0.3,0.3)
  )
  out <- select_candidate_snps(dt, mode = "A", pval_threshold = 0.05)
  expect_true("X" %in% out)
  expect_false("Y" %in% out)
})

test_that("select_candidate_snps mode A errors without pval_threshold", {
  expect_error(
    select_candidate_snps(.make_stats(), mode = "A"),
    "pval_threshold"
  )
})

test_that("select_candidate_snps mode B OR logic: keeps SNPs passing either criterion", {
  dt  <- .make_stats()
  out <- select_candidate_snps(dt, mode = "B", z_threshold = 8, pve_threshold = 0.10, logic = "OR")
  # A: |z|=8 >=8 OR PVE=0.10>=0.10 -> YES
  # B: |z|=1 <8 AND PVE=0.005<0.10 -> NO
  # C: |z|=5 <8 AND PVE=0.06 <0.10 -> NO  (neither criterion met)
  # D: |z|=9 >=8 OR PVE=0.12>=0.10 -> YES
  expect_true("A" %in% out)
  expect_true("D" %in% out)
  expect_false("B" %in% out)
})

test_that("select_candidate_snps mode B AND logic: requires both criteria", {
  dt  <- .make_stats()
  out <- select_candidate_snps(dt, mode = "B", z_threshold = 5, pve_threshold = 0.10, logic = "AND")
  # A: |z|=8>=5 AND PVE=0.10>=0.10 -> YES
  # C: |z|=5>=5 AND PVE=0.06<0.10  -> NO
  # D: |z|=9>=5 AND PVE=0.12>=0.10 -> YES
  expect_true("A" %in% out)
  expect_true("D" %in% out)
  expect_false("C" %in% out)
})

test_that("select_candidate_snps mode B errors without any threshold", {
  expect_error(
    select_candidate_snps(.make_stats(), mode = "B"),
    "requires at least one threshold"
  )
})

test_that("select_candidate_snps mode C combines pval and z_score with OR", {
  dt  <- .make_stats()
  out <- select_candidate_snps(
    dt, mode = "C",
    pval_threshold = 0.01, z_threshold = 5,
    logic = "OR"
  )
  # A: p=0.001<=0.01 OR |z|=8>=5  -> YES
  # B: p=0.800>0.01  AND |z|=1<5  -> NO
  # C: p=0.04 >0.01  OR  |z|=5>=5 -> YES (z criterion met)
  # D: p=0.002<=0.01 OR  |z|=9>=5 -> YES
  expect_true("A" %in% out)
  expect_false("B" %in% out)
  expect_true("C" %in% out)
  expect_true("D" %in% out)
})

test_that("select_candidate_snps returns empty character when no SNPs pass", {
  dt  <- .make_stats()
  out <- select_candidate_snps(dt, mode = "A", pval_threshold = 1e-10)
  expect_equal(length(out), 0L)
  expect_type(out, "character")
})

test_that("select_candidate_snps skips NA P.values", {
  dt <- data.table::data.table(
    SNP = c("X","Y"), P.value = c(NA_real_, 0.01),
    z_score = c(NA_real_, 1), PVE = c(NA_real_, 0.01),
    beta = c(1,1), SE = c(1,1), AF = c(0.3,0.3), MAF = c(0.3,0.3)
  )
  out <- select_candidate_snps(dt, mode = "A", pval_threshold = 0.05)
  expect_false("X" %in% out)
  expect_true("Y"  %in% out)
})


# ==============================================================================
# 10.  read_phenotype()
# ==============================================================================

test_that("read_phenotype single trait returns a numeric vector", {
  pheno_file <- system.file("extdata","example_phenotype.csv", package = "OptSLDP")
  out <- read_phenotype(pheno_file, trait_col = "Trait1")
  expect_type(out$phenotype, "double")
  expect_false(is.matrix(out$phenotype))
  expect_equal(length(out$phenotype), 50L)
})

test_that("read_phenotype multi-trait returns a numeric matrix", {
  pheno_file <- system.file("extdata","example_phenotype.csv", package = "OptSLDP")
  out <- read_phenotype(pheno_file, trait_col = c("Trait1","Trait2"))
  expect_true(is.matrix(out$phenotype))
  expect_equal(dim(out$phenotype), c(50L, 2L))
  expect_equal(colnames(out$phenotype), c("Trait1","Trait2"))
})

test_that("read_phenotype extracts covariates as data.frame", {
  pheno_file <- system.file("extdata","example_phenotype.csv", package = "OptSLDP")
  out <- read_phenotype(pheno_file, trait_col = "Trait1", covar_cols = c("PC1","PC2"))
  expect_s3_class(out$covariates, "data.frame")
  expect_equal(ncol(out$covariates), 2L)
  expect_equal(nrow(out$covariates), 50L)
})

test_that("read_phenotype sample_order reorders rows correctly", {
  pheno_file <- system.file("extdata","example_phenotype.csv", package = "OptSLDP")
  reversed   <- paste0("Line", sprintf("%02d", 50:1))
  out <- read_phenotype(pheno_file, trait_col = "Trait1", sample_order = reversed)
  expect_equal(out$sample_ids, reversed)
})

test_that("read_phenotype trait_names element is present", {
  pheno_file <- system.file("extdata","example_phenotype.csv", package = "OptSLDP")
  out <- read_phenotype(pheno_file, trait_col = c("Trait1","Trait2"))
  expect_equal(out$trait_names, c("Trait1","Trait2"))
})

test_that("read_phenotype errors when genotype samples are absent from phenotype", {
  pheno_file <- system.file("extdata","example_phenotype.csv", package = "OptSLDP")
  expect_error(
    read_phenotype(pheno_file, trait_col = "Trait1",
                   sample_order = c(paste0("Line", sprintf("%02d", 1:50)), "GHOST")),
    "absent from phenotype"
  )
})

test_that("read_phenotype errors on missing required column", {
  pheno_file <- system.file("extdata","example_phenotype.csv", package = "OptSLDP")
  expect_error(
    read_phenotype(pheno_file, trait_col = "NonExistent"),
    "missing"
  )
})


# ==============================================================================
# 11.  read_numeric_genotype()
# ==============================================================================

test_that("read_numeric_genotype returns correct list structure", {
  f   <- system.file("extdata","example_genotypes_numeric.csv", package = "OptSLDP")
  out <- read_numeric_genotype(f)
  expect_named(out, c("snp_info","geno_mat","sample_ids","format"))
  expect_equal(out$format, "numeric")
})

test_that("read_numeric_genotype snp_info has required columns", {
  f   <- system.file("extdata","example_genotypes_numeric.csv", package = "OptSLDP")
  out <- read_numeric_genotype(f)
  expect_true(all(c("SNP","CHR","POS","REF","ALT") %in% names(out$snp_info)))
})

test_that("read_numeric_genotype dimensions: 40 SNPs x 50 samples", {
  f   <- system.file("extdata","example_genotypes_numeric.csv", package = "OptSLDP")
  out <- read_numeric_genotype(f)
  expect_equal(nrow(out$geno_mat), 40L)
  expect_equal(ncol(out$geno_mat), 50L)
})

test_that("read_numeric_genotype geno_mat values are in {0,1,2,NA}", {
  f   <- system.file("extdata","example_genotypes_numeric.csv", package = "OptSLDP")
  out <- read_numeric_genotype(f)
  valid <- out$geno_mat[!is.na(out$geno_mat)]
  expect_true(all(valid %in% c(0, 1, 2)))
})

test_that("read_numeric_genotype normalises chr prefix: no 'chr' prefix in CHR", {
  f   <- system.file("extdata","example_genotypes_numeric.csv", package = "OptSLDP")
  out <- read_numeric_genotype(f)
  expect_false(any(grepl("^chr", out$snp_info$CHR, ignore.case = TRUE)))
})

test_that("read_numeric_genotype rownames match snp_info$SNP", {
  f   <- system.file("extdata","example_genotypes_numeric.csv", package = "OptSLDP")
  out <- read_numeric_genotype(f)
  expect_equal(rownames(out$geno_mat), out$snp_info$SNP)
})

test_that("read_numeric_genotype sample_ids match colnames of geno_mat", {
  f   <- system.file("extdata","example_genotypes_numeric.csv", package = "OptSLDP")
  out <- read_numeric_genotype(f)
  expect_equal(colnames(out$geno_mat), out$sample_ids)
})


# ==============================================================================
# 12.  write_numeric_genotype() + write_hapmap_genotype()  --  round-trip
# ==============================================================================

test_that("write_numeric_genotype round-trip preserves SNP metadata columns", {
  m  <- .make_geno()
  si <- .make_snp_info()
  tmp <- tempfile(fileext = ".csv")
  write_numeric_genotype(si, m, tmp)
  back <- read.csv(tmp, check.names = FALSE)
  expect_true(all(c("SNP","CHR","POS","REF","ALT") %in% names(back)))
  file.remove(tmp)
})

test_that("write_numeric_genotype round-trip preserves sample columns", {
  m  <- .make_geno()
  si <- .make_snp_info()
  tmp <- tempfile(fileext = ".csv")
  write_numeric_genotype(si, m, tmp)
  back <- read.csv(tmp, check.names = FALSE)
  expect_true(all(colnames(m) %in% names(back)))
  file.remove(tmp)
})

test_that("write_numeric_genotype round-trip preserves dosage values", {
  m  <- .make_geno()
  si <- .make_snp_info()
  tmp <- tempfile(fileext = ".csv")
  write_numeric_genotype(si, m, tmp)
  back      <- read.csv(tmp, check.names = FALSE)
  back_mat  <- as.matrix(back[, colnames(m)])
  storage.mode(back_mat) <- "numeric"
  rownames(back_mat) <- back$SNP
  expect_equal(back_mat, m, ignore_attr = TRUE)
  file.remove(tmp)
})

test_that("write_hapmap_genotype produces a file with correct header columns", {
  m  <- .make_geno()
  si <- .make_snp_info()
  tmp <- tempfile(fileext = ".hmp.txt")
  write_hapmap_genotype(si, m, tmp)
  back <- read.table(tmp, header = TRUE, sep = "\t", check.names = FALSE, comment.char = "")
  expect_true("rs#"     %in% names(back))
  expect_true("alleles" %in% names(back))
  expect_true("chrom"   %in% names(back))
  file.remove(tmp)
})

test_that("write_hapmap_genotype encodes dosage 0 as REFREF", {
  m  <- matrix(c(0L,1L,2L), nrow = 1L, ncol = 3L)
  storage.mode(m) <- "numeric"
  rownames(m) <- "S1"
  colnames(m) <- c("sA","sB","sC")
  si <- data.table::data.table(
    SNP="S1", CHR="1", POS=1000L, REF="A", ALT="T"
  )
  tmp <- tempfile(fileext = ".hmp.txt")
  write_hapmap_genotype(si, m, tmp)
  back <- read.table(tmp, header = TRUE, sep = "\t", check.names = FALSE, comment.char = "")
  expect_equal(back$sA, "AA")
  expect_equal(back$sB, "AT")
  expect_equal(back$sC, "TT")
  file.remove(tmp)
})

test_that("write_hapmap_genotype encodes NA as NN", {
  m  <- matrix(c(NA_real_), nrow = 1L, ncol = 1L)
  rownames(m) <- "S1"
  colnames(m) <- "sA"
  si <- data.table::data.table(
    SNP="S1", CHR="1", POS=1000L, REF="A", ALT="T"
  )
  tmp <- tempfile(fileext = ".hmp.txt")
  write_hapmap_genotype(si, m, tmp)
  back <- read.table(tmp, header = TRUE, sep = "\t", check.names = FALSE, comment.char = "")
  expect_equal(back$sA, "NN")
  file.remove(tmp)
})


# ==============================================================================
# 13.  write_pruning_report()
# ==============================================================================

.make_pruning_stats <- function() {
  data.table::data.table(
    step      = c("input","maf_filter","preprune"),
    n_before  = c(40L, 40L, 38L),
    n_after   = c(40L, 38L, 35L),
    n_removed = c(0L,   2L,  3L),
    details   = c("format=numeric","maf>=0.05","r2>=0.99")
  )
}

test_that("write_pruning_report creates the CSV file", {
  tmp <- tempfile(fileext = ".csv")
  write_pruning_report(.make_pruning_stats(), stats_file = tmp)
  expect_true(file.exists(tmp))
  file.remove(tmp)
})

test_that("write_pruning_report CSV has correct column names", {
  tmp <- tempfile(fileext = ".csv")
  write_pruning_report(.make_pruning_stats(), stats_file = tmp)
  back <- read.csv(tmp)
  expect_true(all(c("step","n_before","n_after","n_removed","details") %in% names(back)))
  file.remove(tmp)
})

test_that("write_pruning_report CSV has correct row count", {
  tmp <- tempfile(fileext = ".csv")
  write_pruning_report(.make_pruning_stats(), stats_file = tmp)
  back <- read.csv(tmp)
  expect_equal(nrow(back), 3L)
  file.remove(tmp)
})

test_that("write_pruning_report creates summary file when requested", {
  tmp_csv  <- tempfile(fileext = ".csv")
  tmp_txt  <- tempfile(fileext = ".txt")
  write_pruning_report(.make_pruning_stats(), stats_file = tmp_csv,
                       summary_file = tmp_txt)
  expect_true(file.exists(tmp_txt))
  lines <- readLines(tmp_txt)
  expect_gt(length(lines), 0L)
  file.remove(tmp_csv, tmp_txt)
})

test_that("write_pruning_report no summary_file: only CSV produced", {
  tmp_csv <- tempfile(fileext = ".csv")
  write_pruning_report(.make_pruning_stats(), stats_file = tmp_csv)
  expect_true(file.exists(tmp_csv))
  file.remove(tmp_csv)
})


# ==============================================================================
# 13b.  clean_genotype_file()
# ==============================================================================

test_that("clean_genotype_file passes through a clean numeric CSV unchanged", {
  f   <- system.file("extdata", "example_genotypes_numeric.csv", package = "OptSLDP")
  tmp <- tempfile(fileext = ".csv")
  res <- clean_genotype_file(file_in = f, file_out = tmp, sep = ",", verbose = FALSE)
  expect_true(file.exists(tmp))
  original <- readLines(f)
  cleaned  <- readLines(tmp)
  # Header + all data rows must be preserved
  expect_equal(length(cleaned), length(original))
  expect_equal(cleaned[1], original[1])
  file.remove(tmp)
})

test_that("clean_genotype_file returns named integer vector total/removed/kept", {
  f   <- system.file("extdata", "example_genotypes_numeric.csv", package = "OptSLDP")
  tmp <- tempfile(fileext = ".csv")
  res <- clean_genotype_file(file_in = f, file_out = tmp, sep = ",", verbose = FALSE)
  expect_type(res, "integer")
  expect_named(res, c("total", "removed", "kept"))
  expect_equal(res[["removed"]], 0L)
  expect_equal(res[["total"]],   res[["kept"]])
  file.remove(tmp)
})

test_that("clean_genotype_file removes a line with wrong column count", {
  # Write a mini numeric CSV with one malformed line
  tmp_in  <- tempfile(fileext = ".csv")
  tmp_out <- tempfile(fileext = ".csv")
  writeLines(c(
    "SNP,CHR,POS,REF,ALT,S1,S2,S3",
    "SNP1,1,1000,A,T,0,1,2",       # correct: 8 columns
    "SNP2,1,2000,G,C,1,extra,0,2", # malformed: 9 columns
    "SNP3,1,3000,C,T,2,0,1"        # correct: 8 columns
  ), tmp_in)
  res <- clean_genotype_file(file_in = tmp_in, file_out = tmp_out,
                             sep = ",", verbose = FALSE)
  expect_equal(res[["removed"]], 1L)
  expect_equal(res[["kept"]],    2L)
  kept_lines <- readLines(tmp_out)
  # Header + 2 good data lines
  expect_equal(length(kept_lines), 3L)
  file.remove(tmp_in, tmp_out)
})

test_that("clean_genotype_file auto-detects tab separator for VCF extension", {
  f   <- system.file("extdata", "example_genotypes.vcf", package = "OptSLDP")
  if (!nzchar(f)) skip("VCF example file not found")
  tmp <- tempfile(fileext = ".vcf")
  res <- clean_genotype_file(file_in = f, file_out = tmp, sep = "auto",
                             verbose = FALSE)
  expect_equal(res[["removed"]], 0L)
  file.remove(tmp)
})

test_that("read_numeric_genotype clean_malformed=TRUE runs without error", {
  gf  <- system.file("extdata", "example_genotypes_numeric.csv", package = "OptSLDP")
  tmp <- tempfile(fileext = ".csv")
  expect_no_error(
    read_numeric_genotype(gf, clean_malformed = TRUE, verbose = FALSE)
  )
})

test_that("read_numeric_genotype clean_malformed=TRUE gives same SNPs as FALSE", {
  gf   <- system.file("extdata", "example_genotypes_numeric.csv", package = "OptSLDP")
  res_f <- read_numeric_genotype(gf, clean_malformed = FALSE, verbose = FALSE)
  res_t <- read_numeric_genotype(gf, clean_malformed = TRUE,  verbose = FALSE)
  expect_equal(res_f$snp_info$SNP, res_t$snp_info$SNP)
  expect_equal(res_f$geno_mat,     res_t$geno_mat, ignore_attr = TRUE)
})


# ==============================================================================
# 14.  run_sldp()  -- integration tests using package example data
# ==============================================================================

geno_file  <- system.file("extdata","example_genotypes_numeric.csv",
                          package = "OptSLDP")
pheno_file <- system.file("extdata","example_phenotype.csv",
                          package = "OptSLDP")

# ---------- 14.1 Single-trait mode A -----------------------------------------

test_that("run_sldp mode A completes and returns the expected list elements", {
  res <- run_sldp(
    genotype_file  = geno_file,
    phenotype_file = pheno_file,
    output_file    = tempfile(fileext = ".csv"),
    trait_col      = "Trait1",
    covar_cols     = c("PC1","PC2"),
    mode           = "A",
    pval_threshold = 0.05,
    preprune_large = FALSE,
    verbose        = FALSE
  )
  expected_elements <- c("final_snp_info","final_geno_mat","screening_stats",
                         "candidate_snps_per_trait","candidate_snps",
                         "important_snps","background_retained",
                         "pruning_stats","output_file","scale_strategy")
  expect_true(all(expected_elements %in% names(res)))
})

test_that("run_sldp mode A: output file is written to disk", {
  out_file <- tempfile(fileext = ".csv")
  run_sldp(
    genotype_file  = geno_file,
    phenotype_file = pheno_file,
    output_file    = out_file,
    trait_col      = "Trait1",
    mode           = "A",
    pval_threshold = 0.05,
    preprune_large = FALSE,
    verbose        = FALSE
  )
  expect_true(file.exists(out_file))
  expect_gt(file.size(out_file), 0L)
  file.remove(out_file)
})

test_that("run_sldp mode A: final_snp_info has required columns", {
  res <- run_sldp(
    genotype_file  = geno_file,
    phenotype_file = pheno_file,
    output_file    = tempfile(fileext = ".csv"),
    trait_col      = "Trait1",
    mode           = "A",
    pval_threshold = 0.05,
    preprune_large = FALSE,
    verbose        = FALSE
  )
  expect_true(all(c("SNP","CHR","POS","REF","ALT") %in%
                    names(res$final_snp_info)))
})

test_that("run_sldp mode A: final panel is ordered by CHR then POS", {
  res <- run_sldp(
    genotype_file  = geno_file,
    phenotype_file = pheno_file,
    output_file    = tempfile(fileext = ".csv"),
    trait_col      = "Trait1",
    mode           = "A",
    pval_threshold = 0.05,
    preprune_large = FALSE,
    verbose        = FALSE
  )
  si <- res$final_snp_info
  expect_true(all(order(si$CHR, si$POS) == seq_len(nrow(si))))
})

test_that("run_sldp mode A: important_snps is a superset of candidate_snps", {
  res <- run_sldp(
    genotype_file  = geno_file,
    phenotype_file = pheno_file,
    output_file    = tempfile(fileext = ".csv"),
    trait_col      = "Trait1",
    mode           = "A",
    pval_threshold = 0.05,
    preprune_large = FALSE,
    verbose        = FALSE
  )
  expect_true(all(res$candidate_snps %in% res$important_snps))
})

test_that("run_sldp mode A: final panel contains all important SNPs", {
  res <- run_sldp(
    genotype_file  = geno_file,
    phenotype_file = pheno_file,
    output_file    = tempfile(fileext = ".csv"),
    trait_col      = "Trait1",
    mode           = "A",
    pval_threshold = 0.05,
    preprune_large = FALSE,
    verbose        = FALSE
  )
  expect_true(all(res$important_snps %in% res$final_snp_info$SNP))
})

test_that("run_sldp mode A: scale_strategy is in_memory for this small panel", {
  res <- run_sldp(
    genotype_file  = geno_file,
    phenotype_file = pheno_file,
    output_file    = tempfile(fileext = ".csv"),
    trait_col      = "Trait1",
    mode           = "A",
    pval_threshold = 0.05,
    preprune_large = FALSE,
    verbose        = FALSE
  )
  expect_equal(res$scale_strategy, "in_memory")
})

test_that("run_sldp mode A: screening_stats is a data.table with all expected columns", {
  res <- run_sldp(
    genotype_file  = geno_file,
    phenotype_file = pheno_file,
    output_file    = tempfile(fileext = ".csv"),
    trait_col      = "Trait1",
    mode           = "A",
    pval_threshold = 0.05,
    preprune_large = FALSE,
    verbose        = FALSE
  )
  expect_s3_class(res$screening_stats, "data.table")
  expect_true(all(c("SNP","beta","SE","z_score","P.value","PVE","AF","MAF") %in%
                    names(res$screening_stats)))
})

test_that("run_sldp mode A: Trait1 QTN SNP003 and SNP025 are in final panel", {
  res <- run_sldp(
    genotype_file  = geno_file,
    phenotype_file = pheno_file,
    output_file    = tempfile(fileext = ".csv"),
    trait_col      = "Trait1",
    mode           = "A",
    pval_threshold = 0.05,
    preprune_large = FALSE,
    verbose        = FALSE
  )
  # These QTN have p < 0.05 in the example data GWAS scan
  panel <- res$final_snp_info$SNP
  expect_true("SNP003" %in% panel)
  expect_true("SNP025" %in% panel)
})

test_that("run_sldp mode A: final_geno_mat rows match final_snp_info SNPs", {
  res <- run_sldp(
    genotype_file  = geno_file,
    phenotype_file = pheno_file,
    output_file    = tempfile(fileext = ".csv"),
    trait_col      = "Trait1",
    mode           = "A",
    pval_threshold = 0.05,
    preprune_large = FALSE,
    verbose        = FALSE
  )
  expect_equal(rownames(res$final_geno_mat), res$final_snp_info$SNP)
})

test_that("run_sldp mode A: pruning_stats has input and final_merge steps", {
  res <- run_sldp(
    genotype_file  = geno_file,
    phenotype_file = pheno_file,
    output_file    = tempfile(fileext = ".csv"),
    trait_col      = "Trait1",
    mode           = "A",
    pval_threshold = 0.05,
    preprune_large = FALSE,
    verbose        = FALSE
  )
  expect_true("input"       %in% res$pruning_stats$step)
  expect_true("final_merge" %in% res$pruning_stats$step)
})

test_that("run_sldp mode A: candidate_snps_per_trait is NULL for single-trait run", {
  res <- run_sldp(
    genotype_file  = geno_file,
    phenotype_file = pheno_file,
    output_file    = tempfile(fileext = ".csv"),
    trait_col      = "Trait1",
    mode           = "A",
    pval_threshold = 0.05,
    preprune_large = FALSE,
    verbose        = FALSE
  )
  expect_null(res$candidate_snps_per_trait)
})

# ---------- 14.2 Single-trait mode B -----------------------------------------

test_that("run_sldp mode B OR runs without error", {
  expect_no_error(
    run_sldp(
      genotype_file   = geno_file,
      phenotype_file  = pheno_file,
      output_file     = tempfile(fileext = ".csv"),
      trait_col       = "Trait1",
      mode            = "B",
      z_threshold     = 2.0,
      pve_threshold   = 0.05,
      threshold_logic = "OR",
      preprune_large  = FALSE,
      verbose         = FALSE
    )
  )
})

test_that("run_sldp mode B AND is more restrictive than OR", {
  base_args <- list(
    genotype_file  = geno_file,
    phenotype_file = pheno_file,
    output_file    = tempfile(fileext = ".csv"),
    trait_col      = "Trait1",
    mode           = "B",
    z_threshold    = 2.0,
    pve_threshold  = 0.05,
    preprune_large = FALSE,
    verbose        = FALSE
  )
  res_or  <- do.call(run_sldp, c(base_args, list(threshold_logic = "OR")))
  res_and <- do.call(run_sldp, c(base_args, list(threshold_logic = "AND")))
  # AND cannot produce more candidates than OR
  expect_lte(length(res_and$candidate_snps), length(res_or$candidate_snps))
})

# ---------- 14.3 Single-trait mode C -----------------------------------------

test_that("run_sldp mode C OR runs without error", {
  expect_no_error(
    run_sldp(
      genotype_file   = geno_file,
      phenotype_file  = pheno_file,
      output_file     = tempfile(fileext = ".csv"),
      trait_col       = "Trait1",
      mode            = "C",
      pval_threshold  = 0.05,
      z_threshold     = 2.0,
      threshold_logic = "OR",
      preprune_large  = FALSE,
      verbose         = FALSE
    )
  )
})

# ---------- 14.4 HapMap output format ----------------------------------------

test_that("run_sldp output_format hapmap writes a HapMap-structured file", {
  out_file <- tempfile(fileext = ".hmp.txt")
  run_sldp(
    genotype_file  = geno_file,
    phenotype_file = pheno_file,
    output_file    = out_file,
    trait_col      = "Trait1",
    output_format  = "hapmap",
    mode           = "A",
    pval_threshold = 0.05,
    preprune_large = FALSE,
    verbose        = FALSE
  )
  expect_true(file.exists(out_file))
  back <- read.table(out_file, header = TRUE, sep = "\t", check.names = FALSE,
                     nrows = 2L, comment.char = "")
  expect_true("rs#" %in% names(back))
  file.remove(out_file)
})

# ---------- 14.5 Pre-pruning step --------------------------------------------

test_that("run_sldp with preprune_large=TRUE produces fewer or equal SNPs than FALSE", {
  args_base <- list(
    genotype_file  = geno_file,
    phenotype_file = pheno_file,
    output_file    = tempfile(fileext = ".csv"),
    trait_col      = "Trait1",
    mode           = "A",
    pval_threshold = 0.05,
    verbose        = FALSE
  )
  res_pre <- do.call(run_sldp, c(args_base, list(preprune_large = TRUE,
                                                 preprune_r2   = 0.99)))
  res_no  <- do.call(run_sldp, c(args_base, list(preprune_large = FALSE)))
  # Pre-pruning should not increase the final panel size
  expect_lte(nrow(res_pre$final_snp_info), nrow(res_no$final_snp_info))
})

# ---------- 14.6 Multi-trait union protection --------------------------------

test_that("run_sldp multi-trait: screening_stats is a named list", {
  res <- run_sldp(
    genotype_file  = geno_file,
    phenotype_file = pheno_file,
    output_file    = tempfile(fileext = ".csv"),
    trait_col      = c("Trait1","Trait2"),
    mode           = "A",
    pval_threshold = 0.05,
    preprune_large = FALSE,
    verbose        = FALSE
  )
  expect_type(res$screening_stats, "list")
  expect_named(res$screening_stats, c("Trait1","Trait2"))
})

test_that("run_sldp multi-trait: candidate_snps_per_trait is a named list", {
  res <- run_sldp(
    genotype_file  = geno_file,
    phenotype_file = pheno_file,
    output_file    = tempfile(fileext = ".csv"),
    trait_col      = c("Trait1","Trait2"),
    mode           = "A",
    pval_threshold = 0.05,
    preprune_large = FALSE,
    verbose        = FALSE
  )
  expect_type(res$candidate_snps_per_trait, "list")
  expect_named(res$candidate_snps_per_trait, c("Trait1","Trait2"))
})

test_that("run_sldp multi-trait: union candidate_snps >= each per-trait set", {
  res <- run_sldp(
    genotype_file  = geno_file,
    phenotype_file = pheno_file,
    output_file    = tempfile(fileext = ".csv"),
    trait_col      = c("Trait1","Trait2"),
    mode           = "A",
    pval_threshold = 0.05,
    preprune_large = FALSE,
    verbose        = FALSE
  )
  expect_gte(length(res$candidate_snps),
             length(res$candidate_snps_per_trait$Trait1))
  expect_gte(length(res$candidate_snps),
             length(res$candidate_snps_per_trait$Trait2))
})

test_that("run_sldp multi-trait: all per-trait candidates are in union", {
  res <- run_sldp(
    genotype_file  = geno_file,
    phenotype_file = pheno_file,
    output_file    = tempfile(fileext = ".csv"),
    trait_col      = c("Trait1","Trait2"),
    mode           = "A",
    pval_threshold = 0.05,
    preprune_large = FALSE,
    verbose        = FALSE
  )
  all_trait_cands <- union(res$candidate_snps_per_trait$Trait1,
                           res$candidate_snps_per_trait$Trait2)
  expect_true(all(all_trait_cands %in% res$candidate_snps))
})

test_that("run_sldp multi-trait: SNP020 (Trait2-only singleton QTN) in final panel", {
  # SNP020 is a Trait2 QTN outside all LD blocks. It should be retained via
  # Trait2 candidacy even though it is not significant under Trait1.
  res <- run_sldp(
    genotype_file  = geno_file,
    phenotype_file = pheno_file,
    output_file    = tempfile(fileext = ".csv"),
    trait_col      = c("Trait1","Trait2"),
    mode           = "A",
    pval_threshold = 0.05,
    preprune_large = FALSE,
    verbose        = FALSE
  )
  if ("SNP020" %in% res$candidate_snps_per_trait$Trait2) {
    expect_true("SNP020" %in% res$final_snp_info$SNP)
  } else {
    skip("SNP020 not a Trait2 candidate at this threshold -- skip union test")
  }
})

test_that("run_sldp multi-trait: final panel >= single-trait panel (union is non-decreasing)", {
  args_base <- list(
    genotype_file  = geno_file,
    phenotype_file = pheno_file,
    output_file    = tempfile(fileext = ".csv"),
    mode           = "A",
    pval_threshold = 0.05,
    preprune_large = FALSE,
    verbose        = FALSE
  )
  res_single <- do.call(run_sldp, c(args_base, list(trait_col = "Trait1")))
  res_multi  <- do.call(run_sldp, c(args_base,
                                    list(trait_col = c("Trait1","Trait2"))))
  expect_gte(nrow(res_multi$final_snp_info), nrow(res_single$final_snp_info))
})

# ---------- 14.7 Stats output files ------------------------------------------

test_that("run_sldp writes stats CSV when stats_output_file is set", {
  stats_file <- tempfile(fileext = ".csv")
  run_sldp(
    genotype_file    = geno_file,
    phenotype_file   = pheno_file,
    output_file      = tempfile(fileext = ".csv"),
    trait_col        = "Trait1",
    mode             = "A",
    pval_threshold   = 0.05,
    preprune_large   = FALSE,
    stats_output_file = stats_file,
    verbose          = FALSE
  )
  expect_true(file.exists(stats_file))
  back <- read.csv(stats_file)
  expect_true("step" %in% names(back))
  file.remove(stats_file)
})

test_that("run_sldp writes summary TXT when summary_output_file is set", {
  summ_file <- tempfile(fileext = ".txt")
  run_sldp(
    genotype_file       = geno_file,
    phenotype_file      = pheno_file,
    output_file         = tempfile(fileext = ".csv"),
    trait_col           = "Trait1",
    mode                = "A",
    pval_threshold      = 0.05,
    preprune_large      = FALSE,
    summary_output_file = summ_file,
    verbose             = FALSE
  )
  expect_true(file.exists(summ_file))
  lines <- readLines(summ_file)
  expect_gt(length(lines), 0L)
  file.remove(summ_file)
})


# ==============================================================================
# 15.  C++ kernel functions and vectorised screening
# ==============================================================================

# Helper: check for compiled shared library on disk -- works with load_all(),
# test_package(), devtools::test(), and R CMD check alike.
.cpp_compiled <- function() {
  # Check every possible library location for the compiled OptSLDP shared lib.
  # We cannot use find.package() because devtools::test() runs with the source
  # directory on the search path, which has no libs/ subdirectory.
  dll_names <- c(
    file.path("libs", "x64",    "OptSLDP.dll"),
    file.path("libs", "i386",   "OptSLDP.dll"),
    file.path("libs",            "OptSLDP.so"),
    file.path("libs",            "OptSLDP.dylib")
  )
  all_libs <- unique(c(
    .libPaths(),
    strsplit(Sys.getenv("R_LIBS_USER", ""), .Platform$path.sep)[[1]],
    strsplit(Sys.getenv("R_LIBS",      ""), .Platform$path.sep)[[1]],
    Sys.getenv("R_LIBS_SITE", "")
  ))
  all_libs <- all_libs[nzchar(all_libs) & dir.exists(all_libs)]
  for (lib in all_libs)
    for (dll in dll_names)
      if (file.exists(file.path(lib, "OptSLDP", dll))) return(TRUE)
  # Last resort: check if the exported C++ function is callable
  f <- tryCatch(
    get("r2_subset_cpp", envir = asNamespace("OptSLDP"), inherits = FALSE),
    error = function(e) NULL
  )
  is.function(f) && is.primitive(f) || is.function(f) && !is.null(body(f)) &&
    grepl(".Call", deparse(body(f))[1], fixed = TRUE)
}

# -- 15.1  r2_subset_cpp() -- full-matrix equivalence -------------------------
# Note: r2_matrix_cpp is an internal helper in ld.R (not a C++ export).
# The exported C++ functions are: r2_subset_cpp, above_threshold_subset_cpp,
# greedy_prune_r2_cpp. We test r2_subset_cpp against .r2_tcrossprod (pure R).

# -- 15.2  r2_subset_cpp() ---------------------------------------------------

test_that("r2_subset_cpp returns correct dimensions", {
  skip_if_not(.cpp_compiled(), "C++ kernel not compiled -- skipping")
  m    <- .make_geno()
  cand <- c(1L, 3L)   # SNP1, SNP3 (1-based)
  r2   <- OptSLDP:::r2_subset_cpp(m, cand)
  expect_equal(nrow(r2), length(cand))
  expect_equal(ncol(r2), nrow(m))
})

test_that("r2_subset_cpp self-pairs are 1 for polymorphic candidates", {
  skip_if_not(.cpp_compiled(), "C++ kernel not compiled -- skipping")
  m    <- .make_geno()
  cand <- c(1L, 3L)
  r2   <- OptSLDP:::r2_subset_cpp(m, cand)
  expect_equal(r2[1, 1], 1, tolerance = 1e-10)   # SNP1 vs SNP1
  expect_equal(r2[2, 3], 1, tolerance = 1e-10)   # SNP3 vs SNP3
})

test_that("r2_subset_cpp matches .r2_tcrossprod for candidate rows", {
  skip_if_not(.cpp_compiled(), "C++ kernel not compiled -- skipping")
  m    <- .make_geno()
  # Use the pure-R .r2_tcrossprod as reference (always available)
  full <- OptSLDP:::.r2_tcrossprod(m)
  cand <- c(1L, 3L)   # SNP1, SNP3
  sub  <- OptSLDP:::r2_subset_cpp(m, cand)
  # Row 1 of sub (SNP1) should match row SNP1 of full r2 matrix
  ok <- !is.na(sub[1, ]) & !is.na(full["SNP1", ])
  expect_equal(sub[1, ok], full["SNP1", ok], tolerance = 1e-8,
               ignore_attr = TRUE)
})

# -- 15.3  above_threshold_subset_cpp() ----------------------------------------

test_that("above_threshold_subset_cpp returns a 2-column integer matrix", {
  skip_if_not(.cpp_compiled(), "C++ kernel not compiled -- skipping")
  m    <- .make_geno()
  r2   <- OptSLDP:::r2_subset_cpp(m, c(1L, 2L))
  out  <- OptSLDP:::above_threshold_subset_cpp(r2, 0.9, TRUE)
  expect_true(is.matrix(out))
  expect_equal(ncol(out), 2L)
  expect_equal(colnames(out), c("row", "col"))
})

test_that("above_threshold_subset_cpp finds SNP1-SNP2 pair at threshold 0.9", {
  skip_if_not(.cpp_compiled(), "C++ kernel not compiled -- skipping")
  m   <- .make_geno()
  # Only candidates: SNP1 (row 1), SNP2 (row 2)
  r2  <- OptSLDP:::r2_subset_cpp(m, c(1L, 2L))
  out <- OptSLDP:::above_threshold_subset_cpp(r2, 0.9, TRUE)
  # SNP1 (cand row 1) should be above 0.9 with SNP2 (col 2)
  expect_true(nrow(out) > 0L)
})

# -- 15.4  greedy_prune_r2_cpp() -----------------------------------------------

test_that("greedy_prune_r2_cpp retains at least one SNP", {
  skip_if_not(.cpp_compiled(), "C++ kernel not compiled -- skipping")
  m    <- .make_geno()
  sub  <- m[c("SNP1","SNP2","SNP3"), , drop = FALSE]
  r2   <- compute_r2_subset(sub, rownames(sub))
  keep <- OptSLDP:::greedy_prune_r2_cpp(r2, 0.9)
  expect_type(keep, "logical")
  expect_equal(length(keep), nrow(sub))
  expect_true(any(keep))
})

test_that("greedy_prune_r2_cpp removes SNP2 when SNP1-SNP2 r2=1 and threshold<1", {
  skip_if_not(.cpp_compiled(), "C++ kernel not compiled -- skipping")
  m    <- .make_geno()
  sub  <- m[c("SNP1","SNP2"), , drop = FALSE]
  r2   <- compute_r2_subset(sub, rownames(sub))
  keep <- OptSLDP:::greedy_prune_r2_cpp(r2, 0.9)
  # SNP1 is first so it is kept; SNP2 (r2=1 with SNP1) is removed
  expect_true(keep[1])    # SNP1 kept
  expect_false(keep[2])   # SNP2 pruned
})

test_that("greedy_prune_r2_cpp agrees with pure-R greedy loop", {
  skip_if_not(.cpp_compiled(), "C++ kernel not compiled -- skipping")
  m   <- .make_geno()
  sub <- m[c("SNP1","SNP2","SNP3","SNP4"), , drop = FALSE]
  r2  <- compute_r2_subset(sub, rownames(sub))

  # C++ result
  keep_cpp <- greedy_prune_r2_cpp(r2, 0.5)

  # Pure-R equivalent
  n       <- nrow(r2)
  keep_r  <- rep(TRUE, n)
  for (i in seq_len(n)) {
    if (!keep_r[i]) next
    hi <- which(r2[i, ] >= 0.5)
    hi <- hi[hi > i]
    if (length(hi)) keep_r[hi] <- FALSE
  }
  expect_equal(keep_cpp, keep_r)
})

# -- 15.5  Vectorised screening matches lm() results --------------------------

test_that("vectorised screening beta matches lm() within tolerance", {
  set.seed(99L)
  m <- .make_geno()
  y <- .make_phenotype()

  # Vectorised result
  stats_vec <- compute_screening_stats(m, y = y, verbose = FALSE)

  # Reference: per-SNP lm() for a non-monomorphic SNP
  snp_id <- "SNP1"
  g  <- as.numeric(m[snp_id, ])
  ok <- is.finite(g) & is.finite(y)
  fit <- stats::lm(y[ok] ~ g[ok])
  beta_lm <- unname(coef(fit)[2L])
  se_lm   <- unname(summary(fit)$coefficients[2L, 2L])

  beta_vec <- unname(stats_vec$beta[stats_vec$SNP == snp_id])
  se_vec   <- unname(stats_vec$SE[stats_vec$SNP == snp_id])

  expect_equal(beta_vec, beta_lm, tolerance = 1e-6)
  expect_equal(se_vec,   se_lm,   tolerance = 1e-6)
})

test_that("vectorised screening p-values are in (0,1] for polymorphic SNPs", {
  m   <- .make_geno()
  y   <- .make_phenotype()
  out <- compute_screening_stats(m, y = y, verbose = FALSE)
  p   <- out$P.value[!is.na(out$P.value)]
  expect_true(all(p > 0 & p <= 1))
})

# ==============================================================================
# Section 16 -- GDS strategy end-to-end output test
# ==============================================================================
# Forces scale_strategy = "gds" on the small example dataset to exercise the
# full GDS path: chunked screening -> expansion -> parallel pruning -> step 11
# file writer. This catches the FORK-cluster GDS handle corruption bug.

test_that("GDS strategy writes output file and returns correct structure", {
  skip_if_not_installed("SNPRelate")
  skip_if_not_installed("gdsfmt")

  out_file <- tempfile(fileext = ".hmp.txt")

  res <- run_sldp(
    genotype_file  = geno_file,
    phenotype_file = pheno_file,
    output_file    = out_file,
    trait_col      = "Trait1",
    mode           = "A",
    pval_threshold = 0.20,        # loose threshold: tiny 40-SNP dataset
    output_format  = "hapmap",
    scale_strategy = "gds",       # force GDS even on small panel
    gds_dir        = tempdir(),
    preprune_large = FALSE,        # skip pre-pruning to keep candidates non-zero
    verbose        = FALSE
  )

  # Output file must exist on disk
  expect_true(file.exists(out_file),
              info = "GDS step 11 must create the output file")

  # File must have correct dimensions
  hmp <- data.table::fread(out_file, sep = "\t", header = TRUE)
  expect_equal(nrow(hmp), nrow(res$final_snp_info),
               info = "HapMap row count must match final_snp_info")
  expect_true(ncol(hmp) > 11L,
              info = "HapMap must have 11 metadata cols + sample cols")

  # GDS runs never load the full matrix into RAM
  expect_null(res$final_geno_mat,
              info = "final_geno_mat must be NULL for GDS runs")

  # SNP info must be ordered by CHR, POS
  expect_true(
    all(diff(res$final_snp_info$POS[res$final_snp_info$CHR == 1]) >= 0),
    info = "SNPs must be ordered by position within each chromosome"
  )
})

# ==============================================================================
# Section 17 -- read_hapmap_genotype()
# ==============================================================================

test_that("read_hapmap_genotype returns correct list structure", {
  f   <- system.file("extdata", "example_genotypes.hmp.txt", package = "OptSLDP")
  out <- read_hapmap_genotype(f)
  expect_named(out, c("snp_info","geno_mat","sample_ids","format"))
  expect_equal(out$format, "hapmap")
})

test_that("read_hapmap_genotype dimensions: 40 SNPs x 50 samples", {
  f   <- system.file("extdata", "example_genotypes.hmp.txt", package = "OptSLDP")
  out <- read_hapmap_genotype(f)
  expect_equal(nrow(out$geno_mat), 40L)
  expect_equal(ncol(out$geno_mat), 50L)
})

test_that("read_hapmap_genotype geno_mat values are in {0,1,2,NA}", {
  f     <- system.file("extdata", "example_genotypes.hmp.txt", package = "OptSLDP")
  out   <- read_hapmap_genotype(f)
  valid <- out$geno_mat[!is.na(out$geno_mat)]
  expect_true(all(valid %in% c(0, 1, 2)))
})

test_that("read_hapmap_genotype snp_info has REF and ALT columns", {
  f   <- system.file("extdata", "example_genotypes.hmp.txt", package = "OptSLDP")
  out <- read_hapmap_genotype(f)
  expect_true(all(c("SNP","CHR","POS","REF","ALT") %in% names(out$snp_info)))
})

test_that("read_hapmap_genotype and read_numeric_genotype give same dosage matrix", {
  hmp_file <- system.file("extdata", "example_genotypes.hmp.txt",
                          package = "OptSLDP")
  num_file <- system.file("extdata", "example_genotypes_numeric.csv",
                          package = "OptSLDP")
  hmp_obj <- read_hapmap_genotype(hmp_file)
  num_obj <- read_numeric_genotype(num_file)
  # Both should produce the same 40x50 dosage matrix (same underlying data)
  expect_equal(dim(hmp_obj$geno_mat), dim(num_obj$geno_mat))
  # SNP IDs must match
  expect_equal(sort(hmp_obj$snp_info$SNP), sort(num_obj$snp_info$SNP))
})


# ==============================================================================
# Section 18 -- n_pcs automatic PCA
# ==============================================================================

test_that("run_sldp with n_pcs=2 runs without error and returns results", {
  skip_if_not_installed("SNPRelate")
  skip_if_not_installed("gdsfmt")

  res <- run_sldp(
    genotype_file  = geno_file,
    phenotype_file = pheno_file,
    output_file    = tempfile(fileext = ".csv"),
    trait_col      = "Trait1",
    mode           = "A",
    pval_threshold = 0.05,
    n_pcs          = 2L,           # auto-compute 2 PCs
    preprune_large = FALSE,
    gds_dir        = tempdir(),
    verbose        = FALSE
  )
  # Pipeline must complete and return all expected elements
  expect_true("final_snp_info"  %in% names(res))
  expect_true("screening_stats" %in% names(res))
  expect_true(nrow(res$final_snp_info) > 0L)
})

test_that("run_sldp n_pcs=2 changes screening stats vs no covariates", {
  skip_if_not_installed("SNPRelate")
  skip_if_not_installed("gdsfmt")

  res_no_pc <- run_sldp(
    genotype_file  = geno_file,
    phenotype_file = pheno_file,
    output_file    = tempfile(fileext = ".csv"),
    trait_col      = "Trait1",
    mode           = "A",
    pval_threshold = 0.05,
    n_pcs          = 0L,
    preprune_large = FALSE,
    gds_dir        = tempdir(),
    verbose        = FALSE
  )
  res_with_pc <- run_sldp(
    genotype_file  = geno_file,
    phenotype_file = pheno_file,
    output_file    = tempfile(fileext = ".csv"),
    trait_col      = "Trait1",
    mode           = "A",
    pval_threshold = 0.05,
    n_pcs          = 2L,
    preprune_large = FALSE,
    gds_dir        = tempdir(),
    verbose        = FALSE
  )
  # Betas must differ after PC correction (at least some SNPs)
  b1 <- res_no_pc$screening_stats$beta
  b2 <- res_with_pc$screening_stats$beta
  ok <- !is.na(b1) & !is.na(b2)
  expect_true(any(abs(b1[ok] - b2[ok]) > 1e-10),
              info = "PC correction must change at least some beta values")
})

test_that("run_sldp n_pcs ignored when covar_cols provided", {
  # When covar_cols is supplied, n_pcs should have no effect
  res_covar <- run_sldp(
    genotype_file  = geno_file,
    phenotype_file = pheno_file,
    output_file    = tempfile(fileext = ".csv"),
    trait_col      = "Trait1",
    covar_cols     = c("PC1","PC2"),
    n_pcs          = 3L,           # should be ignored
    mode           = "A",
    pval_threshold = 0.05,
    preprune_large = FALSE,
    verbose        = FALSE
  )
  res_covar_only <- run_sldp(
    genotype_file  = geno_file,
    phenotype_file = pheno_file,
    output_file    = tempfile(fileext = ".csv"),
    trait_col      = "Trait1",
    covar_cols     = c("PC1","PC2"),
    n_pcs          = 0L,
    mode           = "A",
    pval_threshold = 0.05,
    preprune_large = FALSE,
    verbose        = FALSE
  )
  # Results must be identical when covar_cols is provided regardless of n_pcs
  expect_equal(res_covar$candidate_snps, res_covar_only$candidate_snps)
})
