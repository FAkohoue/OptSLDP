# ==============================================================================
# OptSLDP -- sldp.R
# Main pipeline: run_sldp()
#
# Changes in v4
# -------------
# Multi-trait:
#   trait_col accepts a character vector.  When length > 1 the pipeline runs
#   screening and candidate selection separately for every trait, then takes
#   the UNION of candidate SNPs before expansion and background pruning.
#   Every SNP important for any trait is protected.  The return value carries
#   a per-trait breakdown alongside the shared final panel.
#
# Inspiration 5 -- gc() in the background pruning chromosome loop:
#   After processing each chromosome in the in-memory background pruning loop
#   the r^2 matrix is explicitly freed and gc(FALSE) is called, preventing
#   allocator fragmentation from accumulating across 20-30 chromosomes.
# ==============================================================================


# -- Internal: build GDS context -----------------------------------------------

#' @keywords internal
#' @noRd
.build_gds_context <- function(snp_info, geno_mat, strategy,
                               gds_dir, n_cores, verbose) {
  ctx <- list(strategy = strategy, geno_mat = geno_mat,
              genofile = NULL, gds_path = NULL, n_cores = n_cores)

  if (identical(strategy, "gds")) {
    .assert_packages("SNPRelate")
    if (!dir.exists(gds_dir)) dir.create(gds_dir, recursive = TRUE)
    gds_path <- file.path(gds_dir, "sldp_main.gds")
    .write_gds(geno_mat, snp_info, gds_path,
               n_cores = n_cores, verbose = verbose)
    ctx$gds_path <- gds_path
    ctx$genofile <- .open_gds(gds_path, key = "main")
    ctx$geno_mat <- NULL   # free RAM immediately after writing to disk
  }
  ctx
}


# -- Main pipeline --------------------------------------------------------------

#' Run the Selective Linkage Disequilibrium Pruning (SLDP) pipeline
#'
#' End-to-end workflow.  The user supplies exactly **two input files** -- a
#' genotype file and a phenotype file -- plus an output path.  All internal
#' steps (phenotype reading, sample alignment, GWAS screening, candidate
#' selection, important-SNP expansion, background pruning, output writing) are
#' handled automatically.
#'
#' @section Multi-trait analysis:
#' Set `trait_col` to a character vector containing more than one column name
#' from the phenotype file.  The pipeline then:
#' \enumerate{
#'   \item Runs screening (marginal regression + residualisation) **separately**
#'     for each trait.
#'   \item Selects candidate SNPs independently for each trait using the shared
#'     threshold settings.
#'   \item Takes the **union** of all per-trait candidate sets before expansion.
#'   \item Expands the union into the protected set (positional window + LD
#'     neighbours) and runs background pruning once on the remainder.
#' }
#' The final pruned panel therefore retains every SNP that is important for
#' **any** of the requested traits.  The return value includes per-trait
#' screening statistics and per-trait candidate lists alongside the single
#' shared output file.
#'
#' @section Scale strategy:
#' | Strategy    | Trigger (default)      | LD backend                          |
#' |-------------|------------------------|-------------------------------------|
#' | `in_memory` | n_snps <= 200 000       | `cor()` on full matrix in RAM       |
#' | `chunked`   | 200K < n_snps <= 2M    | `cor()` per chromosome in RAM       |
#' | `gds`       | n_snps > 2M            | `snpgdsLDMat()` from disk           |
#'
#' Adjustable via `options(optsldp.thresh_small, optsldp.thresh_medium)`.
#'
#' @param genotype_file       Path to genotype file (numeric CSV, HapMap, VCF).
#' @param phenotype_file      Path to phenotype file.  Must contain a sample ID
#'   column (`sample_col`) and one or more numeric trait columns (`trait_col`).
#' @param output_file         Path for the pruned output genotype file.
#' @param sample_col          Sample ID column in the phenotype file.
#'   Default `"Sample"`.
#' @param trait_col           Trait column name(s).  A single string or a
#'   character vector for multi-trait analysis.  Default `"Trait1"`.
#' @param covar_cols          Covariate column name(s) in the phenotype file.
#'   Covariates are shared across all traits.  `NULL` = no covariates.
#' @param n_pcs               Number of principal components to compute
#'   automatically and use as covariates. Computed **after MAF filtering** from
#'   a chromosome-balanced SNP subset via GRM eigendecomposition -- no LD
#'   pruning pass required. Only used when `covar_cols` is `NULL`. Set to `0`
#'   (default) to disable. Typical values: 3-5 for structured populations.
#' @param pca_method          Eigendecomposition backend for automatic PCA:
#'   \itemize{
#'     \item `"auto"` (default): use `RSpectra::eigs_sym()` if the `RSpectra`
#'       package is installed, otherwise fall back to base `eigen()`. A message
#'       is printed in either case when `verbose = TRUE`.
#'     \item `"grm_eigen"`: always use base `eigen()`. Computes all
#'       eigenvalues of the GRM; reliable on all platforms with no extra
#'       dependencies.
#'     \item `"rspectra"`: always use `RSpectra::eigs_sym()`. Computes only
#'       the top `n_pcs` eigenvalues; faster for large sample sets (> 500
#'       samples). Requires `RSpectra` to be installed:
#'       `install.packages("RSpectra")`. Raises an error if not available.
#'   }
#' @param pca_max_snps        Maximum SNPs used for PCA when the post-MAF SNP
#'   count exceeds 1 M. Default `40000L`. Smaller subsets (20k for 5k-200k
#'   SNPs, 30k for 200k-1M SNPs) are chosen automatically.
#' @param pca_seed            Random seed for chromosome-balanced SNP sampling.
#'   Default `1L`.
#' @param format              Genotype format: `"auto"` (default), `"numeric"`,
#'   `"hapmap"`, or `"vcf"`.
#' @param output_format       Output format: `"numeric"` (default) or
#'   `"hapmap"`.
#' @param mode                Candidate selection mode: `"A"` (p-value),
#'   `"B"` (effect-size), or `"C"` (hybrid).
#' @param scale_strategy      Override auto-detection: `"in_memory"`,
#'   `"chunked"`, or `"gds"`.  `NULL` = auto.
#' @param gds_dir             Directory for intermediate GDS files.
#'   Default `tempdir()`.
#' @param n_cores             Threads for GDS operations.
#'   Default `max(1, detectCores() - 1)`.
#' @param maf_threshold       Minimum MAF. Default `0.05`.
#' @param preprune_large      Apply high-LD pre-pruning. Default `TRUE`.
#' @param preprune_r2         r^2 threshold for pre-pruning. Default `0.99`.
#' @param pval_threshold      P-value upper bound (modes A / C).
#' @param z_threshold         |z-score| lower bound (modes B / C).
#' @param pve_threshold       PVE lower bound (modes B / C).
#' @param threshold_logic     `"AND"` (default) or `"OR"`.
#' @param window_kb           Positional expansion window in kb. Default `50`.
#' @param include_ld_neighbors Include LD neighbours. Default `TRUE`.
#' @param r2_flag             r^2 threshold for neighbour inclusion. Default
#'   `0.90`.
#' @param r2_genome           r^2 threshold for background pruning. Default
#'   `0.80`.
#' @param slide_max_bp        Sliding window (bp) for GDS pruning. Default
#'   `1 000 000`.
#' @param stats_output_file   Optional path for pruning statistics CSV.
#' @param summary_output_file Optional path for plain-text summary.
#' @param clean_malformed     If `TRUE`, stream-clean the genotype file
#'   before reading by removing any lines whose column count does not
#'   match the header.  Works for all accepted formats (numeric CSV,
#'   HapMap, VCF).  Needed for files from NGSEP and other callers that
#'   produce malformed lines.  Adds one extra streaming pass.
#'   Default `FALSE`.
#' @param verbose             Print timestamped progress. Default `TRUE`.
#'
#' @return A named list:
#'   \describe{
#'     \item{`final_snp_info`}{`data.table` of retained SNP metadata.}
#'     \item{`final_geno_mat`}{Retained genotype matrix (`NULL` for GDS runs).}
#'     \item{`screening_stats`}{Single `data.table` (one trait) **or** named
#'       list of `data.table`s (multiple traits).}
#'     \item{`candidate_snps_per_trait`}{Named list of per-trait candidate SNP
#'       IDs (only present in multi-trait runs; `NULL` otherwise).}
#'     \item{`candidate_snps`}{Union of candidate SNPs across all traits (or
#'       the single-trait candidate set).}
#'     \item{`important_snps`}{All protected SNP IDs.}
#'     \item{`background_retained`}{Retained background SNP IDs.}
#'     \item{`pruning_stats`}{Step-by-step count `data.table`.}
#'     \item{`output_file`}{Path to the written output.}
#'     \item{`scale_strategy`}{Strategy used.}
#'   }
#'
#' @export
#' @examples
#' \dontrun{
#' geno_file  <- system.file("extdata", "example_genotypes_numeric.csv",
#'                            package = "OptSLDP")
#' pheno_file <- system.file("extdata", "example_phenotype.csv",
#'                            package = "OptSLDP")
#'
#' # -- Single trait, numeric output (default) ----------------------------------
#' res <- run_sldp(
#'   genotype_file  = geno_file,
#'   phenotype_file = pheno_file,
#'   output_file    = tempfile(fileext = ".csv"),
#'   trait_col      = "Trait1",
#'   mode           = "A",
#'   pval_threshold = 0.05,
#'   output_format  = "numeric"    # default: values are 0/1/2/NA
#' )
#'
#' # -- Single trait, HapMap output ---------------------------------------------
#' res_hmp <- run_sldp(
#'   genotype_file  = geno_file,
#'   phenotype_file = pheno_file,
#'   output_file    = tempfile(fileext = ".hmp.txt"),  # use .hmp.txt extension
#'   trait_col      = "Trait1",
#'   mode           = "A",
#'   pval_threshold = 0.05,
#'   output_format  = "hapmap"     # nucleotide calls: AA/AT/TT/NN
#' )
#'
#' # -- Multi-trait with automatic PCA (population structure correction) -------
#' res <- run_sldp(
#'   genotype_file  = geno_file,
#'   phenotype_file = pheno_file,
#'   output_file    = tempfile(fileext = ".csv"),
#'   trait_col      = c("Trait1", "Trait2"),
#'   n_pcs          = 3L,             # auto-compute 3 PCs from genotypes
#'   mode           = "A",
#'   pval_threshold = 0.05,
#'   output_format  = "numeric"
#' )
#' # Inspect per-trait results
#' names(res$screening_stats)          # "Trait1" "Trait2"
#' res$candidate_snps_per_trait$Trait1
#'
#' # -- Multi-trait with user-supplied PCs ------------------------------------
#' res <- run_sldp(
#'   genotype_file  = geno_file,
#'   phenotype_file = pheno_file,
#'   output_file    = tempfile(fileext = ".csv"),
#'   trait_col      = c("Trait1", "Trait2"),
#'   covar_cols     = c("PC1", "PC2"),  # columns already in phenotype file
#'   mode           = "A",
#'   pval_threshold = 0.05,
#'   output_format  = "numeric"
#' )
#' names(res$screening_stats)          # "Trait1" "Trait2"
#' res$candidate_snps_per_trait$Trait1
#' }
run_sldp <- function(genotype_file,
                     phenotype_file,
                     output_file,
                     sample_col             = "Sample",
                     trait_col              = "Trait1",
                     covar_cols             = NULL,
                     n_pcs                  = 0L,
                     pca_method             = c("auto", "grm_eigen", "rspectra"),
                     pca_max_snps           = 40000L,
                     pca_seed               = 1L,
                     format                 = c("auto", "numeric",
                                                "hapmap", "vcf"),
                     output_format          = c("numeric", "hapmap"),
                     mode                   = c("A", "B", "C"),
                     scale_strategy         = NULL,
                     gds_dir                = tempdir(),
                     n_cores                = max(1L,
                                                  parallel::detectCores() - 1L),
                     maf_threshold          = 0.05,
                     preprune_large         = TRUE,
                     preprune_r2            = 0.99,
                     pval_threshold         = NULL,
                     z_threshold            = NULL,
                     pve_threshold          = NULL,
                     threshold_logic        = c("AND", "OR"),
                     window_kb              = 50,
                     include_ld_neighbors   = TRUE,
                     r2_flag                = 0.90,
                     r2_genome              = 0.80,
                     slide_max_bp           = 1000000L,
                     stats_output_file      = NULL,
                     summary_output_file    = NULL,
                     clean_malformed        = FALSE,
                     verbose                = TRUE) {

  format          <- match.arg(format)
  output_format   <- match.arg(output_format)
  mode            <- match.arg(mode)
  threshold_logic <- match.arg(threshold_logic)
  pca_method      <- match.arg(pca_method)
  trait_col       <- as.character(trait_col)
  multi_trait     <- length(trait_col) > 1L

  pruning_stats <- .new_pruning_stats()
  on.exit(.close_all_gds_handles(), add = TRUE)
  .report_progress("SLDP started", verbose = verbose)

  # -- Step 1: Read genotypes --------------------------------------------------
  .report_progress("[1] Reading genotype data ...", verbose = verbose)
  fmt_resolved <- if (identical(format, "auto")) {
    lower <- tolower(genotype_file)
    if (grepl("\\.vcf(\\.gz)?$", lower)) "vcf"
    else if (grepl("\\.hmp(\\.txt)?$", lower) ||
             grepl("hapmap", lower)) "hapmap"
    else "numeric"
  } else format

  g        <- if (identical(fmt_resolved, "vcf")) {
    read_vcf_genotype(
      file              = genotype_file,
      clean_malformed   = clean_malformed,
      gds_dir           = gds_dir,
      n_cores           = n_cores,
      verbose           = verbose
    )
  } else {
    if (identical(fmt_resolved, "numeric"))
      read_numeric_genotype(genotype_file,
                            clean_malformed = clean_malformed,
                            verbose = verbose)
    else
      read_hapmap_genotype(genotype_file,
                           clean_malformed = clean_malformed,
                           verbose = verbose)
  }
  data.table::setDT(g$snp_info)
  snp_info <- g$snp_info
  geno_mat <- g$geno_mat
  n0       <- nrow(snp_info)
  pruning_stats <- .append_pruning_stat(pruning_stats, "input", n0, n0,
                                        paste0("format=", g$format))
  .report_progress("    Loaded ", n0, " SNPs x ", ncol(geno_mat), " samples",
                   verbose = verbose)

  # -- Step 2: Read and align phenotype ---------------------------------------
  # read_phenotype() handles sample reordering; .harmonise_samples() is the
  # hard safety check that the two ID vectors are now identical.
  .report_progress("[2] Reading phenotype data ...", verbose = verbose)
  pheno_data <- read_phenotype(
    file         = phenotype_file,
    sample_col   = sample_col,
    trait_col    = trait_col,
    covar_cols   = covar_cols,
    sample_order = g$sample_ids
  )
  .harmonise_samples(g$sample_ids, pheno_data$sample_ids)

  phenotype  <- pheno_data$phenotype   # vector (single) or matrix (multi)
  covariates <- pheno_data$covariates

  if (multi_trait) {
    .report_progress("    ", ncol(phenotype), " traits: ",
                     paste(trait_col, collapse = ", "), verbose = verbose)
  } else {
    .report_progress("    Trait: ", trait_col,
                     if (!is.null(covariates))
                       paste0("; ", ncol(covariates), " covariate(s)")
                     else "",
                     verbose = verbose)
  }

  # -- Step 3: Scale strategy & GDS context -----------------------------------
  strategy <- .select_strategy(n0, user_strategy = scale_strategy)
  .report_progress("[3] Scale strategy: ", strategy, verbose = verbose)
  ctx      <- .build_gds_context(snp_info, geno_mat, strategy,
                                 gds_dir = gds_dir, n_cores = n_cores,
                                 verbose = verbose)
  geno_mat <- ctx$geno_mat   # NULL for GDS path

  # -- Step 4: MAF filter -----------------------------------------------------
  .report_progress("[4] MAF filtering (>= ", maf_threshold, ") ...",
                   verbose = verbose)
  n_before <- n0
  maf_res  <- filter_snps_by_maf(snp_info, geno_mat,
                                 maf_threshold = maf_threshold, ctx = ctx)
  data.table::setDT(maf_res$snp_info)
  snp_info     <- maf_res$snp_info
  geno_mat     <- maf_res$geno_mat
  ctx$geno_mat <- geno_mat
  pruning_stats <- .append_pruning_stat(pruning_stats, "maf_filter",
                                        n_before, nrow(snp_info),
                                        paste0("maf>=", maf_threshold))
  .report_progress("    After MAF filter: ", nrow(snp_info), " SNPs",
                   verbose = verbose)

  # -- Step 4b: Fast automatic PCA for population structure correction --------
  # Runs AFTER MAF filtering. Uses a chromosome-balanced random sample of
  # post-MAF SNPs to build a sample x sample GRM, then extracts the top PCs
  # by eigendecomposition. No LD pruning pass required -- just one small matrix
  # extraction from GDS. Strategy:
  #
  #   post-MAF SNPs   |  PCA SNP subset
  #   ----------------|-----------------
  #   < 5 000         |  all SNPs
  #   5 000 - 200 000 |  20 000
  #   200 000 - 1 M   |  30 000
  #   > 1 M           |  40 000  (default pca_max_snps)
  #
  # Backend:
  #   "auto"       -> RSpectra::eigs_sym() if installed, else base eigen()
  #   "rspectra"   -> always RSpectra
  #   "grm_eigen"  -> always base eigen()
  if (n_pcs > 0L && is.null(covariates)) {

    .report_progress("[4b] Computing ", n_pcs,
                     " PCs for population structure (fast GRM method) ...",
                     verbose = verbose)

    total_snps <- nrow(snp_info)

    # -- Choose subset size based on post-MAF SNP count ----------------------
    target_snps <- if (total_snps < 5000L) {
      total_snps
    } else if (total_snps < 200000L) {
      min(20000L, total_snps)
    } else if (total_snps < 1000000L) {
      min(30000L, total_snps)
    } else {
      min(pca_max_snps, total_snps)
    }

    .report_progress("    SNP subset target: ", target_snps,
                     " (from ", total_snps, " post-MAF SNPs)",
                     verbose = verbose)

    # -- Chromosome-balanced sampling ----------------------------------------
    set.seed(pca_seed)
    chr_levels_pca <- unique(snp_info$CHR)
    pca_snps <- unique(unlist(lapply(chr_levels_pca, function(chr) {
      chr_ids <- snp_info$SNP[snp_info$CHR == chr]
      n_take  <- max(1L, round(target_snps * length(chr_ids) / total_snps))
      sample(chr_ids, min(n_take, length(chr_ids)), replace = FALSE)
    }), use.names = FALSE))

    .report_progress("    Chromosome-balanced subset: ", length(pca_snps),
                     " SNPs across ", length(chr_levels_pca),
                     " chromosomes", verbose = verbose)

    # -- Extract genotype subset from GDS or in-memory matrix ----------------
    if (identical(strategy, "gds")) {
      pca_X <- tryCatch(
        .extract_geno_gds(ctx$genofile,
                          snp_ids    = pca_snps,
                          sample_ids = g$sample_ids),
        error = function(e)
          stop("PCA genotype extraction failed: ", e$message, call. = FALSE)
      )
    } else {
      # In-memory / chunked: slice directly from geno_mat
      pca_X <- geno_mat[pca_snps, , drop = FALSE]
    }

    # -- Mean imputation of missing genotypes --------------------------------
    if (anyNA(pca_X)) {
      snp_means_pca <- rowMeans(pca_X, na.rm = TRUE)
      na_idx_pca    <- which(is.na(pca_X), arr.ind = TRUE)
      pca_X[na_idx_pca] <- snp_means_pca[na_idx_pca[, 1L]]
    }

    # -- Build sample x sample GRM -------------------------------------------
    # Xc = centred genotype matrix (samples x SNPs); GRM = Xc %*% t(Xc) / m
    .report_progress("    Building GRM (", ncol(pca_X),
                     " samples x ", nrow(pca_X), " SNPs) ...",
                     verbose = verbose)
    Xc_pca <- scale(t(pca_X), center = TRUE, scale = FALSE)
    rm(pca_X)
    m_pca  <- ncol(Xc_pca)
    K_pca  <- tcrossprod(Xc_pca) / m_pca
    rm(Xc_pca)

    # -- Eigendecomposition backend ------------------------------------------
    # "auto"       -> use RSpectra if installed (faster for large GRMs),
    #                 otherwise fall back to base eigen()
    # "rspectra"   -> require RSpectra; error if not installed
    # "grm_eigen"  -> always use base eigen() regardless of RSpectra availability
    #
    # Install RSpectra for faster PCA on large sample sets (> 500 samples):
    #   install.packages("RSpectra")
    backend_pca <- if (identical(pca_method, "auto")) {
      if (requireNamespace("RSpectra", quietly = TRUE)) {
        .report_progress(
          "    RSpectra available -- using fast partial eigendecomposition",
          verbose = verbose)
        "rspectra"
      } else {
        .report_progress(
          "    RSpectra not installed -- using base eigen(). ",
          "Install RSpectra for faster PCA: install.packages(\"RSpectra\")",
          verbose = verbose)
        "grm_eigen"
      }
    } else if (identical(pca_method, "rspectra")) {
      if (!requireNamespace("RSpectra", quietly = TRUE))
        stop('pca_method = "rspectra" requires the RSpectra package. ',
             'Install it with: install.packages("RSpectra")', call. = FALSE)
      "rspectra"
    } else {
      "grm_eigen"
    }
    .report_progress("    Eigendecomposition backend: ", backend_pca,
                     verbose = verbose)

    pcs <- if (identical(backend_pca, "rspectra")) {
      tryCatch(
        RSpectra::eigs_sym(K_pca, k = n_pcs)$vectors,
        error = function(e) {
          warning("RSpectra::eigs_sym() failed, falling back to base eigen(): ",
                  e$message, call. = FALSE)
          eigen(K_pca, symmetric = TRUE)$vectors[, seq_len(n_pcs), drop = FALSE]
        }
      )
    } else {
      # base eigen() -- computes all eigenvalues but reliable on any platform
      eigen(K_pca, symmetric = TRUE)$vectors[, seq_len(n_pcs), drop = FALSE]
    }
    rm(K_pca)

    # -- Attach PCs as covariates --------------------------------------------
    pc_mat <- pcs
    rownames(pc_mat) <- g$sample_ids
    colnames(pc_mat) <- paste0("PC", seq_len(n_pcs))
    covariates <- as.data.frame(pc_mat)

    .report_progress("    PC1-PC", n_pcs,
                     " added as covariates for screening",
                     verbose = verbose)
  }

  # -- Step 5 (optional): High-LD pre-pruning ---------------------------------
  if (isTRUE(preprune_large)) {
    .report_progress("[5] High-LD pre-pruning (r2 >= ", preprune_r2, ") ...",
                     verbose = verbose)
    n_before <- nrow(snp_info)
    pre      <- preprune_high_ld(snp_info, geno_mat,
                                 r2_pre = preprune_r2,
                                 slide_max_bp = slide_max_bp,
                                 ctx = ctx, verbose = verbose)
    data.table::setDT(pre$snp_info)
    snp_info     <- pre$snp_info
    geno_mat     <- pre$geno_mat
    ctx$geno_mat <- geno_mat
    pruning_stats <- .append_pruning_stat(pruning_stats, "preprune_high_ld",
                                          n_before, nrow(snp_info),
                                          paste0("r2>=", preprune_r2))
    .report_progress("    After pre-pruning: ", nrow(snp_info), " SNPs",
                     verbose = verbose)
  }

  # -- Step 6: Screening statistics -------------------------------------------
  # GDS strategy: stream one chromosome at a time, screen it, accumulate
  # only the statistics (not the genotypes). Peak RAM = one chromosome's
  # genotype matrix (~200-300 MB for 2M SNPs x 204 samples) instead of
  # the full 2.65M x 204 matrix (~4 GB).
  # In-memory/chunked: geno_mat is already in RAM, screen directly.
  .report_progress("[6] Computing screening statistics ...", verbose = verbose)

  if (identical(strategy, "gds")) {
    .report_progress("    Chunked chromosome screening (",
                     length(unique(snp_info$CHR)), " chromosomes, ",
                     "chunk_size=50000) ...", verbose = verbose)

    chr_levels <- unique(snp_info$CHR)
    chunk_size <- 50000L   # SNPs per chunk -- balances cache locality vs GDS overhead

    # Residualise phenotype on covariates ONCE before all chromosome/chunk loops
    y_list <- if (!is.null(phenotype) && is.matrix(phenotype)) {
      lapply(seq_len(ncol(phenotype)), function(j) {
        yj <- as.numeric(phenotype[, j])
        if (!is.null(covariates))
          yj <- stats::residuals(stats::lm(yj ~ .,
                                           data = as.data.frame(covariates)))
        yj
      })
    } else {
      y_single <- as.numeric(phenotype)
      if (!is.null(covariates))
        y_single <- stats::residuals(stats::lm(y_single ~ .,
                                               data = as.data.frame(covariates)))
      list(y_single)
    }
    n_traits_screen <- length(y_list)

    all_stats <- vector("list", n_traits_screen)
    for (j in seq_len(n_traits_screen)) all_stats[[j]] <- vector("list", 0L)

    for (chr_i in chr_levels) {
      chr_snp_ids <- snp_info$SNP[snp_info$CHR == chr_i]
      n_chr       <- length(chr_snp_ids)

      # Split chromosome into chunks for better cache locality
      chunk_starts <- seq(1L, n_chr, by = chunk_size)

      for (cs in chunk_starts) {
        ce        <- min(cs + chunk_size - 1L, n_chr)
        chunk_ids <- chr_snp_ids[cs:ce]

        geno_chunk <- tryCatch(
          .extract_geno_gds(ctx$genofile,
                            snp_ids    = chunk_ids,
                            sample_ids = g$sample_ids),
          error = function(e) NULL
        )
        if (is.null(geno_chunk) || nrow(geno_chunk) == 0L) next

        # Compute genotype variance ONCE per chunk -- reused across all traits
        # Genotype variance is phenotype-independent; no need to recompute per trait
        g_var_cache <- {
          row_means_c <- rowMeans(geno_chunk, na.rm = TRUE)
          g_imp_c     <- geno_chunk
          na_idx_c    <- which(is.na(geno_chunk), arr.ind = TRUE)
          if (nrow(na_idx_c) > 0L)
            g_imp_c[na_idx_c] <- row_means_c[na_idx_c[, 1L]]
          g_dev_c <- sweep(g_imp_c, 1L, row_means_c, `-`)
          rowSums(g_dev_c^2) / (ncol(geno_chunk) - 1L)
        }

        # Screen chunk for each trait, reusing g_var_cache
        for (j in seq_len(n_traits_screen)) {
          chunk_stat <- .screen_single_trait(
            geno_chunk,
            y           = y_list[[j]],
            covar       = NULL,          # already residualised above
            g_var_cache = g_var_cache,   # avoid recomputing variance per trait
            verbose     = FALSE
          )
          all_stats[[j]] <- c(all_stats[[j]], list(chunk_stat))
        }

        rm(geno_chunk, g_var_cache); gc(FALSE)
      }
    }

    # Combine chunk statistics into per-trait data.tables
    geno_screen <- NULL
    if (n_traits_screen == 1L) {
      screen_stats_pre <- data.table::rbindlist(all_stats[[1L]])
    } else {
      trait_names_screen <- if (is.matrix(phenotype)) colnames(phenotype) else trait_col
      screen_stats_pre <- stats::setNames(
        lapply(all_stats, data.table::rbindlist),
        trait_names_screen
      )
    }
    rm(all_stats); gc(FALSE)

  } else {
    geno_screen      <- geno_mat
    screen_stats_pre <- NULL
  }

  # -- Step 7: Candidate selection (single-trait or multi-trait) --------------
  if (multi_trait) {
    # -- Multi-trait branch ----------------------------------------------------
    .report_progress("[7] Multi-trait screening (",
                     length(trait_col), " traits) ...", verbose = verbose)

    if (!is.null(screen_stats_pre)) {
      # GDS path: stats already computed per chromosome in step 6
      .report_progress("  Using pre-computed chromosome screening stats ...",
                       verbose = verbose)
      screening_stats <- screen_stats_pre
      candidate_snps_per_trait <- lapply(trait_col, function(tn) {
        select_candidate_snps(
          screening_stats[[tn]],
          mode            = mode,
          pval_threshold  = pval_threshold,
          z_threshold     = z_threshold,
          pve_threshold   = pve_threshold,
          logic           = threshold_logic
        )
      })
      names(candidate_snps_per_trait) <- trait_col
      candidate_snps <- unique(unlist(candidate_snps_per_trait))
    } else {
      # In-memory/chunked path: screen now
      screen_res <- .screen_all_traits(
        geno_mat        = geno_screen,
        pheno_matrix    = phenotype,
        covar           = covariates,
        mode            = mode,
        pval_threshold  = pval_threshold,
        z_threshold     = z_threshold,
        pve_threshold   = pve_threshold,
        threshold_logic = threshold_logic,
        n_cores         = min(length(trait_col), n_cores),
        verbose         = verbose
      )
      screening_stats          <- screen_res$screening_stats
      candidate_snps_per_trait <- screen_res$candidate_snps_per_trait
      candidate_snps           <- screen_res$candidate_snps_union
    }

    # Log per-trait counts and union total
    for (tn in trait_col) {
      pruning_stats <- .append_pruning_stat(
        pruning_stats,
        paste0("candidates_", tn),
        nrow(snp_info), length(candidate_snps_per_trait[[tn]]),
        paste0("mode=", mode, "; trait=", tn)
      )
    }
    pruning_stats <- .append_pruning_stat(
      pruning_stats, "candidate_union",
      nrow(snp_info), length(candidate_snps),
      paste0("union across ", length(trait_col), " traits")
    )
    .report_progress("    Union candidate SNPs: ", length(candidate_snps),
                     verbose = verbose)

  } else {
    # -- Single-trait branch ---------------------------------------------------
    .report_progress("[7] Computing screening statistics + selecting candidates",
                     " (mode=", mode, ") ...", verbose = verbose)

    if (!is.null(screen_stats_pre)) {
      # GDS path: stats pre-computed per chromosome in step 6
      screening_stats <- screen_stats_pre
    } else {
      # In-memory/chunked path: screen now
      screening_stats <- .screen_single_trait(
        geno_screen, y = phenotype, covar = covariates, verbose = verbose
      )
    }
    candidate_snps <- select_candidate_snps(
      screening_stats,
      mode            = mode,
      pval_threshold  = pval_threshold,
      z_threshold     = z_threshold,
      pve_threshold   = pve_threshold,
      logic           = threshold_logic
    )
    candidate_snps_per_trait <- NULL   # not applicable for single trait

    pruning_stats <- .append_pruning_stat(
      pruning_stats, "candidate_threshold",
      nrow(snp_info), length(candidate_snps),
      paste0("mode=", mode)
    )
    .report_progress("    Candidate SNPs: ", length(candidate_snps),
                     verbose = verbose)
  }

  # -- Step 8: Important SNP expansion ----------------------------------------
  # Works on the union of candidate SNPs -- identical code path for both
  # single-trait and multi-trait.
  .report_progress("[8] Expanding important SNPs (+/-", window_kb, " kb) ...",
                   verbose = verbose)
  important_snps <- expand_important_snps(
    candidate_snps       = candidate_snps,
    snp_info             = snp_info,
    geno_mat             = geno_mat,
    window_kb            = window_kb,
    include_ld_neighbors = include_ld_neighbors,
    r2_flag              = r2_flag,
    ctx                  = ctx,
    verbose              = verbose
  )
  pruning_stats <- .append_pruning_stat(
    pruning_stats, "important_expansion",
    length(candidate_snps), length(important_snps),
    paste0("window_kb=", window_kb, "; r2_flag=", r2_flag)
  )
  .report_progress("    Important SNPs: ", length(important_snps),
                   verbose = verbose)

  # -- Step 9: Background pruning ---------------------------------------------
  .report_progress("[9] Pruning background SNPs (r2 >= ", r2_genome, ") ...",
                   verbose = verbose)
  remaining_snps      <- setdiff(snp_info$SNP, important_snps)
  background_retained <- prune_background_snps(
    remaining_snps = remaining_snps,
    snp_info       = snp_info,
    geno_mat       = geno_mat,
    r2_genome      = r2_genome,
    slide_max_bp   = slide_max_bp,
    ctx            = ctx,
    verbose        = verbose
  )
  pruning_stats <- .append_pruning_stat(
    pruning_stats, "background_prune",
    length(remaining_snps), length(background_retained),
    paste0("r2_genome=", r2_genome)
  )

  # -- Step 10: Merge and order ------------------------------------------------
  final_snps <- unique(c(important_snps, background_retained))
  pruning_stats <- .append_pruning_stat(
    pruning_stats, "final_merge",
    nrow(snp_info), length(final_snps),
    "important + pruned background"
  )
  .report_progress("[10] Final SNP panel: ", length(final_snps), " SNPs",
                   verbose = verbose)

  final_snp_info <- snp_info[match(final_snps, snp_info$SNP), ]
  data.table::setDT(final_snp_info)
  data.table::setorder(final_snp_info, CHR, POS)

  # -- Step 11: Write output ---------------------------------------------------
  # GDS strategy: write directly from GDS chromosome by chromosome -- never
  # load the full final panel into RAM. final_geno_mat is returned as NULL.
  # In-memory/chunked: slice the already-loaded matrix (cheap).
  #
  # IMPORTANT: reopen GDS before writing. The FORK cluster in step 9
  # (parallel background pruning) shares file descriptors with the parent
  # process. After the cluster closes, the parent's GDS handle may be in
  # an invalid state. Reopening guarantees a clean handle for step 11.
  .report_progress("[11] Writing output to: ", output_file, verbose = verbose)

  if (identical(strategy, "gds")) {
    # Reopen GDS with a fresh handle -- safe after FORK workers have exited
    .close_gds(ctx$genofile, key = "main")
    ctx$genofile <- .open_gds(ctx$gds_path, key = "main_write")

    chr_levels <- unique(final_snp_info$CHR)
    first_chr  <- TRUE
    for (chr_i in chr_levels) {
      chr_snp_ids <- final_snp_info$SNP[final_snp_info$CHR == chr_i]
      chr_si      <- final_snp_info[final_snp_info$CHR == chr_i, ]
      geno_chr    <- tryCatch(
        .extract_geno_gds(ctx$genofile,
                          snp_ids    = chr_snp_ids,
                          sample_ids = g$sample_ids),
        error = function(e) {
          warning("Step 11: failed to extract chromosome ", chr_i,
                  ": ", e$message, call. = FALSE)
          NULL
        }
      )
      if (is.null(geno_chr)) {
        warning("Step 11: skipped chromosome ", chr_i,
                " -- no genotypes extracted", call. = FALSE)
        next
      }

      if (identical(output_format, "numeric")) {
        # For first chromosome write header + rows; subsequent append rows only
        chr_df <- data.frame(
          SNP = chr_si$SNP, CHR = chr_si$CHR, POS = chr_si$POS,
          REF = chr_si$REF, ALT = chr_si$ALT,
          as.data.frame(geno_chr, check.names = FALSE),
          check.names = FALSE, stringsAsFactors = FALSE
        )
        data.table::fwrite(chr_df, file = output_file,
                           append = !first_chr, col.names = first_chr,
                           sep = ",", quote = FALSE)
      } else {
        # HapMap: convert dosage to nucleotide calls then append
        geno_hmp <- matrix(NA_character_, nrow = nrow(geno_chr),
                           ncol = ncol(geno_chr),
                           dimnames = dimnames(geno_chr))
        for (i in seq_len(nrow(geno_chr))) {
          r  <- chr_si$REF[i]; a <- chr_si$ALT[i]; gi <- geno_chr[i, ]
          geno_hmp[i, gi == 0]   <- paste0(r, r)
          geno_hmp[i, gi == 1]   <- paste0(r, a)
          geno_hmp[i, gi == 2]   <- paste0(a, a)
          geno_hmp[i, is.na(gi)] <- "NN"
        }
        hmp_df <- data.frame(
          "rs#"        = chr_si$SNP,
          alleles      = paste0(chr_si$REF, "/", chr_si$ALT),
          chrom        = chr_si$CHR,   pos = chr_si$POS,
          strand       = "+",          "assembly#" = NA_character_,
          center       = NA_character_, protLSID   = NA_character_,
          assayLSID    = NA_character_, panelLSID  = NA_character_,
          QCcode       = NA_character_,
          as.data.frame(geno_hmp, check.names = FALSE),
          check.names = FALSE, stringsAsFactors = FALSE
        )
        data.table::fwrite(hmp_df, file = output_file,
                           append = !first_chr, col.names = first_chr,
                           sep = "	", quote = FALSE, na = "NN")
      }
      first_chr <- FALSE
      rm(geno_chr); gc(FALSE)
    }
    final_geno_mat <- NULL   # never loaded into RAM; use extract_final_geno()

  } else {
    # In-memory / chunked: matrix already in RAM, slice and write at once
    final_geno_mat <- geno_mat[final_snp_info$SNP, , drop = FALSE]
    if (identical(output_format, "numeric"))
      write_numeric_genotype(final_snp_info, final_geno_mat, output_file)
    else
      write_hapmap_genotype(final_snp_info, final_geno_mat, output_file)
  }

  # -- Step 12: Pruning report -------------------------------------------------
  if (!is.null(stats_output_file) || !is.null(summary_output_file)) {
    write_pruning_report(
      pruning_stats,
      stats_file   = stats_output_file %||% tempfile(fileext = ".csv"),
      summary_file = summary_output_file
    )
  }

  .report_progress("SLDP finished.", verbose = verbose)

  list(
    final_snp_info           = final_snp_info,
    final_geno_mat           = final_geno_mat,
    screening_stats          = screening_stats,
    candidate_snps_per_trait = candidate_snps_per_trait,
    candidate_snps           = candidate_snps,
    important_snps           = important_snps,
    background_retained      = background_retained,
    pruning_stats            = pruning_stats,
    output_file              = output_file,
    scale_strategy           = strategy
  )
}


#' Extract the final genotype matrix from a GDS-strategy run
#'
#' @param gds_path   Path to the main GDS file.
#' @param snp_ids    SNP IDs to extract.
#' @param sample_ids Optional sample IDs. Default `NULL` (all samples).
#' @return Numeric matrix (SNPs x samples, coded 0/1/2/NA).
#' @export
extract_final_geno <- function(gds_path, snp_ids, sample_ids = NULL) {
  .assert_packages("SNPRelate")
  gf <- .open_gds(gds_path, key = "extract_tmp")
  on.exit(.close_gds(gf, key = "extract_tmp"), add = TRUE)
  .extract_geno_gds(gf, snp_ids = snp_ids, sample_ids = sample_ids)
}
