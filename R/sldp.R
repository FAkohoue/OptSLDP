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
#' # -- Single trait ------------------------------------------------------------
#' res <- run_sldp(
#'   genotype_file  = geno_file,
#'   phenotype_file = pheno_file,
#'   output_file    = tempfile(fileext = ".csv"),
#'   trait_col      = "Trait1",
#'   mode           = "A",
#'   pval_threshold = 0.05
#' )
#'
#' # -- Multi-trait: union protection -------------------------------------------
#' res <- run_sldp(
#'   genotype_file  = geno_file,
#'   phenotype_file = pheno_file,
#'   output_file    = tempfile(fileext = ".csv"),
#'   trait_col      = c("Trait1", "Trait2"),
#'   covar_cols     = c("PC1", "PC2"),
#'   mode           = "A",
#'   pval_threshold = 0.05
#' )
#' # Inspect per-trait results
#' names(res$screening_stats)          # "Trait1" "Trait2"
#' res$candidate_snps_per_trait$Trait1
#' }
run_sldp <- function(genotype_file,
                     phenotype_file,
                     output_file,
                     sample_col             = "Sample",
                     trait_col              = "Trait1",
                     covar_cols             = NULL,
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
                     verbose                = TRUE) {

  format          <- match.arg(format)
  output_format   <- match.arg(output_format)
  mode            <- match.arg(mode)
  threshold_logic <- match.arg(threshold_logic)
  trait_col       <- as.character(trait_col)
  multi_trait     <- length(trait_col) > 1L

  pruning_stats <- .new_pruning_stats()
  on.exit(.close_all_gds_handles(), add = TRUE)
  .report_progress("SLDP started", verbose = verbose)

  # -- Step 1: Read genotypes --------------------------------------------------
  .report_progress("[1] Reading genotype data ...", verbose = verbose)
  g        <- read_genotype(genotype_file, format = format)
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
  # For GDS runs the post-filter subset is extracted into RAM here.
  # It is far smaller than the original 10M+ panel.
  .report_progress("[6] Computing screening statistics ...", verbose = verbose)

  geno_screen <- if (identical(strategy, "gds")) {
    .report_progress("    Extracting ", nrow(snp_info),
                     " post-filter SNPs from GDS ...", verbose = verbose)
    .extract_geno_gds(ctx$genofile,
                      snp_ids    = snp_info$SNP,
                      sample_ids = g$sample_ids)
  } else {
    geno_mat
  }

  # -- Step 7: Candidate selection (single-trait or multi-trait) --------------
  if (multi_trait) {
    # -- Multi-trait branch ----------------------------------------------------
    # Screening runs separately per trait; candidates are unioned.
    .report_progress("[7] Multi-trait screening (",
                     length(trait_col), " traits) ...", verbose = verbose)

    screen_res <- .screen_all_traits(
      geno_mat        = geno_screen,
      pheno_matrix    = phenotype,       # n_samples x n_traits matrix
      covar           = covariates,
      mode            = mode,
      pval_threshold  = pval_threshold,
      z_threshold     = z_threshold,
      pve_threshold   = pve_threshold,
      threshold_logic = threshold_logic,
      verbose         = verbose
    )

    screening_stats          <- screen_res$screening_stats
    candidate_snps_per_trait <- screen_res$candidate_snps_per_trait
    candidate_snps           <- screen_res$candidate_snps_union

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

    screening_stats <- .screen_single_trait(
      geno_screen, y = phenotype, covar = covariates, verbose = verbose
    )
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

  final_geno_mat <- if (identical(strategy, "gds")) {
    .extract_geno_gds(ctx$genofile,
                      snp_ids    = final_snp_info$SNP,
                      sample_ids = g$sample_ids)
  } else {
    geno_mat[final_snp_info$SNP, , drop = FALSE]
  }

  # -- Step 11: Write output ---------------------------------------------------
  .report_progress("[11] Writing output to: ", output_file, verbose = verbose)
  if (identical(output_format, "numeric"))
    write_numeric_genotype(final_snp_info, final_geno_mat, output_file)
  else
    write_hapmap_genotype(final_snp_info, final_geno_mat, output_file)

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
