#' Example data for OptSLDP
#'
#' Synthetic dataset for demonstrating and testing the OptSLDP pipeline.
#' All components are internally consistent: SNP IDs, sample IDs, phenotype
#' vectors, covariate columns, GWAS results, and ground-truth QTN tables all
#' refer to the same set of 40 synthetic markers and 50 synthetic samples, so
#' any combination of components can be passed to package functions without
#' modification.
#'
#' @name example_optsldp
#' @docType data
#' @usage data(example_optsldp)
#' @keywords datasets
#'
#' @format A named list with the following components:
#' \describe{
#'   \item{`snp_info`}{Data frame with columns `SNP`, `CHR`, `POS`, `REF`,
#'     and `ALT`. Marker metadata for 40 SNPs on 2 chromosomes (20 SNPs
#'     each), positions 10,000-200,000 bp.}
#'   \item{`geno_mat`}{Numeric matrix (40 SNPs x 50 samples) of additive
#'     dosage values coded as 0, 1, 2, or `NA`. Approximately 5% of values
#'     are missing at random. Row names are SNP IDs; column names are sample
#'     IDs. The same genotype matrix is also provided in three file formats
#'     in `inst/extdata/`.}
#'   \item{`sample_ids`}{Character vector of 50 sample IDs (`Line01` to
#'     `Line50`). Matches `colnames(geno_mat)` and the `Sample` column in
#'     `pheno_df`.}
#'   \item{`trait1`}{Numeric vector of length 50. Simulated phenotype for
#'     Trait1 (blast resistance proxy, heritability \eqn{h^2 = 0.55}).
#'     QTN: SNP003 (+0.80), SNP015 (-0.50), SNP025 (+0.65).}
#'   \item{`trait2`}{Numeric vector of length 50. Simulated phenotype for
#'     Trait2 (heritability \eqn{h^2 = 0.45}).
#'     QTN: SNP006 (+0.70), SNP028 (-0.60), SNP020 (+0.45). SNP020 is a
#'     singleton outside all LD blocks, specifically included to test
#'     multi-trait union-protection in `run_sldp()`.}
#'   \item{`pheno_matrix`}{Numeric matrix (50 samples x 2 traits) with
#'     column names `Trait1` and `Trait2` and row names matching
#'     `sample_ids`. Pass directly to `compute_screening_stats()` or
#'     to `run_sldp()` via `trait_col = c("Trait1", "Trait2")`.}
#'   \item{`covariates`}{Data frame (50 rows x 2 columns: `PC1`, `PC2`).
#'     Simulated principal-component covariates. Pass to `run_sldp()` via
#'     `covar_cols = c("PC1", "PC2")`.}
#'   \item{`pheno_df`}{Data frame (50 rows x 5 columns: `Sample`, `Trait1`,
#'     `Trait2`, `PC1`, `PC2`). The complete phenotype table; matches
#'     `inst/extdata/example_phenotype.csv` exactly.}
#'   \item{`gwas_trait1`}{Data frame with columns `SNP`, `CHR`, `POS`,
#'     `beta`, `SE`, and `P.value`. Pre-computed marginal OLS results for
#'     Trait1, sorted by ascending p-value. For diagnostic use only; not
#'     required by `run_sldp()`.}
#'   \item{`gwas_trait2`}{Data frame with the same structure as
#'     `gwas_trait1` for Trait2.}
#'   \item{`qtn_trait1`}{Data frame with columns `SNP` and `effect`.
#'     Ground-truth causal variants for Trait1 (SNP003, SNP015, SNP025).
#'     For testing only -- unknown in real analyses.}
#'   \item{`qtn_trait2`}{Data frame with columns `SNP` and `effect`.
#'     Ground-truth causal variants for Trait2 (SNP006, SNP028, SNP020).
#'     For testing only -- unknown in real analyses.}
#'   \item{`genotype_numeric`}{Data frame representation of the genotype
#'     matrix in the standard numeric-dosage format (columns: `SNP`, `CHR`,
#'     `POS`, `REF`, `ALT`, then one column per sample). Matches
#'     `inst/extdata/example_genotypes_numeric.csv` exactly.}
#'   \item{`description`}{Character scalar. Plain-text summary of the
#'     dataset design and key parameters.}
#'   \item{`n_snps`}{Integer scalar. Number of SNPs (40).}
#'   \item{`n_samples`}{Integer scalar. Number of samples (50).}
#'   \item{`seed`}{Integer scalar. Random seed used to generate the dataset
#'     (42L). Reproducibility: re-run
#'     `source("data-raw/generate_example_data.R")` with this seed to
#'     regenerate identical data.}
#' }
#'
#' @details
#' All data are synthetic and generated solely for use in package examples,
#' unit tests, and the introductory vignette. No real breeding trial data are
#' included.
#'
#' ## Genotype structure
#'
#' The 40-SNP panel contains three LD blocks and a set of singletons:
#'
#' | Block | SNPs | Chromosome | Target \eqn{r^2} |
#' |-------|------|-----------|-----------|
#' | A | SNP001-SNP008 | 1 | ~ 0.80 |
#' | B | SNP012-SNP018 | 1 | ~ 0.60 |
#' | C | SNP023-SNP030 | 2 | ~ 0.80 |
#' | Singletons | remaining SNPs | 1 & 2 | low |
#'
#' Approximately 5% of dosage values are set to `NA` at random, simulating
#' typical low-coverage sequencing or genotyping array missingness.
#'
#' ## Genetic architecture
#'
#' Trait1 and Trait2 share no causal variants (no QTN overlap), so their
#' union is strictly larger than either individual set. SNP020 -- the
#' Trait2-only singleton QTN -- lies outside all three LD blocks, which
#' means it is retained in the final panel only if the multi-trait
#' union-protection logic operates correctly.
#'
#' ## Illustrative workflow
#'
#' ```r
#' data("example_optsldp", package = "OptSLDP")
#' x <- example_optsldp
#'
#' # Quick single-trait screen using the pre-loaded matrix
#' stats <- compute_screening_stats(x$geno_mat, y = x$trait1, verbose = FALSE)
#' head(stats[order(stats$P.value), ])
#'
#' # Multi-trait run via file interface
#' geno_file  <- system.file("extdata", "example_genotypes_numeric.csv",
#'                            package = "OptSLDP")
#' pheno_file <- system.file("extdata", "example_phenotype.csv",
#'                            package = "OptSLDP")
#'
#' res <- run_sldp(
#'   genotype_file  = geno_file,
#'   phenotype_file = pheno_file,
#'   output_file    = tempfile(fileext = ".csv"),
#'   trait_col      = c("Trait1", "Trait2"),
#'   covar_cols     = c("PC1", "PC2"),
#'   mode           = "A",
#'   pval_threshold = 0.05,
#'   verbose        = FALSE
#' )
#'
#' # Verify QTN retention
#' all(x$qtn_trait1$SNP %in% res$final_snp_info$SNP)
#' all(x$qtn_trait2$SNP %in% res$final_snp_info$SNP)
#' ```
#'
#' @source Synthetically generated for package documentation and testing.
#'   Regenerate with `source("data-raw/generate_example_data.R")` from the
#'   package root (seed = 42).
#'
#' @seealso
#'   `run_sldp()` for the end-to-end pipeline,
#'   `compute_screening_stats()` for marginal SNP screening,
#'   `select_candidate_snps()` for candidate selection modes.
"example_optsldp"
