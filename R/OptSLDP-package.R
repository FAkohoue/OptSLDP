#' OptSLDP: An Optimized Selective Linkage Disequilibrium Pruning Pipeline
#'
#' `OptSLDP` implements an optimized and extended version of the Selective
#' Linkage Disequilibrium Pruning (SLDP) pipeline for genomic prediction panel
#' construction. The package builds on the algorithm of Zhu et al. (2023) and
#' addresses four concrete limitations of the original implementation:
#' memory-safe chunked genotype reading, explicit per-chromosome garbage
#' collection, multi-trait union protection, and covariate-adjusted marginal
#' screening.
#'
#' @details
#' ## Core idea
#'
#' Standard LD pruning (e.g., PLINK `--indep-pairwise`) removes markers from
#' correlated pairs without regard to whether either marker is associated with
#' the trait. This discards statistical evidence around QTL regions. SLDP
#' addresses this by partitioning the marker set before pruning:
#'
#' - **Important SNPs** are markers that are statistically associated with at
#'   least one trait, plus all markers in strong LD with them within a
#'   positional window. These are retained in full.
#' - **Background SNPs** are the remainder. Standard greedy LD pruning is
#'   applied to this subset to control redundancy.
#'
#' The final panel is the union of all important SNPs and the retained
#' background markers.
#'
#' ## OptSLDP improvements over original SLDP
#'
#' | Improvement | Description |
#' |---|---|
#' | Chunked reading | Genotype files are read in 50 K-row chunks with a pre-allocated matrix; peak RAM equals one chunk rather than twice the file |
#' | Per-chromosome `gc()` | `gc(FALSE)` is called after each chromosome's r^2 matrix is released in pre-pruning and background-pruning loops |
#' | Multi-trait union protection | Screening and candidate selection run independently per trait; the union of all per-trait candidate sets drives expansion and protection |
#' | Covariate-adjusted screening | Phenotypes are residualised on covariates once before the SNP scan, equivalent to including covariates in every marginal regression |
#'
#' ## Automatic scale strategy
#'
#' The pipeline selects one of three LD backends automatically based on the
#' post-filter SNP count:
#'
#' | Strategy | Default trigger | LD backend |
#' |---|---|---|
#' | `in_memory` | n <= 200,000 | `cor()` on full matrix in RAM |
#' | `chunked` | 200,000 < n <= 2,000,000 | `cor()` per chromosome |
#' | `gds` | n > 2,000,000 | `snpgdsLDMat()` streaming from disk |
#'
#' Thresholds are adjustable via
#' `options(optsldp.thresh_small, optsldp.thresh_medium)`.
#'
#' ## Main pipeline
#'
#' `run_sldp()` is the primary user-facing function. It accepts a genotype
#' file and a phenotype file and produces a pruned output genotype file in
#' 12 sequential steps: genotype reading, phenotype alignment, scale strategy
#' selection, MAF filtering, optional high-LD pre-pruning, marginal screening,
#' candidate selection, important-SNP expansion, background pruning, panel
#' merge, output writing, and optional pruning-statistics report.
#'
#' @section Exported functions:
#'
#' **Main pipeline**
#'
#' | Function | Role |
#' |---|---|
#' | `run_sldp()` | End-to-end pipeline: two files in, one pruned file out |
#' | `extract_final_geno()` | Extract pruned genotype matrix from a GDS-strategy run |
#'
#' **Input / output**
#'
#' | Function | Role |
#' |---|---|
#' | `read_genotype()` | Auto-dispatch genotype reader (numeric / HapMap / VCF) |
#' | `read_numeric_genotype()` | Read numeric dosage CSV/TXT |
#' | `read_hapmap_genotype()` | Read HapMap nucleotide format |
#' | `read_vcf_genotype()` | Read VCF v4.2 (plain or bgzipped) |
#' | `read_phenotype()` | Read phenotype file; single-trait or multi-trait |
#' | `write_numeric_genotype()` | Write numeric dosage output |
#' | `write_hapmap_genotype()` | Write HapMap output |
#' | `write_pruning_report()` | Write step-by-step pruning statistics |
#'
#' **Quality control**
#'
#' | Function | Role |
#' |---|---|
#' | `clean_genotype_file()` | Remove malformed lines from any genotype file |
#' | `compute_maf()` | Compute allele frequency and MAF |
#' | `filter_snps_by_maf()` | Remove SNPs below MAF threshold |
#' | `preprune_high_ld()` | Chromosome-wise removal of near-duplicate SNPs |
#'
#' **Screening and candidate selection**
#'
#' | Function | Role |
#' |---|---|
#' | `compute_screening_stats()` | Marginal SNP screening with optional covariate adjustment |
#' | `select_candidate_snps()` | Select candidates by p-value, effect size, or hybrid mode |
#'
#' **LD computation and pruning**
#'
#' | Function | Role |
#' |---|---|
#' | `compute_r2_subset()` | Pairwise r^2 for a SNP subset |
#' | `find_ld_neighbors()` | LD neighbours of a focal SNP |
#' | `expand_important_snps()` | Expand candidates into the protected set |
#' | `prune_background_snps()` | Greedy LD pruning of background markers |
#'
#' @section Example data:
#' `example_optsldp` is a synthetic dataset containing a 40-SNP x 50-sample
#' genotype matrix, two simulated traits with known QTN architecture,
#' pre-computed GWAS results, and ground-truth causal variant tables. Load
#' with:
#'
#' ```r
#' data("example_optsldp", package = "OptSLDP")
#' ```
#'
#' @section Dependencies:
#' **data.table** (>= 1.14.0) is used for all tabular data operations.
#' **parallel** provides `detectCores()` for the GDS backend thread count.
#' **stats** provides `lm()`, `cor()`, and `residuals()`.
#'
#' Optional Bioconductor packages extend the package to large panels:
#' **SNPRelate** and **gdsfmt** enable the GDS streaming backend for panels
#' with more than 2 million SNPs; **VariantAnnotation**, **GenomeInfoDb**,
#' **SummarizedExperiment**, **Biostrings**, and **S4Vectors** enable VCF
#' input.
#'
#' @references
#' Zhu D, Zhao Y, Zhang R, Wu H, Cai G, Wu Z, Wang Y, Hu X (2023). Genomic
#' prediction based on selective linkage disequilibrium pruning of
#' low-coverage whole-genome sequence variants in a pure Duroc population.
#' \emph{Genetics Selection Evolution}, \strong{55}(1), 72.
#' \doi{10.1186/s12711-023-00843-w}
#'
#' Akohoue F, Herrera CC, Balanta SJC, et al. (2026). Enhancing genomic
#' prediction ability of blast resistance using genome-wide association
#' study-derived marker weights in two rice (\emph{Oryza sativa} L.)
#' populations. \emph{Theoretical and Applied Genetics}, \strong{139}, 52.
#' \doi{10.1007/s00122-026-05159-z}
#'
#' @useDynLib OptSLDP, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom stats sd
#' @keywords internal
"_PACKAGE"
