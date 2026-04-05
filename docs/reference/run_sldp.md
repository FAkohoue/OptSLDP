# Run the Selective Linkage Disequilibrium Pruning (SLDP) pipeline

End-to-end workflow. The user supplies exactly **two input files** – a
genotype file and a phenotype file – plus an output path. All internal
steps (phenotype reading, sample alignment, GWAS screening, candidate
selection, important-SNP expansion, background pruning, output writing)
are handled automatically.

## Usage

``` r
run_sldp(
  genotype_file,
  phenotype_file,
  output_file,
  sample_col = "Sample",
  trait_col = "Trait1",
  covar_cols = NULL,
  n_pcs = 0L,
  format = c("auto", "numeric", "hapmap", "vcf"),
  output_format = c("numeric", "hapmap"),
  mode = c("A", "B", "C"),
  scale_strategy = NULL,
  gds_dir = tempdir(),
  n_cores = max(1L, parallel::detectCores() - 1L),
  maf_threshold = 0.05,
  preprune_large = TRUE,
  preprune_r2 = 0.99,
  pval_threshold = NULL,
  z_threshold = NULL,
  pve_threshold = NULL,
  threshold_logic = c("AND", "OR"),
  window_kb = 50,
  include_ld_neighbors = TRUE,
  r2_flag = 0.9,
  r2_genome = 0.8,
  slide_max_bp = 1000000L,
  stats_output_file = NULL,
  summary_output_file = NULL,
  clean_malformed = FALSE,
  verbose = TRUE
)
```

## Arguments

- genotype_file:

  Path to genotype file (numeric CSV, HapMap, VCF).

- phenotype_file:

  Path to phenotype file. Must contain a sample ID column (`sample_col`)
  and one or more numeric trait columns (`trait_col`).

- output_file:

  Path for the pruned output genotype file.

- sample_col:

  Sample ID column in the phenotype file. Default `"Sample"`.

- trait_col:

  Trait column name(s). A single string or a character vector for
  multi-trait analysis. Default `"Trait1"`.

- covar_cols:

  Covariate column name(s) in the phenotype file. Covariates are shared
  across all traits. `NULL` = no covariates.

- n_pcs:

  Number of principal components to compute automatically from the
  genotype data and use as covariates. Computed via `snpgdsPCA()` from
  SNPRelate on a LD-pruned SNP subset. Only used when `covar_cols` is
  `NULL`. Set to `0` (default) to disable. Typical values: 3-5 for
  structured populations.

- format:

  Genotype format: `"auto"` (default), `"numeric"`, `"hapmap"`, or
  `"vcf"`.

- output_format:

  Output format: `"numeric"` (default) or `"hapmap"`.

- mode:

  Candidate selection mode: `"A"` (p-value), `"B"` (effect-size), or
  `"C"` (hybrid).

- scale_strategy:

  Override auto-detection: `"in_memory"`, `"chunked"`, or `"gds"`.
  `NULL` = auto.

- gds_dir:

  Directory for intermediate GDS files. Default
  [`tempdir()`](https://rdrr.io/r/base/tempfile.html).

- n_cores:

  Threads for GDS operations. Default `max(1, detectCores() - 1)`.

- maf_threshold:

  Minimum MAF. Default `0.05`.

- preprune_large:

  Apply high-LD pre-pruning. Default `TRUE`.

- preprune_r2:

  r^2 threshold for pre-pruning. Default `0.99`.

- pval_threshold:

  P-value upper bound (modes A / C).

- z_threshold:

  \|z-score\| lower bound (modes B / C).

- pve_threshold:

  PVE lower bound (modes B / C).

- threshold_logic:

  `"AND"` (default) or `"OR"`.

- window_kb:

  Positional expansion window in kb. Default `50`.

- include_ld_neighbors:

  Include LD neighbours. Default `TRUE`.

- r2_flag:

  r^2 threshold for neighbour inclusion. Default `0.90`.

- r2_genome:

  r^2 threshold for background pruning. Default `0.80`.

- slide_max_bp:

  Sliding window (bp) for GDS pruning. Default `1 000 000`.

- stats_output_file:

  Optional path for pruning statistics CSV.

- summary_output_file:

  Optional path for plain-text summary.

- clean_malformed:

  If `TRUE`, stream-clean the genotype file before reading by removing
  any lines whose column count does not match the header. Works for all
  accepted formats (numeric CSV, HapMap, VCF). Needed for files from
  NGSEP and other callers that produce malformed lines. Adds one extra
  streaming pass. Default `FALSE`.

- verbose:

  Print timestamped progress. Default `TRUE`.

## Value

A named list:

- `final_snp_info`:

  `data.table` of retained SNP metadata.

- `final_geno_mat`:

  Retained genotype matrix (`NULL` for GDS runs).

- `screening_stats`:

  Single `data.table` (one trait) **or** named list of `data.table`s
  (multiple traits).

- `candidate_snps_per_trait`:

  Named list of per-trait candidate SNP IDs (only present in multi-trait
  runs; `NULL` otherwise).

- `candidate_snps`:

  Union of candidate SNPs across all traits (or the single-trait
  candidate set).

- `important_snps`:

  All protected SNP IDs.

- `background_retained`:

  Retained background SNP IDs.

- `pruning_stats`:

  Step-by-step count `data.table`.

- `output_file`:

  Path to the written output.

- `scale_strategy`:

  Strategy used.

## Multi-trait analysis

Set `trait_col` to a character vector containing more than one column
name from the phenotype file. The pipeline then:

1.  Runs screening (marginal regression + residualisation)
    **separately** for each trait.

2.  Selects candidate SNPs independently for each trait using the shared
    threshold settings.

3.  Takes the **union** of all per-trait candidate sets before
    expansion.

4.  Expands the union into the protected set (positional window + LD
    neighbours) and runs background pruning once on the remainder.

The final pruned panel therefore retains every SNP that is important for
**any** of the requested traits. The return value includes per-trait
screening statistics and per-trait candidate lists alongside the single
shared output file.

## Scale strategy

|  |  |  |
|----|----|----|
| Strategy | Trigger (default) | LD backend |
| `in_memory` | n_snps \<= 200 000 | [`cor()`](https://rdrr.io/r/stats/cor.html) on full matrix in RAM |
| `chunked` | 200K \< n_snps \<= 2M | [`cor()`](https://rdrr.io/r/stats/cor.html) per chromosome in RAM |
| `gds` | n_snps \> 2M | `snpgdsLDMat()` from disk |

Adjustable via `options(optsldp.thresh_small, optsldp.thresh_medium)`.

## Examples

``` r
if (FALSE) { # \dontrun{
geno_file  <- system.file("extdata", "example_genotypes_numeric.csv",
                           package = "OptSLDP")
pheno_file <- system.file("extdata", "example_phenotype.csv",
                           package = "OptSLDP")

# -- Single trait, numeric output (default) ----------------------------------
res <- run_sldp(
  genotype_file  = geno_file,
  phenotype_file = pheno_file,
  output_file    = tempfile(fileext = ".csv"),
  trait_col      = "Trait1",
  mode           = "A",
  pval_threshold = 0.05,
  output_format  = "numeric"    # default: values are 0/1/2/NA
)

# -- Single trait, HapMap output ---------------------------------------------
res_hmp <- run_sldp(
  genotype_file  = geno_file,
  phenotype_file = pheno_file,
  output_file    = tempfile(fileext = ".hmp.txt"),  # use .hmp.txt extension
  trait_col      = "Trait1",
  mode           = "A",
  pval_threshold = 0.05,
  output_format  = "hapmap"     # nucleotide calls: AA/AT/TT/NN
)

# -- Multi-trait with automatic PCA (population structure correction) -------
res <- run_sldp(
  genotype_file  = geno_file,
  phenotype_file = pheno_file,
  output_file    = tempfile(fileext = ".csv"),
  trait_col      = c("Trait1", "Trait2"),
  n_pcs          = 3L,             # auto-compute 3 PCs from genotypes
  mode           = "A",
  pval_threshold = 0.05,
  output_format  = "numeric"
)
# Inspect per-trait results
names(res$screening_stats)          # "Trait1" "Trait2"
res$candidate_snps_per_trait$Trait1

# -- Multi-trait with user-supplied PCs ------------------------------------
res <- run_sldp(
  genotype_file  = geno_file,
  phenotype_file = pheno_file,
  output_file    = tempfile(fileext = ".csv"),
  trait_col      = c("Trait1", "Trait2"),
  covar_cols     = c("PC1", "PC2"),  # columns already in phenotype file
  mode           = "A",
  pval_threshold = 0.05,
  output_format  = "numeric"
)
names(res$screening_stats)          # "Trait1" "Trait2"
res$candidate_snps_per_trait$Trait1
} # }
```
