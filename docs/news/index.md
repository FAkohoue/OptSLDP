# Changelog

## OptSLDP 0.1.0

### Performance improvements

- **C++ screening kernel** (`screen_chunk_cpp()`, `src/ld_kernel.cpp`) –
  single-pass OLS per SNP chunk in compiled C++. Genotype variance is
  computed once per chunk and reused across all traits, eliminating
  repeated `rowSums(g_dev^2)` calls. Pure-R fallback preserved. Typical
  gain: 2-4x over R matrix algebra.

- **SNP chunking in step 6** – GDS screening now processes 50,000-SNP
  chunks within each chromosome instead of loading entire chromosomes.
  Improves CPU cache locality and reduces peak allocation size per
  chunk.

- **Chromosome-streaming screening (step 6)** – GDS strategy extracts
  and screens one chunk at a time, accumulating only the statistics
  table (not genotypes). Peak RAM drops from ~4 GB (full 2.65M x 204
  extraction) to ~300 MB (one 50k-SNP chunk). Step 7 uses pre-computed
  statistics – no genotype access needed at candidate selection.

- **Chromosome-streaming output writing (step 11)** – final panel
  written chromosome by chromosome from GDS using `fwrite(append=TRUE)`.
  The full final genotype matrix is never loaded into RAM.
  `$final_geno_mat` is `NULL` for GDS runs; use
  [`extract_final_geno()`](https://FAkohoue.github.io/OptSLDP/reference/extract_final_geno.md)
  to recover it if needed.

- **Parallel background pruning (step 9)** – chromosomes pruned
  simultaneously using a FORK cluster on Linux (one independent
  read-only GDS handle per worker). With 8 cores and 11 chromosomes,
  step 9 time drops from ~66 minutes to ~8-10 minutes. Falls back to
  sequential on Windows.

- **`method = "corr"` enforced in `snpgdsPruneLD`** – both
  `.prune_background_chr_gds()` and `.preprune_chr_gds()` now pass
  `method = "corr"` to `snpgdsLDpruning`, ensuring Pearson r is used
  consistently with the `r2_genome` and `r2_pre` thresholds everywhere
  in the pipeline. The previous default (`"composite"`, Weir 1979) is a
  different estimator that diverges from Pearson r at low MAF.

### Bug fixes

- **`r2_matrix_cpp` reference removed from `ld.R`** – the dead
  `tryCatch(r2_matrix_cpp, ...)` branch in `.r2_tcrossprod()` referenced
  a function that was never exported. Removed;
  [`stats::sd`](https://rdrr.io/r/stats/sd.html) called explicitly.

- **`RcppArmadillo` removed from `Imports`** – compile-time only;
  belongs in `LinkingTo` only. Previously triggered an `R CMD check`
  NOTE.

- **`start` and `end` added to
  [`globalVariables()`](https://rdrr.io/r/utils/globalVariables.html)**
  – data.table column names used in `setkey()` inside `pruning.R`;
  previously triggered an `R CMD check` NOTE.

- **`output_format = "hapmap"` example fixed** –
  [`run_sldp()`](https://FAkohoue.github.io/OptSLDP/reference/run_sldp.md)
  roxygen example now explicitly shows `output_format = "hapmap"` with
  `.hmp.txt` file extension. Previous example omitted `output_format`,
  silently defaulting to `"numeric"`.

- **README broken R syntax fixed** –
  `res_mt\`candidate_snps_per_trait\`Trait1`backtick notation replaced with correct`res_mt\\candidate_snps_per_trait\\Trait1\`
  dollar-sign access throughout README and vignette.

- **README LaTeX rendering fixed** – `\left\{` and `\right\}` replaced
  with plain `\{` and `\}` throughout. GitHub MathJax cannot parse
  curly-brace delimiters with `\left`/`\right`. Statistics table
  replaced with plain-text formulas to avoid complex `\widehat`, `\frac`
  inside markdown table cells.

### Infrastructure

- **`utils_cpp.R` rewritten** – `.get_cpp_fun()` / `.cpp_available()`
  helper uses `asNamespace("OptSLDP")` for robust function lookup under
  both
  [`devtools::load_all()`](https://rdrr.io/pkg/devtools/man/load_all.html)
  and installed builds. `.screen_chunk()` wrapper added for
  `screen_chunk_cpp()`. All four C++ wrappers now have pure-R fallbacks
  that produce correct results (not just empty matrices).

- **`screen_chunk_cpp` registered** – added to `src/RcppExports.cpp` and
  `R/RcppExports.R`. Run
  [`Rcpp::compileAttributes()`](https://rdrr.io/pkg/Rcpp/man/compileAttributes.html)
  after installing to pick up the new export.

- **`@importFrom stats sd`** added to `OptSLDP-package.R` and
  `@importFrom stats sd` added to `ld.R` so `NAMESPACE` is complete.

------------------------------------------------------------------------

## OptSLDP 0.1.0

### Performance

- **Vectorised OLS screening** – `.screen_single_trait()` uses matrix
  algebra (`sweep` + `tcrossprod` via BLAS DGEMM) instead of per-SNP
  [`lm()`](https://rdrr.io/r/stats/lm.html) calls. Speedup: 50-200x for
  the screening step on large panels.

- **Parallel multi-trait screening** – traits are screened
  simultaneously using a FORK cluster on Linux when `n_cores > 1`,
  halving screening time for two-trait runs. Falls back to sequential on
  Windows.

- **C++ LD kernel** (`src/ld_kernel.cpp`, RcppArmadillo) – three
  compiled functions eliminate R interpreter overhead from all LD hot
  paths:

  - `r2_subset_cpp()`: candidate-row-only r^2 via BLAS DGEMM.
  - `above_threshold_subset_cpp()`: sparse (row, col) scan in C++.
  - `greedy_prune_r2_cpp()`: greedy forward-selection pruning in C++.
    Pure-R fallbacks provided for systems without a C++ compiler.

- **Batched chromosome LD expansion** –
  [`expand_important_snps()`](https://FAkohoue.github.io/OptSLDP/reference/expand_important_snps.md)
  extracts genotypes once per chromosome via `.extract_geno_gds()` and
  uses
  [`data.table::foverlaps()`](https://rdrr.io/pkg/data.table/man/foverlaps.html)
  for O(n log n) positional window joins. Reduces GDS reads from
  n_candidates (thousands) to n_chromosomes (~11 for rice).

- **Skip-on-rerun GDS caching** – VCF-to-GDS conversion and
  `sldp_main.gds` write check for existing files and skip if present.
  Reruns skip ~8 minutes of one-time preparation.

- **[`tcrossprod()`](https://rdrr.io/r/base/crossprod.html) replaces
  [`stats::cor()`](https://rdrr.io/r/stats/cor.html)** – all in-memory
  LD computations use mean-imputation +
  [`tcrossprod()`](https://rdrr.io/r/base/crossprod.html) (BLAS DGEMM)
  via `.r2_tcrossprod()`.

### New functions

- **[`clean_genotype_file()`](https://FAkohoue.github.io/OptSLDP/reference/clean_genotype_file.md)**
  – stream-clean any genotype file (numeric CSV, HapMap, VCF/bgzipped
  VCF) by removing lines with unexpected column counts.

- **[`run_sldp()`](https://FAkohoue.github.io/OptSLDP/reference/run_sldp.md)
  gains `clean_malformed` parameter** – pass `TRUE` to clean before
  reading. Temporary file deleted automatically after reading.

### Infrastructure

- `Rcpp` and `RcppArmadillo` added to `Imports` and `LinkingTo`.
- `R/utils_cpp.R` – R-level wrappers around C++ functions with pure-R
  fallbacks.
- `tools` added to `Imports`.

------------------------------------------------------------------------

First public release. See README for full feature list and references.
