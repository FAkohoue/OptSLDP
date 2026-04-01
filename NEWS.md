# OptSLDP 0.1.0

First public release of `OptSLDP`, an optimized and extended implementation of
the Selective Linkage Disequilibrium Pruning (SLDP) pipeline for genomic
prediction panel construction. The package builds on the algorithm of
Zhu et al. (2023, *Genetics Selection Evolution*, **55**(1), 72,
<doi:10.1186/s12711-023-00843-w>) and is further motivated by the
GWAS-derived marker weighting framework of Akohoue et al. (2026,
*Theoretical and Applied Genetics*, **139**, 52,
<doi:10.1007/s00122-026-05159-z>).

---

## New functions

### Main pipeline

* `run_sldp()` — end-to-end pipeline accepting two input files (a genotype
  file and a phenotype file) and one output path. Internally handles phenotype
  reading, sample alignment, MAF filtering, optional high-LD pre-pruning,
  marginal screening, candidate selection, important-SNP expansion, background
  LD pruning, and output writing across 12 sequential steps. Scale strategy
  (in-memory / chunked / GDS) is selected automatically from the post-filter
  SNP count or overridden explicitly.

* `extract_final_geno()` — extracts the final retained genotype matrix from a
  SNPRelate GDS file after a GDS-strategy `run_sldp()` run, for cases where
  the in-memory matrix is needed post-pipeline without re-running the full
  workflow.

### I/O — genotype readers

* `read_genotype()` — auto-dispatch wrapper that detects the genotype file
  format from the file extension (`.vcf` / `.vcf.gz` → VCF; `.hmp.txt` →
  HapMap; everything else → numeric dosage) and calls the appropriate
  format-specific reader.

* `read_numeric_genotype()` — reads files in the package's standard tabular
  format (`SNP, CHR, POS, REF, ALT, sample1, sample2, ...`, values 0/1/2/NA).
  Files with more than `chunk_threshold` rows (default 200 000) are processed
  by a two-pass chunked strategy: a header-only scan pre-allocates the full
  output matrix, then data are filled in `chunk_rows`-row slices (default
  50 000) so that peak RAM equals one chunk rather than twice the full file.

* `read_hapmap_genotype()` — reads standard HapMap format, converting
  two-character nucleotide calls (`"AA"`, `"AT"`, `"TT"`, `"NN"`) to additive
  dosage (0/1/2/NA). Uses the same two-pass chunked strategy as
  `read_numeric_genotype()` for large files.

* `read_vcf_genotype()` — reads VCF v4.2 files (plain or bgzipped) via
  `VariantAnnotation`. Accepts phased (`0|1`) and unphased (`0/1`) GT fields.
  Multi-allelic sites use the first ALT allele; missing calls (`./.`) become
  `NA`. Requires Bioconductor packages `VariantAnnotation`, `GenomeInfoDb`,
  `SummarizedExperiment`, `Biostrings`, and `S4Vectors`.

### I/O — phenotype reader

* `read_phenotype()` — reads a delimited phenotype file, extracts one or more
  trait columns and optional shared covariates, and aligns rows to the
  genotype sample order. A single trait column returns a plain numeric vector
  (backward compatible); a character vector of column names returns a named
  numeric matrix (n_samples × n_traits) for multi-trait analysis. Duplicate
  sample IDs, mismatched sample sets, and non-numeric trait values all produce
  informative errors or warnings.

### I/O — writers

* `write_numeric_genotype()` — writes the pruned genotype matrix and SNP
  metadata in the standard `SNP, CHR, POS, REF, ALT, sample...` CSV format.

* `write_hapmap_genotype()` — converts the pruned dosage matrix back to
  two-character nucleotide calls and writes a tab-delimited HapMap file.

* `write_pruning_report()` — writes the step-by-step pruning statistics
  `data.table` to a CSV file, and optionally writes a human-readable
  plain-text summary showing SNP counts at every pipeline step.

### Quality control

* `compute_maf()` — computes ALT allele frequency (AF) and minor allele
  frequency (MAF) for every SNP from the dosage matrix. AF is estimated as
  `rowMeans(geno_mat, na.rm = TRUE) / 2`; MAF as `min(AF, 1 - AF)`.
  Scale-aware: delegates to `snpgdsSNPRateFreq()` when a GDS context is
  active.

* `filter_snps_by_maf()` — removes SNPs whose MAF is below `maf_threshold`
  (default 0.05) or is `NA`. Returns filtered `snp_info`, `geno_mat`, and
  a MAF table. Scale-aware.

* `preprune_high_ld()` — chromosome-wise greedy removal of near-duplicate
  SNPs at a high r² threshold (default 0.99) before the main pipeline.
  Intended to catch collapsed haplotypes and duplicate probes, not biological
  LD. Scale-aware. Explicit `rm()` + `gc(FALSE)` after each chromosome's r²
  matrix releases allocator pressure across chromosome passes.

### Screening and candidate selection

* `compute_screening_stats()` — fits a marginal linear regression of the
  phenotype on each SNP. When covariates are supplied, the phenotype is
  residualised on them once before the scan (equivalent to including
  covariates in every regression but computed in a single `lm()` call).
  Returns per-SNP `beta`, `SE`, `z_score`, `P.value`, `PVE`, `AF`, and `MAF`.
  Accepts a numeric vector (single trait) or a numeric matrix (multi-trait),
  returning a single `data.table` or a named list of `data.table`s
  respectively.

* `select_candidate_snps()` — applies one of three selection modes to a
  screening statistics table:
  * **Mode A** — p-value threshold: `p ≤ τ_p`
  * **Mode B** — effect-size criteria: `|z| ≥ τ_z` and/or `R² ≥ τ_pve`,
    combined by `"AND"` or `"OR"` logic
  * **Mode C** — hybrid: any combination of modes A and B criteria

### LD computation and pruning

* `compute_r2_subset()` — computes the pairwise r² matrix for a SNP subset.
  In-memory / chunked path uses `cor()` with pairwise-complete observations;
  GDS path uses `snpgdsLDMat()` which streams genotypes from disk. Scale-aware.

* `find_ld_neighbors()` — returns the IDs of all SNPs from a candidate set
  whose pairwise r² with a focal SNP meets or exceeds `r2_threshold`. Used
  internally during important-SNP expansion. Scale-aware.

* `expand_important_snps()` — expands the candidate set into a protected set
  in two layers: (1) all SNPs within ±`window_kb` kb on the same chromosome,
  (2) optionally, SNPs in that window with r² ≥ `r2_flag` against the focal
  candidate. The protected set is the union over all candidates. Scale-aware.

* `prune_background_snps()` — greedy forward-selection LD pruning of all
  SNPs outside the protected set, applied chromosome-by-chromosome in genomic
  position order. Discards SNP j > i if r²(i, j) ≥ `r2_genome` (default
  0.80). Explicit `rm()` + `gc(FALSE)` after each chromosome's r² matrix.
  Scale-aware; GDS path delegates to `snpgdsPruneLD()`.

---

## Key design features

### Multi-trait union protection

`run_sldp()` accepts `trait_col` as a character vector of any length.
When multiple traits are specified:

1. Marginal screening and candidate selection run independently for each
   trait (separate phenotype residualisation, separate p-values, separate
   candidate sets).
2. The union of all per-trait candidate sets is formed:
   `C* = C_1 ∪ C_2 ∪ ... ∪ C_K`.
3. The union set is expanded into the protected set and background pruning
   runs once on the remainder.

Every SNP important for any trait is therefore retained in the final panel.
The return value carries per-trait screening statistics (`$screening_stats`,
a named list) and per-trait candidate vectors
(`$candidate_snps_per_trait`, a named list) alongside the shared final panel.

### Memory-safe chunked genotype reading

Genotype files exceeding `chunk_threshold` rows (default 200 000) are read by
a two-pass strategy. Pass 1 reads only the header to determine dimensions and
pre-allocates `matrix(NA_real_, nrow, ncol)`. Pass 2 reads `chunk_rows`
(default 50 000) rows at a time and writes each block directly into the
pre-allocated matrix slice. Peak RAM equals one chunk (~400 MB for a 50 K ×
3 K panel) rather than twice the file size (the standard `fread` + `as.matrix`
path). `gc(FALSE)` is called after each chunk is released.

### Explicit per-chromosome garbage collection

`gc(FALSE)` is called after every chromosome's r² matrix is released in
`preprune_high_ld()` and `prune_background_snps()`. This prevents allocator
fragmentation from accumulating across 20–30 chromosome passes — an issue
that causes apparent memory leaks on large panels without explicitly triggering
collection between passes.

### Automatic scale strategy selection

The pipeline selects one of three LD backends based on the post-filter SNP
count, with site-wide adjustable thresholds:

| Strategy | Default trigger | LD backend |
|----------|----------------|-----------|
| `in_memory` | n ≤ 200 000 | `cor()` on full matrix in RAM |
| `chunked` | 200 000 < n ≤ 2 000 000 | `cor()` per chromosome |
| `gds` | n > 2 000 000 | `snpgdsLDMat()` streaming from disk |

Thresholds are configurable via
`options(optsldp.thresh_small, optsldp.thresh_medium)`.

The GDS backend writes the genotype matrix to a binary SNPRelate GDS file
once, streams all LD queries from disk, and closes all file handles
automatically via `on.exit()` and a package-level handle registry.

### Chromosome name normalisation

All three genotype readers normalise chromosome names at read time by stripping
the leading `"chr"` prefix (case-insensitive), so `"chr1"` and `"1"` are
treated as the same chromosome throughout the pipeline without requiring manual
harmonisation.

### Covariate-adjusted screening

When `covar_cols` are specified, the phenotype is projected onto the covariate
space once before the SNP scan:
`y* = y - X_c (X_c' X_c)^{-1} X_c' y`
This is equivalent to including covariates in every marginal regression but
avoids refitting covariates for each of the thousands of SNPs, giving an
O(n_covariates) saving over the naïve approach.

---

## Example data

The package ships five files in `inst/extdata/` encoding the same underlying
40-SNP × 50-sample synthetic dataset in all supported formats:

* `example_genotypes_numeric.csv` — numeric dosage (0/1/2/NA)
* `example_genotypes.hmp.txt` — HapMap nucleotide calls
* `example_genotypes.vcf` — VCF v4.2; SNP001–SNP010 phased (`|`),
  SNP011–SNP040 unphased (`/`)
* `example_phenotype.csv` — two traits (`Trait1`, `Trait2`) and two
  covariates (`PC1`, `PC2`) for 50 samples
* `example_gwas.csv` — pre-computed GWAS results for `Trait1` (diagnostic
  only; not required by `run_sldp()`)

The dataset has a known genetic architecture designed to exercise all pipeline
branches:

* Three LD blocks: Block A (SNP001–SNP008, chr1, r² ≈ 0.80), Block B
  (SNP012–SNP018, chr1, r² ≈ 0.60), Block C (SNP023–SNP030, chr2, r² ≈ 0.80)
* ~5% missing values across all three file formats
* **Trait1** QTN: SNP003 (+0.80), SNP015 (−0.50), SNP025 (+0.65); h² = 0.55;
  one QTN per LD block
* **Trait2** QTN: SNP006 (+0.70), SNP028 (−0.60), SNP020 (+0.45); h² = 0.45;
  SNP020 is a singleton outside all LD blocks, specifically to test that
  union-protection retains it via the Trait2 candidate set

---

## Package infrastructure

* `data-raw/generate_example_data.R` — reproducible script (seed 42) that
  generates all five `inst/extdata/` files and the `data/example_optsldp.rda`
  binary dataset from first principles, with 22 self-consistency assertions.
* `vignettes/OptSLDP-introduction.Rmd` — full introductory vignette covering
  all functions, all three candidate selection modes, single-trait and
  multi-trait workflows, scale strategies, output writers, and memory notes;
  includes 43 named code chunks and 14 LaTeX display-math blocks.
* `pkgdown/extra_head.html` — custom `<head>` content for the pkgdown site:
  Google Fonts preconnect hints, Open Graph / Twitter Card meta tags, MathJax 3
  delimiter configuration, bslib-aligned CSS refinements, and a
  copy-to-clipboard button for code blocks.
* `_pkgdown.yml` — Bootstrap 5 pkgdown site configuration with Inter /
  Fira Code fonts, navy-and-blue colour scheme, and a five-section function
  reference index.
* `.github/workflows/pkgdown.yaml` — GitHub Actions workflow that installs
  Bioconductor dependencies, builds the pkgdown site with
  `build_site_github_pages()`, and deploys to GitHub Pages via
  `actions/deploy-pages@v4`.

---

## References

Zhu D, Zhao Y, Zhang R, Wu H, Cai G, Wu Z, Wang Y, Hu X (2023). Genomic
prediction based on selective linkage disequilibrium pruning of low-coverage
whole-genome sequence variants in a pure Duroc population. *Genetics Selection
Evolution*, **55**(1), 72. <https://doi.org/10.1186/s12711-023-00843-w>

Akohoue F, Herrera CC, Balanta SJC, et al. (2026). Enhancing genomic
prediction ability of blast resistance using genome-wide association
study-derived marker weights in two rice (*Oryza sativa* L.) populations.
*Theoretical and Applied Genetics*, **139**, 52.
<https://doi.org/10.1007/s00122-026-05159-z>
