# Chromosome-wise high-LD pre-pruning

Performs a fast pre-pruning pass at a very high LD threshold (default
`r^2 >= 0.99`) to remove near-duplicate SNPs before the main SLDP
pipeline.

## Usage

``` r
preprune_high_ld(
  snp_info,
  geno_mat,
  r2_pre = 0.99,
  slide_max_bp = 500000L,
  ctx = NULL,
  verbose = TRUE
)
```

## Arguments

- snp_info:

  Marker metadata `data.table`.

- geno_mat:

  Numeric genotype matrix **or** `NULL` (GDS path).

- r2_pre:

  High-LD removal threshold. Default `0.99`.

- slide_max_bp:

  Sliding window in bp for the GDS path. Default `500 000`.

- ctx:

  Optional GDS context list. Default `NULL`.

- verbose:

  Logical. Default `TRUE`.

## Value

A named list with `snp_info`, `geno_mat` (or `NULL`), and `kept_snps`.

## Details

- **In-memory** (`strategy = "in_memory"` / `"chunked"`): greedy
  position-order loop using
  [`compute_r2_subset()`](https://FAkohoue.github.io/OptSLDP/reference/compute_r2_subset.md).

- **GDS** (`strategy = "gds"`): uses `snpgdsLDpruning()` which processes
  genotypes in streaming blocks – fast enough for 10M+ SNP panels.
