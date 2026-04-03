# Compute pairwise LD (r^2) for a subset of SNPs

Dispatches to an in-memory
[`tcrossprod()`](https://rdrr.io/r/base/crossprod.html) BLAS computation
or a disk-backed `snpgdsLDMat()` call depending on the scale strategy
embedded in `ctx`. When no context is supplied the function falls back
to the pure in-memory path, preserving backward compatibility for direct
calls.

## Usage

``` r
compute_r2_subset(geno_mat, snp_ids, ctx = NULL)
```

## Arguments

- geno_mat:

  Numeric genotype matrix (SNPs x samples, coded 0/1/2/NA) **or** `NULL`
  when `ctx$strategy == "gds"`. When provided directly (no `ctx`), this
  is the only data source used.

- snp_ids:

  Character vector of SNP IDs. All IDs must be present in either
  `rownames(geno_mat)` or the GDS file.

- ctx:

  Optional GDS context list produced by `.build_gds_context()`. When
  `NULL` (default) the function uses `geno_mat` directly.

## Value

Symmetric numeric r^2 matrix with row and column names equal to
`snp_ids`.

## Examples

``` r
if (FALSE) { # \dontrun{
# Direct in-memory call (no context needed for small datasets)
r2 <- compute_r2_subset(geno_mat, snp_ids = c("rs1", "rs2", "rs3"))

# Context-aware call (large datasets)
r2 <- compute_r2_subset(NULL, snp_ids = chr_snps, ctx = gds_ctx)
} # }
```
