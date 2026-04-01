# Compute allele frequency and minor allele frequency for every SNP

For the in-memory path, estimates AF from the dosage matrix as
`rowMeans / 2`. For the GDS path, delegates to `snpgdsSNPRateFreq()`
which streams genotypes from disk without loading the full matrix.

## Usage

``` r
compute_maf(geno_mat, ctx = NULL)
```

## Arguments

- geno_mat:

  Numeric genotype matrix (SNPs x samples, coded 0/1/2/NA) **or** `NULL`
  when `ctx$strategy == "gds"`.

- ctx:

  Optional GDS context list. Default `NULL`.

## Value

A `data.table` with columns `SNP`, `AF`, and `MAF`. When using the GDS
path, `SNP` values are the SNP IDs stored in the file.

## Examples

``` r
if (FALSE) { # \dontrun{
maf_dt <- compute_maf(geno_mat)
} # }
```
