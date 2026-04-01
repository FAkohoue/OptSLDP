# Compute marginal screening statistics

Fits a marginal linear regression of the phenotype on each SNP. When
covariates are supplied, `y` is first residualised on them so the scan
is covariate-adjusted.

## Usage

``` r
compute_screening_stats(geno_mat, y, covar = NULL, verbose = TRUE)
```

## Arguments

- geno_mat:

  Numeric genotype matrix (SNPs x samples, coded 0/1/2/NA).

- y:

  Numeric phenotype vector **or** numeric matrix (n_samples x n_traits)
  aligned to `colnames(geno_mat)`.

- covar:

  Optional covariate matrix / data.frame (nrow == nrow(y)). Applied to
  every trait independently. Default `NULL`.

- verbose:

  Print progress messages. Default `TRUE`.

## Value

Single `data.table` (vector `y`) **or** named list of `data.table`s
(matrix `y`). Each table has columns: SNP, beta, SE, z_score, P.value,
PVE, AF, MAF.

## Multi-trait input

When `y` is a numeric matrix (n_samples x n_traits) the function is
called once per column and returns a **named list** of per-trait
`data.table`s. Column names of `y` become the names of the list. When
`y` is a plain numeric vector the function returns a single `data.table`
as before (backward compatible).

## Examples

``` r
if (FALSE) { # \dontrun{
# Single trait
stats <- compute_screening_stats(geno_mat, y = pheno_vec)

# Multiple traits -- returns named list
stats_list <- compute_screening_stats(geno_mat, y = pheno_matrix)
stats_list[["Yield"]]
} # }
```
