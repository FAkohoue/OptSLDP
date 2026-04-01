# Read phenotype and covariate data (supports single or multiple traits)

Reads a phenotype file, extracts one or more trait columns plus optional
covariates, and aligns rows to the genotype sample order. All sample
matching is handled internally – the caller never touches row alignment.

## Usage

``` r
read_phenotype(
  file,
  sample_col = "Sample",
  trait_col = "Trait1",
  covar_cols = NULL,
  sample_order = NULL,
  sep = NULL
)
```

## Arguments

- file:

  Path to a delimited text file (CSV, TSV, ...).

- sample_col:

  Column holding sample identifiers. Default `"Sample"`.

- trait_col:

  Trait column name(s). A single string or a character vector for
  multiple traits. Default `"Trait1"`.

- covar_cols:

  Character vector of covariate column names. `NULL` (default) means no
  covariates. Covariates are shared across all traits.

- sample_order:

  Character vector giving the target sample order – i.e., `g$sample_ids`
  from
  [`read_genotype()`](https://FAkohoue.github.io/OptSLDP/reference/read_genotype.md).
  Rows are reordered to match. `NULL` skips reordering.

- sep:

  Separator. `NULL` triggers auto-detection.

## Value

A named list:

- `phenotype`:

  Numeric vector (single trait) **or** numeric matrix (n_samples x
  n_traits) for multiple traits.

- `trait_names`:

  Character vector of returned trait name(s).

- `covariates`:

  Shared covariate `data.frame`, or `NULL`.

- `sample_ids`:

  Sample IDs in the aligned order.

## Multi-trait behaviour

- **Single trait** (`trait_col` is a length-1 string) – returns
  `$phenotype` as a plain numeric vector. Backward compatible with all
  existing single-trait code.

- **Multiple traits** (`trait_col` is a character vector of length \> 1)
  – returns `$phenotype` as a numeric matrix with dimensions (n_samples
  x n_traits), column names equal to `trait_col`, row names equal to the
  aligned sample IDs.

## Examples

``` r
if (FALSE) { # \dontrun{
geno <- read_genotype("geno.csv")

# Single trait
p1 <- read_phenotype("pheno.csv", sample_order = geno$sample_ids)
stopifnot(is.numeric(p1$phenotype) && !is.matrix(p1$phenotype))

# Multiple traits with shared PC covariates
p2 <- read_phenotype(
  "pheno.csv",
  trait_col    = c("BlastLeaf", "BlastPanicle", "Yield"),
  covar_cols   = c("PC1", "PC2"),
  sample_order = geno$sample_ids
)
stopifnot(is.matrix(p2$phenotype))   # n_samples x 3
} # }
```
