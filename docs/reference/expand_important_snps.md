# Expand candidate SNPs into the full set of important (protected) SNPs

Builds the protected marker set in two layers:

1.  **Positional window** – every SNP within `window_kb` kb of each
    candidate (same chromosome) is added.

2.  **LD neighborhood** (optional) – SNPs in the window with r^2 \>=
    `r2_flag` against the focal candidate are also added.

## Usage

``` r
expand_important_snps(
  candidate_snps,
  snp_info,
  geno_mat,
  window_kb = 50,
  include_ld_neighbors = TRUE,
  r2_flag = 0.9,
  ctx = NULL,
  verbose = TRUE
)
```

## Arguments

- candidate_snps:

  Character vector of candidate SNP IDs.

- snp_info:

  Marker metadata `data.table`.

- geno_mat:

  Numeric genotype matrix **or** `NULL` (GDS path).

- window_kb:

  Positional expansion window in kb (each side). Default `50`.

- include_ld_neighbors:

  Logical; include LD neighbors. Default `TRUE`.

- r2_flag:

  r^2 threshold for neighbor inclusion. Default `0.90`.

- ctx:

  Optional GDS context list. Default `NULL`.

- verbose:

  Logical. Default `TRUE`.

## Value

Character vector of unique important SNP IDs (superset of
`candidate_snps`).

## Details

The LD computation dispatches to the GDS backend automatically when a
GDS context (`ctx`) is supplied.
