# Find LD neighbors for a focal SNP

Returns all SNPs from `candidate_ids` whose pairwise r^2 with
`focal_snp` meets or exceeds `r2_threshold`. Scale-aware: uses GDS
streaming when a context with `strategy == "gds"` is supplied.

## Usage

``` r
find_ld_neighbors(
  focal_snp,
  geno_mat,
  candidate_ids,
  r2_threshold = 0.9,
  ctx = NULL
)
```

## Arguments

- focal_snp:

  SNP ID of the focal marker.

- geno_mat:

  Numeric genotype matrix **or** `NULL` (GDS path).

- candidate_ids:

  Character vector of SNP IDs to evaluate (may include `focal_snp`; it
  is always excluded from the result).

- r2_threshold:

  Minimum r^2 for neighbor classification. Default `0.90`.

- ctx:

  Optional GDS context list. Default `NULL`.

## Value

Character vector of neighbor SNP IDs (excluding `focal_snp`).

## Examples

``` r
if (FALSE) { # \dontrun{
nbrs <- find_ld_neighbors("rs100", geno_mat, window_snps, r2_threshold = 0.90)
nbrs <- find_ld_neighbors("rs100", NULL,     window_snps, r2_threshold = 0.90,
                          ctx = gds_ctx)
} # }
```
