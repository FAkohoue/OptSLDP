# LD-prune background (non-important) SNPs chromosome-wise

Applies greedy LD pruning to all SNPs outside the protected set.

## Usage

``` r
prune_background_snps(
  remaining_snps,
  snp_info,
  geno_mat,
  r2_genome = 0.8,
  slide_max_bp = 1000000L,
  ctx = NULL,
  verbose = TRUE
)
```

## Arguments

- remaining_snps:

  Character vector of non-protected SNP IDs to prune.

- snp_info:

  Marker metadata `data.table`.

- geno_mat:

  Numeric genotype matrix **or** `NULL` (GDS path).

- r2_genome:

  Pruning r^2 threshold. Default `0.80`.

- slide_max_bp:

  Sliding window in bp for the GDS path. Default `1 000 000`.

- ctx:

  Optional GDS context list. Default `NULL`.

- verbose:

  Logical. Default `TRUE`.

## Value

Character vector of retained SNP IDs.

## Details

- **In-memory** (`strategy = "in_memory"` / `"chunked"`): loads the
  chromosome sub-matrix and applies a forward-selection greedy loop.

- **GDS** (`strategy = "gds"`): delegates to `snpgdsLDpruning()` which
  processes genotypes in streaming blocks – the chromosome sub-matrix is
  never fully materialised in RAM.
