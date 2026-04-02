# Filter SNPs by minor allele frequency

Removes SNPs whose MAF is below `maf_threshold` or is `NA`. When a GDS
context is active the filtering is performed entirely in the GDS layer
without loading genotypes into RAM; when in-memory, the matrix is subset
in place.

## Usage

``` r
filter_snps_by_maf(snp_info, geno_mat, maf_threshold = 0.05, ctx = NULL)
```

## Arguments

- snp_info:

  Marker metadata `data.table`.

- geno_mat:

  Numeric genotype matrix **or** `NULL` (GDS path).

- maf_threshold:

  Minimum MAF. Default `0.05`.

- ctx:

  Optional GDS context list. Default `NULL`.

## Value

A named list:

- `snp_info`:

  Filtered metadata, preserving original order.

- `geno_mat`:

  Filtered matrix (`NULL` for GDS path).

- `maf_table`:

  `data.table` of SNP / AF / MAF for retained markers.
