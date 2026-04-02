# Extract the final genotype matrix from a GDS-strategy run

Extract the final genotype matrix from a GDS-strategy run

## Usage

``` r
extract_final_geno(gds_path, snp_ids, sample_ids = NULL)
```

## Arguments

- gds_path:

  Path to the main GDS file.

- snp_ids:

  SNP IDs to extract.

- sample_ids:

  Optional sample IDs. Default `NULL` (all samples).

## Value

Numeric matrix (SNPs x samples, coded 0/1/2/NA).
