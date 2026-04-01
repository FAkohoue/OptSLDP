# Read VCF genotype data

Reads a VCF or bgzipped VCF using `VariantAnnotation` and converts GT
fields (phased or unphased) to additive dosage (0/1/2).

## Usage

``` r
read_vcf_genotype(file)
```

## Arguments

- file:

  Path to a `.vcf` or `.vcf.gz` file.

## Value

Named list: `snp_info`, `geno_mat`, `sample_ids`, `format`.
