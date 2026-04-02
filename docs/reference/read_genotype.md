# Dispatch genotype reading by format

Dispatch genotype reading by format

## Usage

``` r
read_genotype(file, format = c("auto", "numeric", "hapmap", "vcf"))
```

## Arguments

- file:

  Path to the genotype file.

- format:

  One of `"auto"` (default), `"numeric"`, `"hapmap"`, `"vcf"`.

## Value

Named list from the format-specific reader.
