# Read VCF genotype data

Reads a VCF or bgzipped VCF and converts GT fields to additive dosage
(0/1/2/NA). Three code paths are selected automatically:

## Usage

``` r
read_vcf_genotype(
  file,
  vcf_snp_threshold = 50000L,
  clean_malformed = FALSE,
  gds_dir = tempdir(),
  n_cores = 1L,
  verbose = TRUE
)
```

## Arguments

- file:

  Path to a `.vcf` or `.vcf.gz` file.

- vcf_snp_threshold:

  SNP count above which the SNPRelate path is used. Default `50000`.

- clean_malformed:

  Logical. If `TRUE`, stream-clean the VCF before conversion to remove
  lines with unexpected column counts. Adds one extra pass over the file
  but is necessary for some variant callers. Default `FALSE`.

- gds_dir:

  Directory for the temporary GDS file. Default
  [`tempdir()`](https://rdrr.io/r/base/tempfile.html).

- n_cores:

  Threads passed to `snpgdsGetGeno()`. Default `1`.

- verbose:

  Print progress messages. Default `TRUE`.

## Value

Named list: `snp_info`, `geno_mat`, `sample_ids`, `format`.

## Details

- **Small VCF** (fewer than `vcf_snp_threshold` variants, default 50
  000):
  [`VariantAnnotation::readVcf()`](https://rdrr.io/pkg/VariantAnnotation/man/readVcf-methods.html)
  in one pass.

- **Large VCF** (\>= `vcf_snp_threshold` variants):
  [`SNPRelate::snpgdsVCF2GDS()`](https://rdrr.io/pkg/SNPRelate/man/snpgdsVCF2GDS.html)
  streams the VCF to a temporary GDS binary, then dosage is read
  chromosome-by-chromosome. Peak RAM equals one chromosome of dosage.

- **Large VCF with cleaning** (`clean_malformed = TRUE`): before
  conversion, `clean_vcf()` streams through the file and removes any
  lines whose tab-column count does not match the header. Needed for VCF
  files from NGSEP and some other callers that embed "/" in FORMAT
  fields.

All three paths return an identical list so the rest of the pipeline is
unaffected.
