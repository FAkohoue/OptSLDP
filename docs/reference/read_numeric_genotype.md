# Read numeric dosage genotype data

Reads a file in the package's standard format:
`SNP, CHR, POS, REF, ALT, sample1, sample2, ...`.

## Usage

``` r
read_numeric_genotype(
  file,
  sep = NULL,
  check_names = FALSE,
  chunk_rows = 50000L,
  chunk_threshold = 200000L,
  clean_malformed = FALSE,
  verbose = TRUE
)
```

## Arguments

- file:

  Path to a delimited file (`.csv`, `.txt`, etc.).

- sep:

  Field separator. `NULL` triggers auto-detection.

- check_names:

  Sanitise column names. Default `FALSE`.

- chunk_rows:

  Rows per chunk for the chunked path. Default `50 000`.

- chunk_threshold:

  Row count above which chunked reading is used. Default `200 000`.

- clean_malformed:

  Logical. If `TRUE`, stream-clean the file before reading by removing
  lines with unexpected column counts. Default `FALSE`.

- verbose:

  Print progress messages when cleaning. Default `TRUE`.

## Value

Named list: `snp_info`, `geno_mat`, `sample_ids`, `format`.

## Details

Files with more than `chunk_threshold` rows are read in chunks of
`chunk_rows` rows so that peak RAM equals one chunk rather than two full
copies of the file. Smaller files use a single `fread()` call.
