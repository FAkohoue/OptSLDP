# Read HapMap genotype data

Reads a HapMap file and converts nucleotide calls (`"AA"`, `"AT"`, ...)
to additive dosage (0/1/2). Large files (\> `chunk_threshold` rows) are
read in chunks; the dosage conversion is applied per chunk so peak RAM
never exceeds one chunk.

## Usage

``` r
read_hapmap_genotype(
  file,
  chunk_rows = 50000L,
  chunk_threshold = 200000L,
  clean_malformed = FALSE,
  verbose = TRUE
)
```

## Arguments

- file:

  Path to a HapMap text file.

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
