# Clean a malformed genotype file by removing lines with wrong column counts

Streams through any accepted genotype file format (numeric CSV, HapMap,
VCF / bgzipped VCF) in 50 000-line chunks and writes only lines whose
delimiter-separated column count matches the header to a new output
file. The separator is detected automatically from the file extension
and the first non-comment header line:

## Usage

``` r
clean_genotype_file(
  file_in,
  file_out,
  sep = "auto",
  chunk_size = 50000L,
  verbose = TRUE
)
```

## Arguments

- file_in:

  Path to the input genotype file (plain or `.gz`).

- file_out:

  Path for the cleaned output file. If it ends in `.gz` the output is
  gzip-compressed.

- sep:

  Field separator: `","`, `"\t"`, or `"auto"` (default). `"auto"`
  detects from the file extension (`vcf` / `hmp` -\> tab, everything
  else -\> comma).

- chunk_size:

  Lines per chunk. Default `50000`.

- verbose:

  Print progress every 500 000 data lines. Default `TRUE`.

## Value

Invisibly: a named integer vector with elements `total`, `removed`, and
`kept`.

## Details

|         |           |
|---------|-----------|
| Format  | Separator |
| VCF     | tab       |
| HapMap  | tab       |
| Numeric | comma     |

VCF comment lines (starting with `#`) are always passed through
unchanged; the column count is established from the `#CHROM` line. For
numeric and HapMap files the first line is the header and sets the
expected count.
