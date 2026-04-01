# Write HapMap genotype file

Converts numeric dosage values (0/1/2/NA) to nucleotide calls (`"AA"`,
`"AT"`, `"TT"`, `"NN"`) and writes a tab-delimited HapMap file.

## Usage

``` r
write_hapmap_genotype(snp_info, geno_mat, file)
```

## Arguments

- snp_info:

  A data.frame or data.table of SNP metadata with columns `SNP`, `CHR`,
  `POS`, `REF`, `ALT`.

- geno_mat:

  Numeric matrix of dosage values (0/1/2/NA), SNPs x samples.

- file:

  Output file path.

## Value

Invisibly returns `file`.
