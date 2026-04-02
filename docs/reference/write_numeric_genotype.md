# Write numeric dosage genotype file

Write numeric dosage genotype file

## Usage

``` r
write_numeric_genotype(snp_info, geno_mat, file, sep = ",")
```

## Arguments

- snp_info:

  A data.frame or data.table of SNP metadata with columns `SNP`, `CHR`,
  `POS`, `REF`, `ALT`.

- geno_mat:

  Numeric matrix of dosage values (0/1/2/NA), SNPs x samples.

- file:

  Output file path.

- sep:

  Field separator. Default `","`.

## Value

Invisibly returns `file`.
