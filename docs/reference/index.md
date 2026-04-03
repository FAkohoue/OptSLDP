# Package index

## Main pipeline

End-to-end SLDP workflow. Two files in, one pruned genotype file out.

- [`run_sldp()`](https://FAkohoue.github.io/OptSLDP/reference/run_sldp.md)
  : Run the Selective Linkage Disequilibrium Pruning (SLDP) pipeline
- [`extract_final_geno()`](https://FAkohoue.github.io/OptSLDP/reference/extract_final_geno.md)
  : Extract the final genotype matrix from a GDS-strategy run

## Input / output

Readers for numeric dosage, HapMap, and VCF genotype formats; phenotype
reader with multi-trait and sample-alignment support; writers for
numeric dosage and HapMap output; pruning report writer.

- [`read_genotype()`](https://FAkohoue.github.io/OptSLDP/reference/read_genotype.md)
  : Dispatch genotype reading by format
- [`clean_genotype_file()`](https://FAkohoue.github.io/OptSLDP/reference/clean_genotype_file.md)
  : Clean a malformed genotype file by removing lines with wrong column
  counts
- [`read_numeric_genotype()`](https://FAkohoue.github.io/OptSLDP/reference/read_numeric_genotype.md)
  : Read numeric dosage genotype data
- [`read_hapmap_genotype()`](https://FAkohoue.github.io/OptSLDP/reference/read_hapmap_genotype.md)
  : Read HapMap genotype data
- [`read_vcf_genotype()`](https://FAkohoue.github.io/OptSLDP/reference/read_vcf_genotype.md)
  : Read VCF genotype data
- [`read_phenotype()`](https://FAkohoue.github.io/OptSLDP/reference/read_phenotype.md)
  : Read phenotype and covariate data (supports single or multiple
  traits)
- [`write_numeric_genotype()`](https://FAkohoue.github.io/OptSLDP/reference/write_numeric_genotype.md)
  : Write numeric dosage genotype file
- [`write_hapmap_genotype()`](https://FAkohoue.github.io/OptSLDP/reference/write_hapmap_genotype.md)
  : Write HapMap genotype file
- [`write_pruning_report()`](https://FAkohoue.github.io/OptSLDP/reference/write_pruning_report.md)
  : Write pruning statistics report

## Quality control

Minor allele frequency computation and filtering; high-LD pre-pruning to
remove near-duplicate markers before the main pipeline.

- [`compute_maf()`](https://FAkohoue.github.io/OptSLDP/reference/compute_maf.md)
  : Compute allele frequency and minor allele frequency for every SNP
- [`filter_snps_by_maf()`](https://FAkohoue.github.io/OptSLDP/reference/filter_snps_by_maf.md)
  : Filter SNPs by minor allele frequency
- [`preprune_high_ld()`](https://FAkohoue.github.io/OptSLDP/reference/preprune_high_ld.md)
  : Chromosome-wise high-LD pre-pruning

## Screening and candidate selection

Marginal SNP screening with optional covariate adjustment; three
candidate selection modes (p-value, effect size, hybrid).

- [`compute_screening_stats()`](https://FAkohoue.github.io/OptSLDP/reference/compute_screening_stats.md)
  : Compute marginal screening statistics
- [`select_candidate_snps()`](https://FAkohoue.github.io/OptSLDP/reference/select_candidate_snps.md)
  : Select candidate SNPs by screening threshold

## LD computation and pruning

Scale-aware pairwise r² computation; important-SNP expansion by
positional window and LD neighbourhood; background LD pruning.

- [`compute_r2_subset()`](https://FAkohoue.github.io/OptSLDP/reference/compute_r2_subset.md)
  : Compute pairwise LD (r^2) for a subset of SNPs
- [`find_ld_neighbors()`](https://FAkohoue.github.io/OptSLDP/reference/find_ld_neighbors.md)
  : Find LD neighbors for a focal SNP
- [`expand_important_snps()`](https://FAkohoue.github.io/OptSLDP/reference/expand_important_snps.md)
  : Expand candidate SNPs into the full set of important (protected)
  SNPs
- [`prune_background_snps()`](https://FAkohoue.github.io/OptSLDP/reference/prune_background_snps.md)
  : LD-prune background (non-important) SNPs chromosome-wise

## Example data

Synthetic dataset for demonstrating and testing the pipeline. Contains
genotype matrix, phenotype vectors, GWAS results, and ground-truth QTN
tables for 40 SNPs and 50 samples.

- [`example_optsldp`](https://FAkohoue.github.io/OptSLDP/reference/example_optsldp.md)
  : Example data for OptSLDP

## Example data

Synthetic dataset for demonstrating and testing the pipeline. Contains
genotype matrix, phenotype vectors, GWAS results, and ground-truth QTN
tables for 40 SNPs and 50 samples.

- [`example_optsldp`](https://FAkohoue.github.io/OptSLDP/reference/example_optsldp.md)
  : Example data for OptSLDP
