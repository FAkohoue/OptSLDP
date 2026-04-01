# OptSLDP: An Optimized Selective Linkage Disequilibrium Pruning Pipeline

`OptSLDP` implements an optimized and extended version of the Selective
Linkage Disequilibrium Pruning (SLDP) pipeline for genomic prediction
panel construction. The package builds on the algorithm of Zhu et al.
(2023) and addresses four concrete limitations of the original
implementation: memory-safe chunked genotype reading, explicit
per-chromosome garbage collection, multi-trait union protection, and
covariate-adjusted marginal screening.

## Details

### Core idea

Standard LD pruning (e.g., PLINK `--indep-pairwise`) removes markers
from correlated pairs without regard to whether either marker is
associated with the trait. This discards statistical evidence around QTL
regions. SLDP addresses this by partitioning the marker set before
pruning:

- **Important SNPs** are markers that are statistically associated with
  at least one trait, plus all markers in strong LD with them within a
  positional window. These are retained in full.

- **Background SNPs** are the remainder. Standard greedy LD pruning is
  applied to this subset to control redundancy.

The final panel is the union of all important SNPs and the retained
background markers.

### OptSLDP improvements over original SLDP

|  |  |
|----|----|
| Improvement | Description |
| Chunked reading | Genotype files are read in 50 K-row chunks with a pre-allocated matrix; peak RAM equals one chunk rather than twice the file |
| Per-chromosome [`gc()`](https://rdrr.io/r/base/gc.html) | `gc(FALSE)` is called after each chromosome's r^2 matrix is released in pre-pruning and background-pruning loops |
| Multi-trait union protection | Screening and candidate selection run independently per trait; the union of all per-trait candidate sets drives expansion and protection |
| Covariate-adjusted screening | Phenotypes are residualised on covariates once before the SNP scan, equivalent to including covariates in every marginal regression |

### Automatic scale strategy

The pipeline selects one of three LD backends automatically based on the
post-filter SNP count:

|  |  |  |
|----|----|----|
| Strategy | Default trigger | LD backend |
| `in_memory` | n \<= 200,000 | [`cor()`](https://rdrr.io/r/stats/cor.html) on full matrix in RAM |
| `chunked` | 200,000 \< n \<= 2,000,000 | [`cor()`](https://rdrr.io/r/stats/cor.html) per chromosome |
| `gds` | n \> 2,000,000 | `snpgdsLDMat()` streaming from disk |

Thresholds are adjustable via
`options(optsldp.thresh_small, optsldp.thresh_medium)`.

### Main pipeline

[`run_sldp()`](https://FAkohoue.github.io/OptSLDP/reference/run_sldp.md)
is the primary user-facing function. It accepts a genotype file and a
phenotype file and produces a pruned output genotype file in 12
sequential steps: genotype reading, phenotype alignment, scale strategy
selection, MAF filtering, optional high-LD pre-pruning, marginal
screening, candidate selection, important-SNP expansion, background
pruning, panel merge, output writing, and optional pruning-statistics
report.

## Exported functions

**Main pipeline**

|  |  |
|----|----|
| Function | Role |
| [`run_sldp()`](https://FAkohoue.github.io/OptSLDP/reference/run_sldp.md) | End-to-end pipeline: two files in, one pruned file out |
| [`extract_final_geno()`](https://FAkohoue.github.io/OptSLDP/reference/extract_final_geno.md) | Extract pruned genotype matrix from a GDS-strategy run |

**Input / output**

|  |  |
|----|----|
| Function | Role |
| [`read_genotype()`](https://FAkohoue.github.io/OptSLDP/reference/read_genotype.md) | Auto-dispatch genotype reader (numeric / HapMap / VCF) |
| [`read_numeric_genotype()`](https://FAkohoue.github.io/OptSLDP/reference/read_numeric_genotype.md) | Read numeric dosage CSV/TXT |
| [`read_hapmap_genotype()`](https://FAkohoue.github.io/OptSLDP/reference/read_hapmap_genotype.md) | Read HapMap nucleotide format |
| [`read_vcf_genotype()`](https://FAkohoue.github.io/OptSLDP/reference/read_vcf_genotype.md) | Read VCF v4.2 (plain or bgzipped) |
| [`read_phenotype()`](https://FAkohoue.github.io/OptSLDP/reference/read_phenotype.md) | Read phenotype file; single-trait or multi-trait |
| [`write_numeric_genotype()`](https://FAkohoue.github.io/OptSLDP/reference/write_numeric_genotype.md) | Write numeric dosage output |
| [`write_hapmap_genotype()`](https://FAkohoue.github.io/OptSLDP/reference/write_hapmap_genotype.md) | Write HapMap output |
| [`write_pruning_report()`](https://FAkohoue.github.io/OptSLDP/reference/write_pruning_report.md) | Write step-by-step pruning statistics |

**Quality control**

|  |  |
|----|----|
| Function | Role |
| [`compute_maf()`](https://FAkohoue.github.io/OptSLDP/reference/compute_maf.md) | Compute allele frequency and MAF |
| [`filter_snps_by_maf()`](https://FAkohoue.github.io/OptSLDP/reference/filter_snps_by_maf.md) | Remove SNPs below MAF threshold |
| [`preprune_high_ld()`](https://FAkohoue.github.io/OptSLDP/reference/preprune_high_ld.md) | Chromosome-wise removal of near-duplicate SNPs |

**Screening and candidate selection**

|  |  |
|----|----|
| Function | Role |
| [`compute_screening_stats()`](https://FAkohoue.github.io/OptSLDP/reference/compute_screening_stats.md) | Marginal SNP screening with optional covariate adjustment |
| [`select_candidate_snps()`](https://FAkohoue.github.io/OptSLDP/reference/select_candidate_snps.md) | Select candidates by p-value, effect size, or hybrid mode |

**LD computation and pruning**

|  |  |
|----|----|
| Function | Role |
| [`compute_r2_subset()`](https://FAkohoue.github.io/OptSLDP/reference/compute_r2_subset.md) | Pairwise r^2 for a SNP subset |
| [`find_ld_neighbors()`](https://FAkohoue.github.io/OptSLDP/reference/find_ld_neighbors.md) | LD neighbours of a focal SNP |
| [`expand_important_snps()`](https://FAkohoue.github.io/OptSLDP/reference/expand_important_snps.md) | Expand candidates into the protected set |
| [`prune_background_snps()`](https://FAkohoue.github.io/OptSLDP/reference/prune_background_snps.md) | Greedy LD pruning of background markers |

## Example data

`example_optsldp` is a synthetic dataset containing a 40-SNP x 50-sample
genotype matrix, two simulated traits with known QTN architecture,
pre-computed GWAS results, and ground-truth causal variant tables. Load
with:

    data("example_optsldp", package = "OptSLDP")

## Dependencies

**data.table** (\>= 1.14.0) is used for all tabular data operations.
**parallel** provides `detectCores()` for the GDS backend thread count.
**stats** provides [`lm()`](https://rdrr.io/r/stats/lm.html),
[`cor()`](https://rdrr.io/r/stats/cor.html), and
[`residuals()`](https://rdrr.io/r/stats/residuals.html).

Optional Bioconductor packages extend the package to large panels:
**SNPRelate** and **gdsfmt** enable the GDS streaming backend for panels
with more than 2 million SNPs; **VariantAnnotation**, **GenomeInfoDb**,
**SummarizedExperiment**, **Biostrings**, and **S4Vectors** enable VCF
input.

## References

Zhu D, Zhao Y, Zhang R, Wu H, Cai G, Wu Z, Wang Y, Hu X (2023). Genomic
prediction based on selective linkage disequilibrium pruning of
low-coverage whole-genome sequence variants in a pure Duroc population.
*Genetics Selection Evolution*, **55**(1), 72.
[doi:10.1186/s12711-023-00843-w](https://doi.org/10.1186/s12711-023-00843-w)

Akohoue F, Herrera CC, Balanta SJC, et al. (2026). Enhancing genomic
prediction ability of blast resistance using genome-wide association
study-derived marker weights in two rice (*Oryza sativa* L.)
populations. *Theoretical and Applied Genetics*, **139**, 52.
[doi:10.1007/s00122-026-05159-z](https://doi.org/10.1007/s00122-026-05159-z)

## See also

Useful links:

- <https://github.com/FAkohoue/OptSLDP>

- Report bugs at <https://github.com/FAkohoue/OptSLDP/issues>

## Author

**Maintainer**: Félicien Akohoue <akohoue.f@gmail.com>
