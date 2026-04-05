## R CMD check results

0 errors | 0 warnings | 0 notes

Checked on:
- Windows 11 x64, R 4.5.0 (local)
- Ubuntu 20.04 x86_64, R 4.4.1 (CGIAR server)

---

## Test suite

- 212 tests passing, 0 failures, 0 skipped
- Test file: `tests/testthat/test-basic.R`
- Covers: I/O (all 3 formats), MAF filtering, pre-pruning, screening
  (single- and multi-trait), candidate selection (modes A/B/C),
  LD expansion, background pruning, full pipeline (single- and
  multi-trait), C++ kernel functions, GDS backend utilities

---

## Package description

OptSLDP implements an optimized Selective Linkage Disequilibrium Pruning
pipeline for genomic prediction panel construction. The package:

- Reads genotype data in numeric dosage, HapMap, or VCF format
- Performs MAF filtering, optional high-LD pre-pruning, covariate-adjusted
  marginal screening (GLM + PCA), and candidate selection (modes A/B/C)
- Automatic fast PCA after MAF filtering: chromosome-balanced GRM
  eigendecomposition; optional RSpectra backend; no LD pruning pass required
- Expands candidate SNPs into a protected set via positional windows and
  LD neighbours, then applies greedy background pruning
- Scales automatically from small in-memory panels to whole-genome
  sequencing datasets (>10M SNPs) via a GDS/SNPRelate disk-backed backend
- Includes a compiled C++ LD kernel (RcppArmadillo) and C++ screening
  kernel for performance-critical steps

---

## Dependencies

### CRAN packages (Imports)
- data.table (>= 1.14.0)
- parallel
- Rcpp (>= 1.0.0)
- stats, tools (base R)

### Bioconductor packages (Suggests)
- SNPRelate (>= 1.26.0) -- required only for scale_strategy = "gds"
  (datasets > 2M SNPs). Not required for standard use.
- gdsfmt (>= 1.28.0) -- dependency of SNPRelate
- RSpectra -- optional; faster PCA eigendecomposition (CRAN)
- VariantAnnotation, SummarizedExperiment, GenomeInfoDb,
  Biostrings, S4Vectors -- required only for VCF input format

All Bioconductor packages are in Suggests, not Imports. The package
works fully without them for datasets below the GDS threshold and
for non-VCF input formats. Users are prompted with an informative
error if they attempt to use a feature that requires a missing
Bioconductor package.

---

## Notes on non-standard files

- `src/Makevars.win` -- required for BLAS/LAPACK linking on Windows
  (`PKG_LIBS = $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)`)
- `pkgdown/extra_head.html` -- pkgdown site customisation; not part of
  the installed package
- `inst/extdata/` -- five example data files used in vignette and tests

---

## Downstream usage

This package is not yet on CRAN. This is the first submission.
It is in active use for genomic prediction panel construction in
rice breeding programs at CGIAR.
