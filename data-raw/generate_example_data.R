# ==============================================================================
# data-raw/generate_example_data.R
#
# Generates all example data used by the OptSLDP package:
#   inst/extdata/example_genotypes_numeric.csv  — numeric dosage format
#   inst/extdata/example_genotypes.hmp.txt      — HapMap format
#   inst/extdata/example_genotypes.vcf          — VCF v4.2 (uncompressed)
#   inst/extdata/example_phenotype.csv          — two traits + covariates
#   inst/extdata/example_gwas.csv               — GWAS for Trait1 (diagnostic)
#   data/example_optsldp.rda                       — binary R dataset
#
# Design principles
# -----------------
# Genotype files
#   All three formats encode the SAME underlying genotype matrix so tests can
#   verify format-reader equivalence directly.  Three LD blocks give the
#   pipeline meaningful pruning decisions to make.
#
# Two traits with independent genetic architectures
#   Trait1  QTN: SNP003 (+0.80), SNP015 (-0.50), SNP025 (+0.65)  h2=0.55
#           SNPs sit in LD blocks A, B, C respectively.
#   Trait2  QTN: SNP006 (+0.70), SNP028 (-0.60), SNP020 (+0.45)  h2=0.45
#           SNP006 is in block A, SNP028 is in block C, SNP020 is a
#           singleton outside ALL LD blocks.  This ensures the union-
#           protection logic must protect SNP020 even though it is never
#           a candidate under Trait1 alone.
#   There is NO exact QTN overlap between the two traits (different SNP IDs),
#   so union(qtn1, qtn2) is strictly larger than either set individually.
#
# Phenotype file
#   Columns: Sample, Trait1, Trait2, PC1, PC2.
#   This matches the run_sldp() docstring example:
#     trait_col = c("Trait1", "Trait2"), covar_cols = c("PC1", "PC2")
#
# Backward compatibility
#   The binary .rda retains trait1/trait2 as plain vectors (single-trait use)
#   alongside pheno_matrix (n_samples x 2) for multi-trait use.
#   The old $phenotype element is no longer present; use $trait1 instead.
#
# Run from the package root to regenerate all example data:
#   source("data-raw/generate_example_data.R")
# ==============================================================================

set.seed(42L)

dir.create("data-raw/output", showWarnings = FALSE, recursive = TRUE)
dir.create("inst/extdata",    showWarnings = FALSE, recursive = TRUE)
dir.create("data",            showWarnings = FALSE, recursive = TRUE)


# ==============================================================================
# 1.  PARAMETERS
# ==============================================================================

n_snps    <- 40L
n_samples <- 50L
samples   <- sprintf("Line%02d", seq_len(n_samples))
snps      <- sprintf("SNP%03d", seq_len(n_snps))

# Chr 1: SNP001-SNP020  (positions 10 000 - 200 000)
# Chr 2: SNP021-SNP040  (positions 10 000 - 200 000)
chr <- c(rep("1", 20L), rep("2", 20L))
pos <- rep(as.integer(seq(10000L, 200000L, length.out = 20L)), times = 2L)

bases <- c("A", "C", "G", "T")
ref   <- sample(bases, n_snps, replace = TRUE)
alt   <- vapply(ref, function(r) sample(setdiff(bases, r), 1L), character(1L))


# ==============================================================================
# 2.  GENOTYPE MATRIX WITH REALISTIC LD BLOCKS
# ==============================================================================
# Block A: SNP001-SNP008  chr1  r2~0.80
# Block B: SNP012-SNP018  chr1  r2~0.60
# Block C: SNP023-SNP030  chr2  r2~0.80
# SNP019, SNP020, SNP021, SNP022, SNP031-SNP040: singletons (low LD)

.make_block <- function(n_block, n_samp, r2_target, maf_seed = 0.3) {
  h   <- rbinom(n_samp, 1L, maf_seed)
  mat <- matrix(NA_real_, nrow = n_block, ncol = n_samp)
  for (i in seq_len(n_block)) {
    noise    <- rbinom(n_samp, 1L, 1 - sqrt(r2_target))
    mat[i, ] <- pmin((h + noise) %% 2L + rbinom(n_samp, 1L, 0.1), 2L)
  }
  mat
}

geno_mat <- matrix(
  sample(0:2, n_snps * n_samples, replace = TRUE, prob = c(0.40, 0.40, 0.20)),
  nrow = n_snps, ncol = n_samples,
  dimnames = list(snps, samples)
)
geno_mat[1:8,   ] <- .make_block(8L, n_samples, 0.80, 0.30)  # Block A
geno_mat[12:18, ] <- .make_block(7L, n_samples, 0.60, 0.25)  # Block B
geno_mat[23:30, ] <- .make_block(8L, n_samples, 0.80, 0.35)  # Block C

geno_mat <- matrix(
  as.integer(pmin(pmax(round(geno_mat), 0L), 2L)),
  nrow = n_snps, ncol = n_samples,
  dimnames = list(snps, samples)
)
storage.mode(geno_mat) <- "numeric"

# ~5% missing values distributed randomly
miss_idx           <- sample(length(geno_mat),
                             size = round(0.05 * length(geno_mat)))
geno_mat[miss_idx] <- NA_real_


# ==============================================================================
# 3.  SHARED HELPER: compute genetic value from a QTN list
# ==============================================================================
.genetic_value <- function(qtn_snp_ids, effects) {
  gv <- rep(0, n_samples)
  for (k in seq_along(qtn_snp_ids)) {
    g          <- geno_mat[qtn_snp_ids[k], ]
    g[is.na(g)] <- mean(g, na.rm = TRUE)   # mean-impute for phenotype only
    gv         <- gv + effects[k] * g
  }
  gv
}


# ==============================================================================
# 4.  TRAIT 1  — QTN in blocks A (SNP003), B (SNP015), C (SNP025)
# ==============================================================================
qtn1_snps    <- c("SNP003", "SNP015", "SNP025")
qtn1_effects <- c( 0.80,    -0.50,     0.65)
h2_trait1    <- 0.55

gv1 <- .genetic_value(qtn1_snps, qtn1_effects)
# Polygenic background: small random effects on all non-QTN SNPs
for (s in setdiff(snps, qtn1_snps)) {
  g          <- geno_mat[s, ]
  g[is.na(g)] <- mean(g, na.rm = TRUE)
  gv1        <- gv1 + rnorm(1L, 0, 0.05) * g
}
sigma_e1   <- sqrt(var(gv1) * (1 - h2_trait1) / h2_trait1)
trait1_vec <- gv1 + rnorm(n_samples, sd = sigma_e1)


# ==============================================================================
# 5.  TRAIT 2  — QTN in block A (SNP006), block C (SNP028), singleton (SNP020)
# ==============================================================================
# SNP020 is the 20th SNP on chr1 (position 200 000) — a singleton that lies
# completely outside the three LD blocks.  This forces the union-protection
# step to add SNP020 to the important-SNP set on behalf of Trait2 even though
# Trait1 never nominates it as a candidate.
qtn2_snps    <- c("SNP006", "SNP028", "SNP020")
qtn2_effects <- c( 0.70,    -0.60,     0.45)
h2_trait2    <- 0.45

gv2 <- .genetic_value(qtn2_snps, qtn2_effects)
for (s in setdiff(snps, qtn2_snps)) {
  g          <- geno_mat[s, ]
  g[is.na(g)] <- mean(g, na.rm = TRUE)
  gv2        <- gv2 + rnorm(1L, 0, 0.05) * g
}
sigma_e2   <- sqrt(var(gv2) * (1 - h2_trait2) / h2_trait2)
trait2_vec <- gv2 + rnorm(n_samples, sd = sigma_e2)


# ==============================================================================
# 6.  COMBINED PHENOTYPE DATA FRAME AND MATRIX
# ==============================================================================
pheno_df <- data.frame(
  Sample = samples,
  Trait1 = round(trait1_vec, 4L),
  Trait2 = round(trait2_vec, 4L),
  PC1    = round(rnorm(n_samples), 4L),
  PC2    = round(rnorm(n_samples), 4L),
  stringsAsFactors = FALSE
)

# Named matrix form for direct use in compute_screening_stats()
pheno_matrix           <- cbind(Trait1 = trait1_vec, Trait2 = trait2_vec)
rownames(pheno_matrix) <- samples


# ==============================================================================
# 7.  GWAS RESULTS — marginal OLS for both traits
# ==============================================================================
.run_gwas <- function(y_vec) {
  rows <- lapply(seq_len(n_snps), function(i) {
    g  <- geno_mat[i, ]
    ok <- !is.na(g) & !is.na(y_vec)
    if (sum(ok) < 5L || var(g[ok]) == 0) {
      return(data.frame(SNP = snps[i], CHR = chr[i], POS = pos[i],
                        beta = NA_real_, SE = NA_real_, P.value = NA_real_,
                        stringsAsFactors = FALSE))
    }
    sm <- summary(lm(y_vec[ok] ~ g[ok]))$coefficients
    data.frame(SNP     = snps[i], CHR = chr[i], POS = pos[i],
               beta    = round(sm[2L, 1L], 5L),
               SE      = round(sm[2L, 2L], 5L),
               P.value = sm[2L, 4L],
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  out[order(out$P.value), ]
}

gwas_trait1 <- .run_gwas(trait1_vec)
gwas_trait2 <- .run_gwas(trait2_vec)
gwas_df     <- gwas_trait1   # backward-compatible single name


# ==============================================================================
# 8.  WRITE: NUMERIC CSV
# ==============================================================================
geno_df <- data.frame(
  SNP = snps, CHR = chr, POS = as.integer(pos), REF = ref, ALT = alt,
  as.data.frame(geno_mat, check.names = FALSE),
  check.names = FALSE, stringsAsFactors = FALSE
)
write.csv(geno_df, "inst/extdata/example_genotypes_numeric.csv",
          row.names = FALSE, quote = FALSE)
message("  [ok] inst/extdata/example_genotypes_numeric.csv")


# ==============================================================================
# 9.  WRITE: HAPMAP
# ==============================================================================
hmp_geno <- matrix(NA_character_, nrow = n_snps, ncol = n_samples,
                   dimnames = list(snps, samples))
for (i in seq_len(n_snps)) {
  gi                               <- geno_mat[i, ]
  hmp_geno[i, gi == 0 & !is.na(gi)] <- paste0(ref[i], ref[i])
  hmp_geno[i, gi == 1 & !is.na(gi)] <- paste0(ref[i], alt[i])
  hmp_geno[i, gi == 2 & !is.na(gi)] <- paste0(alt[i], alt[i])
  hmp_geno[i, is.na(gi)]             <- "NN"
}
hapmap_df <- data.frame(
  "rs#"       = snps,    alleles    = paste0(ref, "/", alt),
  chrom       = chr,     pos        = as.integer(pos),
  strand      = "+",     "assembly#" = NA_character_,
  center      = NA_character_, protLSID  = NA_character_,
  assayLSID   = NA_character_, panelLSID = NA_character_,
  QCcode      = NA_character_,
  as.data.frame(hmp_geno, check.names = FALSE),
  check.names = FALSE, stringsAsFactors = FALSE
)
write.table(hapmap_df, "inst/extdata/example_genotypes.hmp.txt",
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NN")
message("  [ok] inst/extdata/example_genotypes.hmp.txt")


# ==============================================================================
# 10. WRITE: VCF v4.2
# ==============================================================================
# SNP001-SNP010: phased GT (|)    SNP011-SNP040: unphased GT (/)
.dosage_to_gt <- function(d, phased) {
  sep <- if (phased) "|" else "/"
  vapply(d, function(x) {
    if (is.na(x)) return("./.")
    switch(as.character(x),
           "0" = paste0("0", sep, "0"),
           "1" = paste0("0", sep, "1"),
           "2" = paste0("1", sep, "1"),
           "./.")
  }, character(1L))
}
vcf_header <- c(
  "##fileformat=VCFv4.2",
  paste0("##fileDate=", format(Sys.Date(), "%Y%m%d")),
  "##source=sldp package example data generator",
  paste0("##reference=synthetic_genome_", format(Sys.Date(), "%Y")),
  '##FILTER=<ID=PASS,Description="All filters passed">',
  '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
  paste(c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT",
          samples), collapse = "\t")
)
vcf_body <- vapply(seq_len(n_snps), function(i) {
  gt_vec <- .dosage_to_gt(geno_mat[i, ], phased = i <= 10L)
  paste(chr[i], as.integer(pos[i]), snps[i], ref[i], alt[i],
        ".", "PASS", ".", "GT", paste(gt_vec, collapse = "\t"), sep = "\t")
}, character(1L))
writeLines(c(vcf_header, vcf_body), "inst/extdata/example_genotypes.vcf")
message("  [ok] inst/extdata/example_genotypes.vcf")


# ==============================================================================
# 11. WRITE: PHENOTYPE CSV  (Trait1, Trait2, PC1, PC2)
# ==============================================================================
write.csv(pheno_df, "inst/extdata/example_phenotype.csv",
          row.names = FALSE, quote = FALSE)
message("  [ok] inst/extdata/example_phenotype.csv")


# ==============================================================================
# 12. WRITE: GWAS CSV  (Trait1 only — diagnostic, not required by pipeline)
# ==============================================================================
write.csv(gwas_df, "inst/extdata/example_gwas.csv",
          row.names = FALSE, quote = FALSE)
message("  [ok] inst/extdata/example_gwas.csv")


# ==============================================================================
# 13. BINARY R DATASET
# ==============================================================================
example_optsldp <- list(

  # Marker metadata
  snp_info = data.frame(SNP = snps, CHR = chr, POS = as.integer(pos),
                        REF = ref, ALT = alt, stringsAsFactors = FALSE),

  # Genotype matrix (SNPs x samples, coded 0/1/2/NA)
  geno_mat   = geno_mat,
  sample_ids = samples,

  # Single-trait phenotype vectors (backward compatible, one per trait)
  trait1     = trait1_vec,
  trait2     = trait2_vec,

  # Multi-trait phenotype matrix (n_samples x 2, colnames = Trait1/Trait2)
  # Pass to compute_screening_stats(geno_mat, y = pheno_matrix) or
  # run_sldp(..., trait_col = c("Trait1","Trait2"))
  pheno_matrix = pheno_matrix,

  # Shared covariates (PC1, PC2)
  covariates = pheno_df[, c("PC1", "PC2"), drop = FALSE],

  # Full phenotype data frame (matches example_phenotype.csv exactly)
  pheno_df = pheno_df,

  # Per-trait pre-computed GWAS results (diagnostic use)
  gwas_trait1 = gwas_trait1,
  gwas_trait2 = gwas_trait2,

  # Ground-truth QTN tables (testing only — unknown in real analyses)
  # Use to assert that the pipeline retains causal SNPs in the final panel.
  qtn_trait1 = data.frame(SNP = qtn1_snps, effect = qtn1_effects,
                          stringsAsFactors = FALSE),
  qtn_trait2 = data.frame(SNP = qtn2_snps, effect = qtn2_effects,
                          stringsAsFactors = FALSE),

  # Numeric data.frame representation (matches CSV on disk)
  genotype_numeric = geno_df,

  description = paste0(
    "Synthetic OptSLDP example dataset (v4). ",
    n_snps, " SNPs x ", n_samples, " samples, 2 chromosomes. ",
    "Trait1 QTN: ", paste(qtn1_snps, collapse="/"),
    " (h2=", h2_trait1, "). ",
    "Trait2 QTN: ", paste(qtn2_snps, collapse="/"),
    " (h2=", h2_trait2, "; SNP020 outside all LD blocks). ",
    "Three LD blocks r2=0.60-0.80. ~5% missing."
  ),
  n_snps    = n_snps,
  n_samples = n_samples,
  seed      = 42L
)

save(example_optsldp, file = "data-raw/output/example_optsldp.RData")
save(example_optsldp, file = "data/example_optsldp.rda",
     compress = "xz", compression_level = 9L)
message("  [ok] data-raw/output/example_optsldp.RData")
message("  [ok] data/example_optsldp.rda")


# ==============================================================================
# 14. SELF-CONSISTENCY CHECKS
# ==============================================================================

# -- Genotype dimensions and value range --------------------------------------
stopifnot(
  nrow(geno_mat) == n_snps,
  ncol(geno_mat) == n_samples,
  all(geno_mat[!is.na(geno_mat)] %in% c(0, 1, 2))
)

# -- Phenotype dimensions -----------------------------------------------------
stopifnot(
  length(trait1_vec)  == n_samples,
  length(trait2_vec)  == n_samples,
  nrow(pheno_matrix)  == n_samples,
  ncol(pheno_matrix)  == 2L,
  identical(colnames(pheno_matrix), c("Trait1", "Trait2")),
  identical(rownames(pheno_matrix), samples),
  nrow(pheno_df)      == n_samples,
  all(c("Sample","Trait1","Trait2","PC1","PC2") %in% names(pheno_df))
)

# -- GWAS dimensions ----------------------------------------------------------
stopifnot(
  nrow(gwas_trait1) == n_snps,
  nrow(gwas_trait2) == n_snps,
  all(c("SNP","CHR","POS","beta","SE","P.value") %in% names(gwas_trait1)),
  all(c("SNP","CHR","POS","beta","SE","P.value") %in% names(gwas_trait2))
)

# -- QTN SNPs present in genotype matrix --------------------------------------
stopifnot(
  all(qtn1_snps %in% rownames(geno_mat)),
  all(qtn2_snps %in% rownames(geno_mat))
)

# -- Trait2 singleton QTN is outside all LD blocks ----------------------------
# Block A = SNP001-SNP008, Block B = SNP012-SNP018; SNP020 must be in neither.
stopifnot(
  !("SNP020" %in% snps[1:8]),
  !("SNP020" %in% snps[12:18])
)

# -- Trait2 block QTN sit inside the correct LD blocks ------------------------
stopifnot(
  "SNP006" %in% snps[1:8],    # Block A
  "SNP028" %in% snps[23:30]   # Block C
)

# -- No exact QTN overlap between traits (union is strictly larger) ------------
stopifnot(
  length(intersect(qtn1_snps, qtn2_snps)) == 0L,
  length(union(qtn1_snps, qtn2_snps)) == length(qtn1_snps) + length(qtn2_snps)
)

# -- All output files exist on disk -------------------------------------------
stopifnot(
  file.exists("inst/extdata/example_genotypes_numeric.csv"),
  file.exists("inst/extdata/example_genotypes.hmp.txt"),
  file.exists("inst/extdata/example_genotypes.vcf"),
  file.exists("inst/extdata/example_phenotype.csv"),
  file.exists("inst/extdata/example_gwas.csv"),
  file.exists("data/example_optsldp.rda")
)

# -- Binary dataset: all required fields present ------------------------------
required_fields <- c(
  "snp_info", "geno_mat", "sample_ids",
  "trait1", "trait2", "pheno_matrix",
  "covariates", "pheno_df",
  "gwas_trait1", "gwas_trait2",
  "qtn_trait1", "qtn_trait2",
  "genotype_numeric"
)
stopifnot(all(required_fields %in% names(example_optsldp)))

# -- Binary dataset: dimensional contracts ------------------------------------
stopifnot(
  length(example_optsldp$trait1)     == n_samples,
  length(example_optsldp$trait2)     == n_samples,
  nrow(example_optsldp$pheno_matrix) == n_samples,
  ncol(example_optsldp$pheno_matrix) == 2L,
  nrow(example_optsldp$qtn_trait1)   == length(qtn1_snps),
  nrow(example_optsldp$qtn_trait2)   == length(qtn2_snps),
  nrow(example_optsldp$gwas_trait1)  == n_snps,
  nrow(example_optsldp$gwas_trait2)  == n_snps
)

# -- HapMap: all missing encoded as "NN", not R NA ----------------------------
stopifnot(!any(is.na(hapmap_df[, samples])))

message("\nAll self-consistency checks passed.")
message(sprintf("  SNPs: %d | Samples: %d | missing: ~5%%", n_snps, n_samples))
message(sprintf("  Trait1 QTN: %s (h2=%.2f)",
                paste(qtn1_snps, collapse=", "), h2_trait1))
message(sprintf("  Trait2 QTN: %s (h2=%.2f) — SNP020 outside all LD blocks",
                paste(qtn2_snps, collapse=", "), h2_trait2))
message("Example data generated successfully.")
