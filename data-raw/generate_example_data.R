# ==============================================================================
# data-raw/generate_example_data.R
#
# Generates all example data used by the OptSLDP package:
#   inst/extdata/example_genotypes_numeric.csv  -- numeric dosage format
#   inst/extdata/example_genotypes.hmp.txt      -- HapMap format
#   inst/extdata/example_genotypes.vcf          -- VCF v4.2 (uncompressed)
#   inst/extdata/example_phenotype.csv          -- two traits + covariates
#   inst/extdata/example_gwas.csv               -- GWAS for Trait1 (diagnostic)
#   data/example_optsldp.rda                    -- binary R dataset
#
# Run from the package root:
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

chr <- c(rep("1", 20L), rep("2", 20L))
pos <- rep(as.integer(seq(10000L, 200000L, length.out = 20L)), times = 2L)

bases <- c("A", "C", "G", "T")
ref   <- sample(bases, n_snps, replace = TRUE)
alt   <- vapply(ref, function(r) sample(setdiff(bases, r), 1L), character(1L))

# ==============================================================================
# 2.  GENOTYPE MATRIX WITH REALISTIC LD BLOCKS
# ==============================================================================

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
geno_mat[1:8,   ] <- .make_block(8L, n_samples, 0.80, 0.30)
geno_mat[12:18, ] <- .make_block(7L, n_samples, 0.60, 0.25)
geno_mat[23:30, ] <- .make_block(8L, n_samples, 0.80, 0.35)

geno_mat <- matrix(
  as.integer(pmin(pmax(round(geno_mat), 0L), 2L)),
  nrow = n_snps, ncol = n_samples,
  dimnames = list(snps, samples)
)
storage.mode(geno_mat) <- "numeric"

miss_idx           <- sample(length(geno_mat),
                             size = round(0.05 * length(geno_mat)))
geno_mat[miss_idx] <- NA_real_

# ==============================================================================
# 3.  HELPER: genetic value from QTN list
# ==============================================================================

.genetic_value <- function(qtn_snp_ids, effects) {
  gv <- rep(0, n_samples)
  for (k in seq_along(qtn_snp_ids)) {
    g           <- geno_mat[qtn_snp_ids[k], ]
    g[is.na(g)] <- mean(g, na.rm = TRUE)
    gv          <- gv + effects[k] * g
  }
  gv
}

# ==============================================================================
# 4.  TRAIT 1 -- QTN in blocks A (SNP003), B (SNP015), C (SNP025)
# ==============================================================================

qtn1_snps    <- c("SNP003", "SNP015", "SNP025")
qtn1_effects <- c( 0.80,    -0.50,     0.65)
h2_trait1    <- 0.55

gv1 <- .genetic_value(qtn1_snps, qtn1_effects)
for (s in setdiff(snps, qtn1_snps)) {
  g           <- geno_mat[s, ]
  g[is.na(g)] <- mean(g, na.rm = TRUE)
  gv1         <- gv1 + rnorm(1L, 0, 0.05) * g
}
sigma_e1   <- sqrt(var(gv1) * (1 - h2_trait1) / h2_trait1)
trait1_vec <- gv1 + rnorm(n_samples, sd = sigma_e1)

# ==============================================================================
# 5.  TRAIT 2 -- QTN in block A (SNP006), block C (SNP028), singleton (SNP020)
# ==============================================================================

qtn2_snps    <- c("SNP006", "SNP028", "SNP020")
qtn2_effects <- c( 0.70,    -0.60,     0.45)
h2_trait2    <- 0.45

gv2 <- .genetic_value(qtn2_snps, qtn2_effects)
for (s in setdiff(snps, qtn2_snps)) {
  g           <- geno_mat[s, ]
  g[is.na(g)] <- mean(g, na.rm = TRUE)
  gv2         <- gv2 + rnorm(1L, 0, 0.05) * g
}
sigma_e2   <- sqrt(var(gv2) * (1 - h2_trait2) / h2_trait2)
trait2_vec <- gv2 + rnorm(n_samples, sd = sigma_e2)

# ==============================================================================
# 6.  PHENOTYPE DATA FRAME AND MATRIX
# ==============================================================================

pheno_df <- data.frame(
  Sample = samples,
  Trait1 = round(trait1_vec, 4L),
  Trait2 = round(trait2_vec, 4L),
  PC1    = round(rnorm(n_samples), 4L),
  PC2    = round(rnorm(n_samples), 4L),
  stringsAsFactors = FALSE
)

pheno_matrix           <- cbind(Trait1 = trait1_vec, Trait2 = trait2_vec)
rownames(pheno_matrix) <- samples

# ==============================================================================
# 7.  GWAS RESULTS
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
gwas_df     <- gwas_trait1

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
  gi                                <- geno_mat[i, ]
  hmp_geno[i, gi == 0 & !is.na(gi)] <- paste0(ref[i], ref[i])
  hmp_geno[i, gi == 1 & !is.na(gi)] <- paste0(ref[i], alt[i])
  hmp_geno[i, gi == 2 & !is.na(gi)] <- paste0(alt[i], alt[i])
  hmp_geno[i, is.na(gi)]             <- "NN"
}
hapmap_df <- data.frame(
  "rs#"        = snps,    alleles     = paste0(ref, "/", alt),
  chrom        = chr,     pos         = as.integer(pos),
  strand       = "+",     "assembly#" = NA_character_,
  center       = NA_character_, protLSID  = NA_character_,
  assayLSID    = NA_character_, panelLSID = NA_character_,
  QCcode       = NA_character_,
  as.data.frame(hmp_geno, check.names = FALSE),
  check.names = FALSE, stringsAsFactors = FALSE
)
write.table(hapmap_df, "inst/extdata/example_genotypes.hmp.txt",
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NN")
message("  [ok] inst/extdata/example_genotypes.hmp.txt")

# ==============================================================================
# 10. WRITE: VCF v4.2
# ==============================================================================

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
  "##source=OptSLDP package example data",
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
# 11. WRITE: PHENOTYPE CSV
# ==============================================================================

write.csv(pheno_df, "inst/extdata/example_phenotype.csv",
          row.names = FALSE, quote = FALSE)
message("  [ok] inst/extdata/example_phenotype.csv")

# ==============================================================================
# 12. WRITE: GWAS CSV
# ==============================================================================

write.csv(gwas_df, "inst/extdata/example_gwas.csv",
          row.names = FALSE, quote = FALSE)
message("  [ok] inst/extdata/example_gwas.csv")

# ==============================================================================
# 13. BINARY R DATASET
# ==============================================================================

example_optsldp <- list(
  snp_info         = data.frame(SNP = snps, CHR = chr, POS = as.integer(pos),
                                REF = ref, ALT = alt, stringsAsFactors = FALSE),
  geno_mat         = geno_mat,
  sample_ids       = samples,
  trait1           = trait1_vec,
  trait2           = trait2_vec,
  pheno_matrix     = pheno_matrix,
  covariates       = pheno_df[, c("PC1", "PC2"), drop = FALSE],
  pheno_df         = pheno_df,
  gwas_trait1      = gwas_trait1,
  gwas_trait2      = gwas_trait2,
  qtn_trait1       = data.frame(SNP = qtn1_snps, effect = qtn1_effects,
                                stringsAsFactors = FALSE),
  qtn_trait2       = data.frame(SNP = qtn2_snps, effect = qtn2_effects,
                                stringsAsFactors = FALSE),
  genotype_numeric = geno_df,
  description      = paste0(
    "Synthetic OptSLDP example dataset. ",
    n_snps, " SNPs x ", n_samples, " samples, 2 chromosomes. ",
    "Trait1 QTN: ", paste(qtn1_snps, collapse = "/"),
    " (h2=", h2_trait1, "). ",
    "Trait2 QTN: ", paste(qtn2_snps, collapse = "/"),
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
message("  [ok] data/example_optsldp.rda")

# ==============================================================================
# 14. SELF-CONSISTENCY CHECKS
# ==============================================================================

stopifnot(
  nrow(geno_mat) == n_snps,
  ncol(geno_mat) == n_samples,
  all(geno_mat[!is.na(geno_mat)] %in% c(0, 1, 2)),
  length(trait1_vec)  == n_samples,
  length(trait2_vec)  == n_samples,
  nrow(pheno_matrix)  == n_samples,
  ncol(pheno_matrix)  == 2L,
  identical(colnames(pheno_matrix), c("Trait1", "Trait2")),
  nrow(gwas_trait1)   == n_snps,
  nrow(gwas_trait2)   == n_snps,
  all(qtn1_snps %in% rownames(geno_mat)),
  all(qtn2_snps %in% rownames(geno_mat)),
  !("SNP020" %in% snps[1:8]),
  !("SNP020" %in% snps[12:18]),
  "SNP006" %in% snps[1:8],
  "SNP028" %in% snps[23:30],
  length(intersect(qtn1_snps, qtn2_snps)) == 0L,
  file.exists("inst/extdata/example_genotypes_numeric.csv"),
  file.exists("inst/extdata/example_genotypes.hmp.txt"),
  file.exists("inst/extdata/example_genotypes.vcf"),
  file.exists("inst/extdata/example_phenotype.csv"),
  file.exists("inst/extdata/example_gwas.csv"),
  file.exists("data/example_optsldp.rda"),
  all(c("snp_info","geno_mat","sample_ids","trait1","trait2",
        "pheno_matrix","covariates","pheno_df",
        "gwas_trait1","gwas_trait2","qtn_trait1","qtn_trait2",
        "genotype_numeric") %in% names(example_optsldp)),
  !any(is.na(hapmap_df[, samples]))
)

message("\nAll self-consistency checks passed.")
message(sprintf("  SNPs: %d | Samples: %d | missing ~5%%", n_snps, n_samples))
message(sprintf("  Trait1 QTN: %s (h2=%.2f)", paste(qtn1_snps, collapse=", "), h2_trait1))
message(sprintf("  Trait2 QTN: %s (h2=%.2f)", paste(qtn2_snps, collapse=", "), h2_trait2))
message("Example data generated successfully.")
