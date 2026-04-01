# ==============================================================================
# io.R
# Genotype I/O, phenotype I/O, sample harmonisation, and format dispatch.
#
# Changes in v4
# -------------
# Inspiration 2 -- chunked genotype readers:
#   read_numeric_genotype() and read_hapmap_genotype() pre-allocate the full
#   output matrix from a one-line header scan and fill it in fixed-size row
#   chunks (default 50 000 rows).  Peak RAM ~ one chunk, never two full copies
#   of the file.  Files under `chunk_threshold` rows still use the fast
#   single-pass fread() path.
#
# Inspiration 5 -- gc(FALSE) after every chunk:
#   Each chunk helper calls gc(FALSE) immediately after releasing the chunk
#   data.table.  This prevents fragmentation from accumulating across the
#   dozens of chromosome-level passes that happen later in the pipeline.
#
# Multi-trait read_phenotype():
#   trait_col now accepts a character vector.  Single string -> returns a plain
#   numeric vector (backward compatible).  Character vector -> returns a named
#   numeric matrix (samples x traits) aligned to sample_order.
# ==============================================================================


# -- Internal: low-level chunked matrix builder --------------------------------

#' Fill a pre-allocated genotype matrix from a file in row chunks
#'
#' Two-pass strategy:
#'   Pass 1 -- header-only fread (nrows = 0) -> column names -> dimensions.
#'   Pass 2 -- successive fread(nrows = chunk_rows, skip = ...) -> fill matrix.
#'
#' Each chunk is released and gc(FALSE) is called before reading the next one
#' so that R's memory allocator can reclaim the chunk data.table immediately.
#'
#' @param file       Path to the delimited file.
#' @param meta_cols  Names of the left-hand metadata columns (excluded from the
#'   numeric block).
#' @param sep        Separator string passed to fread(). Default `"auto"`.
#' @param chunk_rows Rows per chunk. Default `50 000`.
#' @param na_strings Strings treated as NA. Default `c("NA", "")`.
#'
#' @return List with `meta` (data.table of metadata columns), `geno_mat`
#'   (pre-allocated numeric matrix, SNPs x samples), and `sample_ids`.
#' @keywords internal
#' @noRd
.read_dosage_chunked <- function(file, meta_cols,
                                  sep        = "auto",
                                  chunk_rows = 50000L,
                                  na_strings = c("NA", "")) {
  .assert_packages("data.table")

  # -- Pass 1: header scan -- zero data rows read -------------------------------
  hdr        <- data.table::fread(file, sep = sep, nrows = 0L,
                                   check.names = FALSE, data.table = TRUE)
  all_cols   <- names(hdr)
  sample_ids <- setdiff(all_cols, meta_cols)
  n_samples  <- length(sample_ids)
  if (n_samples == 0L)
    stop("No sample columns found after removing meta_cols.", call. = FALSE)

  # Total row count (read first column only -- minimal RAM)
  n_rows <- nrow(data.table::fread(file, sep = sep, select = 1L,
                                    header = TRUE, data.table = TRUE))

  # -- Pre-allocate the full output matrix -------------------------------------
  geno_mat           <- matrix(NA_real_, nrow = n_rows, ncol = n_samples)
  colnames(geno_mat) <- sample_ids

  n_chunks   <- ceiling(n_rows / chunk_rows)
  meta_list  <- vector("list", n_chunks)
  row_filled <- 0L
  chunk_idx  <- 0L

  # -- Pass 2: fill chunks -----------------------------------------------------
  repeat {
    if (row_filled >= n_rows) break
    chunk_idx <- chunk_idx + 1L

    chunk <- data.table::fread(
      file,
      sep        = sep,
      nrows      = chunk_rows,
      skip       = row_filled + 1L,   # +1 to jump over the header line
      header     = FALSE,
      col.names  = all_cols,
      check.names = FALSE,
      na.strings = na_strings,
      data.table = TRUE
    )
    if (nrow(chunk) == 0L) break

    r_start <- row_filled + 1L
    r_end   <- row_filled + nrow(chunk)

    # Numeric block -> matrix slice
    geno_block <- as.matrix(chunk[, sample_ids, with = FALSE])
    storage.mode(geno_block) <- "numeric"
    geno_mat[r_start:r_end, ] <- geno_block
    rm(geno_block)

    meta_list[[chunk_idx]] <- chunk[, meta_cols, with = FALSE]

    row_filled <- r_end
    rm(chunk)
    gc(FALSE)   # <- Inspiration 5: release allocator pressure between chunks
  }

  meta_dt <- data.table::rbindlist(meta_list[seq_len(chunk_idx)],
                                    use.names = TRUE)

  list(meta       = meta_dt,
       geno_mat   = geno_mat,
       sample_ids = sample_ids)
}


# -- Genotype readers -----------------------------------------------------------

#' Read numeric dosage genotype data
#'
#' Reads a file in the package's standard format:
#' `SNP, CHR, POS, REF, ALT, sample1, sample2, ...`.
#'
#' Files with more than `chunk_threshold` rows are read in chunks of
#' `chunk_rows` rows so that peak RAM equals one chunk rather than two full
#' copies of the file.  Smaller files use a single `fread()` call.
#'
#' @param file            Path to a delimited file (`.csv`, `.txt`, etc.).
#' @param sep             Field separator. `NULL` triggers auto-detection.
#' @param check_names     Sanitise column names. Default `FALSE`.
#' @param chunk_rows      Rows per chunk for the chunked path. Default `50 000`.
#' @param chunk_threshold Row count above which chunked reading is used.
#'   Default `200 000`.
#'
#' @return Named list: `snp_info`, `geno_mat`, `sample_ids`, `format`.
#' @export
read_numeric_genotype <- function(file, sep = NULL, check_names = FALSE,
                                   chunk_rows      = 50000L,
                                   chunk_threshold = 200000L) {
  .assert_packages("data.table")
  sep_val       <- sep %||% "auto"
  required_cols <- c("SNP", "CHR", "POS", "REF", "ALT")

  # Quick row probe (reads only the first column)
  n_rows_probe <- nrow(data.table::fread(file, sep = sep_val, select = 1L,
                                          header = TRUE, data.table = TRUE))

  if (n_rows_probe > chunk_threshold) {
    # -- Chunked path ----------------------------------------------------------
    res        <- .read_dosage_chunked(file, meta_cols = required_cols,
                                        sep = sep_val, chunk_rows = chunk_rows)
    meta_dt    <- res$meta
    geno_mat   <- res$geno_mat
    sample_ids <- res$sample_ids
    rownames(geno_mat) <- as.character(meta_dt$SNP)

    snp_info <- data.table::copy(meta_dt)
    snp_info[, SNP := as.character(SNP)]
    snp_info[, CHR := .normalise_chr(CHR)]
    snp_info[, POS := as.integer(POS)]
    snp_info[, REF := as.character(REF)]
    snp_info[, ALT := as.character(ALT)]

  } else {
    # -- Single-pass path ------------------------------------------------------
    dt   <- data.table::fread(file, sep = sep_val, check.names = check_names,
                               data.table = TRUE)
    miss <- setdiff(required_cols, names(dt))
    if (length(miss))
      stop("Numeric genotype file missing columns: ",
           paste(miss, collapse = ", "), call. = FALSE)

    sample_ids <- setdiff(names(dt), required_cols)
    geno_mat   <- as.matrix(dt[, sample_ids, with = FALSE])
    storage.mode(geno_mat) <- "numeric"
    rownames(geno_mat) <- as.character(dt$SNP)
    colnames(geno_mat) <- sample_ids

    snp_info <- data.table::copy(dt[, required_cols, with = FALSE])
    snp_info[, SNP := as.character(SNP)]
    snp_info[, CHR := .normalise_chr(CHR)]
    snp_info[, POS := as.integer(POS)]
    snp_info[, REF := as.character(REF)]
    snp_info[, ALT := as.character(ALT)]

    rm(dt); gc(FALSE)
  }

  .validate_genotype_object(snp_info, geno_mat)

  list(snp_info   = snp_info,
       geno_mat   = geno_mat,
       sample_ids = sample_ids,
       format     = "numeric")
}


#' Read HapMap genotype data
#'
#' Reads a HapMap file and converts nucleotide calls (`"AA"`, `"AT"`, ...) to
#' additive dosage (0/1/2).  Large files (> `chunk_threshold` rows) are read
#' in chunks; the dosage conversion is applied per chunk so peak RAM never
#' exceeds one chunk.
#'
#' @param file            Path to a HapMap text file.
#' @param chunk_rows      Rows per chunk for the chunked path. Default `50 000`.
#' @param chunk_threshold Row count above which chunked reading is used.
#'   Default `200 000`.
#'
#' @return Named list: `snp_info`, `geno_mat`, `sample_ids`, `format`.
#' @export
read_hapmap_genotype <- function(file,
                                  chunk_rows      = 50000L,
                                  chunk_threshold = 200000L) {
  .assert_packages("data.table")

  hmp_meta_std <- c("alleles", "chrom", "pos", "strand", "assembly#",
                    "center", "protLSID", "assayLSID", "panelLSID", "QCcode")

  # -- Header scan -------------------------------------------------------------
  hdr       <- data.table::fread(file, nrows = 0L, check.names = FALSE,
                                  data.table = TRUE)
  first_col <- names(hdr)[1L]
  req_hmp   <- c(first_col, "alleles", "chrom", "pos")
  miss      <- setdiff(req_hmp, names(hdr))
  if (length(miss))
    stop("HapMap file missing columns: ", paste(miss, collapse = ", "),
         call. = FALSE)

  meta_cols  <- intersect(c(first_col, hmp_meta_std), names(hdr))
  sample_ids <- setdiff(names(hdr), meta_cols)
  all_cols   <- names(hdr)
  n_samples  <- length(sample_ids)

  # -- Helper: convert one raw HapMap block -> dosage matrix -------------------
  # snp_info_block must have columns REF and ALT aligned to chunk rows.
  .block_to_dosage <- function(raw_chunk, ref_vec, alt_vec) {
    n_snp <- nrow(raw_chunk)
    mat   <- matrix(NA_real_, nrow = n_snp, ncol = n_samples,
                    dimnames = list(NULL, sample_ids))
    for (i in seq_len(n_snp)) {
      ri <- ref_vec[i]; ai <- alt_vec[i]
      gi <- as.character(unlist(raw_chunk[i, sample_ids, with = FALSE]))
      mat[i, ] <- vapply(gi, function(g) {
        if (is.na(g) || g %in% c("N", "NN", "--") || nchar(g) != 2L)
          return(NA_real_)
        a2 <- strsplit(g, "", fixed = TRUE)[[1L]]
        if (any(!a2 %in% c(ri, ai))) return(NA_real_)
        sum(a2 == ai)
      }, numeric(1L))
    }
    mat
  }

  n_rows_probe <- nrow(data.table::fread(file, select = 1L, header = TRUE,
                                          data.table = TRUE))

  if (n_rows_probe > chunk_threshold) {
    # -- Chunked path ----------------------------------------------------------
    geno_mat  <- matrix(NA_real_, nrow = n_rows_probe, ncol = n_samples,
                        dimnames = list(NULL, sample_ids))
    si_list   <- vector("list", ceiling(n_rows_probe / chunk_rows))
    row_filled <- 0L; chunk_idx <- 0L

    repeat {
      if (row_filled >= n_rows_probe) break
      chunk_idx <- chunk_idx + 1L

      chunk <- data.table::fread(
        file, nrows = chunk_rows, skip = row_filled + 1L,
        header = FALSE, col.names = all_cols,
        check.names = FALSE, data.table = TRUE
      )
      if (nrow(chunk) == 0L) break

      al    <- strsplit(as.character(chunk[["alleles"]]), "/", fixed = TRUE)
      ref_v <- vapply(al, function(x) x[1L] %||% NA_character_, character(1L))
      alt_v <- vapply(al, function(x) x[2L] %||% NA_character_, character(1L))

      r_start <- row_filled + 1L
      r_end   <- row_filled + nrow(chunk)
      geno_mat[r_start:r_end, ] <- .block_to_dosage(chunk, ref_v, alt_v)

      si_list[[chunk_idx]] <- data.table::data.table(
        SNP = as.character(chunk[[first_col]]),
        CHR = .normalise_chr(chunk[["chrom"]]),
        POS = as.integer(chunk[["pos"]]),
        REF = ref_v,
        ALT = alt_v
      )
      row_filled <- r_end
      rm(chunk, al, ref_v, alt_v); gc(FALSE)   # <- Inspiration 5
    }

    snp_info           <- data.table::rbindlist(si_list[seq_len(chunk_idx)])
    rownames(geno_mat) <- snp_info$SNP

  } else {
    # -- Single-pass path ------------------------------------------------------
    dt <- data.table::fread(file, check.names = FALSE, data.table = TRUE)
    al    <- strsplit(as.character(dt[["alleles"]]), "/", fixed = TRUE)
    ref_v <- vapply(al, function(x) x[1L] %||% NA_character_, character(1L))
    alt_v <- vapply(al, function(x) x[2L] %||% NA_character_, character(1L))

    snp_info <- data.table::data.table(
      SNP = as.character(dt[[first_col]]),
      CHR = .normalise_chr(dt[["chrom"]]),
      POS = as.integer(dt[["pos"]]),
      REF = ref_v,
      ALT = alt_v
    )
    geno_mat <- .block_to_dosage(dt, ref_v, alt_v)
    rownames(geno_mat) <- snp_info$SNP
    rm(dt, al, ref_v, alt_v); gc(FALSE)
  }

  .validate_genotype_object(snp_info, geno_mat)

  list(snp_info   = snp_info,
       geno_mat   = geno_mat,
       sample_ids = sample_ids,
       format     = "hapmap")
}


#' Read VCF genotype data
#'
#' Reads a VCF or bgzipped VCF using `VariantAnnotation` and converts GT
#' fields (phased or unphased) to additive dosage (0/1/2).
#'
#' @param file Path to a `.vcf` or `.vcf.gz` file.
#' @return Named list: `snp_info`, `geno_mat`, `sample_ids`, `format`.
#' @export
read_vcf_genotype <- function(file) {
  .assert_packages(c("VariantAnnotation", "SummarizedExperiment",
                     "GenomeInfoDb", "Biostrings", "S4Vectors"))

  vcf <- VariantAnnotation::readVcf(file)
  gt  <- VariantAnnotation::geno(vcf)$GT
  if (is.null(gt))
    stop("VCF does not contain GT genotype calls.", call. = FALSE)

  gt_to_dosage <- function(x) {
    x <- gsub("\\|", "/", x)
    vapply(strsplit(x, "/", fixed = TRUE), function(a) {
      if (length(a) != 2L || any(a %in% c(".", NA_character_))) return(NA_real_)
      aa <- suppressWarnings(as.integer(a))
      if (any(is.na(aa))) return(NA_real_)
      sum(aa)
    }, numeric(1L))
  }

  geno_mat <- apply(gt, 2L, gt_to_dosage)
  if (!is.matrix(geno_mat)) geno_mat <- matrix(geno_mat, ncol = ncol(gt))

  rr      <- SummarizedExperiment::rowRanges(vcf)
  alt_vec <- vapply(VariantAnnotation::alt(vcf),
                    function(x) as.character(x[[1L]]), character(1L))
  snp_ids <- names(rr)
  empty   <- is.na(snp_ids) | snp_ids == "" | snp_ids == "."
  if (any(empty))
    snp_ids[empty] <- paste0(as.character(GenomeInfoDb::seqnames(rr)[empty]),
                              ":", S4Vectors::start(rr)[empty])

  snp_info <- data.table::data.table(
    SNP = as.character(snp_ids),
    CHR = .normalise_chr(as.character(GenomeInfoDb::seqnames(rr))),
    POS = as.integer(S4Vectors::start(rr)),
    REF = as.character(VariantAnnotation::ref(vcf)),
    ALT = alt_vec
  )
  rownames(geno_mat) <- snp_info$SNP
  colnames(geno_mat) <- colnames(gt)
  storage.mode(geno_mat) <- "numeric"

  .validate_genotype_object(snp_info, geno_mat)

  list(snp_info   = snp_info,
       geno_mat   = geno_mat,
       sample_ids = colnames(geno_mat),
       format     = "vcf")
}


#' Dispatch genotype reading by format
#'
#' @param file   Path to the genotype file.
#' @param format One of `"auto"` (default), `"numeric"`, `"hapmap"`, `"vcf"`.
#' @return Named list from the format-specific reader.
#' @export
read_genotype <- function(file,
                          format = c("auto", "numeric", "hapmap", "vcf")) {
  format <- match.arg(format)
  if (identical(format, "auto")) {
    lower  <- tolower(file)
    format <- if (grepl("\\.vcf(\\.gz)?$", lower)) "vcf"
              else if (grepl("\\.hmp(\\.txt)?$", lower) ||
                       grepl("hapmap", lower))          "hapmap"
              else                                       "numeric"
  }
  switch(format,
    numeric = read_numeric_genotype(file),
    hapmap  = read_hapmap_genotype(file),
    vcf     = read_vcf_genotype(file)
  )
}


# -- Phenotype reader (multi-trait) ---------------------------------------------

#' Read phenotype and covariate data (supports single or multiple traits)
#'
#' Reads a phenotype file, extracts one or more trait columns plus optional
#' covariates, and aligns rows to the genotype sample order.  All sample
#' matching is handled internally -- the caller never touches row alignment.
#'
#' @section Multi-trait behaviour:
#' * **Single trait** (`trait_col` is a length-1 string) -- returns
#'   `$phenotype` as a plain numeric vector.  Backward compatible with all
#'   existing single-trait code.
#' * **Multiple traits** (`trait_col` is a character vector of length > 1) --
#'   returns `$phenotype` as a numeric matrix with dimensions
#'   (n_samples x n_traits), column names equal to `trait_col`, row names
#'   equal to the aligned sample IDs.
#'
#' @param file         Path to a delimited text file (CSV, TSV, ...).
#' @param sample_col   Column holding sample identifiers. Default `"Sample"`.
#' @param trait_col    Trait column name(s).  A single string or a character
#'   vector for multiple traits.  Default `"Trait1"`.
#' @param covar_cols   Character vector of covariate column names.
#'   `NULL` (default) means no covariates.  Covariates are shared across
#'   all traits.
#' @param sample_order Character vector giving the target sample order -- i.e.,
#'   `g$sample_ids` from `read_genotype()`.  Rows are reordered to match.
#'   `NULL` skips reordering.
#' @param sep          Separator.  `NULL` triggers auto-detection.
#'
#' @return A named list:
#'   \describe{
#'     \item{`phenotype`}{Numeric vector (single trait) **or** numeric matrix
#'       (n_samples x n_traits) for multiple traits.}
#'     \item{`trait_names`}{Character vector of returned trait name(s).}
#'     \item{`covariates`}{Shared covariate `data.frame`, or `NULL`.}
#'     \item{`sample_ids`}{Sample IDs in the aligned order.}
#'   }
#' @export
#' @examples
#' \dontrun{
#' geno <- read_genotype("geno.csv")
#'
#' # Single trait
#' p1 <- read_phenotype("pheno.csv", sample_order = geno$sample_ids)
#' stopifnot(is.numeric(p1$phenotype) && !is.matrix(p1$phenotype))
#'
#' # Multiple traits with shared PC covariates
#' p2 <- read_phenotype(
#'   "pheno.csv",
#'   trait_col    = c("BlastLeaf", "BlastPanicle", "Yield"),
#'   covar_cols   = c("PC1", "PC2"),
#'   sample_order = geno$sample_ids
#' )
#' stopifnot(is.matrix(p2$phenotype))   # n_samples x 3
#' }
read_phenotype <- function(file,
                            sample_col   = "Sample",
                            trait_col    = "Trait1",
                            covar_cols   = NULL,
                            sample_order = NULL,
                            sep          = NULL) {
  .assert_packages("data.table")
  trait_col <- as.character(trait_col)   # guard against factor input

  dt <- data.table::fread(file, sep = sep %||% "auto", data.table = TRUE)

  # -- Column validation --------------------------------------------------------
  need <- unique(c(sample_col, trait_col, covar_cols))
  miss <- setdiff(need, names(dt))
  if (length(miss))
    stop("Phenotype file missing column(s): ",
         paste(miss, collapse = ", "), call. = FALSE)

  dt[[sample_col]] <- as.character(dt[[sample_col]])
  if (anyDuplicated(dt[[sample_col]]))
    stop("Duplicate sample IDs in column '", sample_col, "'.", call. = FALSE)

  # Coerce all trait columns to numeric; warn on silent NA introduction
  for (tc in trait_col) {
    raw         <- dt[[tc]]
    dt[[tc]]    <- suppressWarnings(as.numeric(raw))
    n_lost      <- sum(is.na(dt[[tc]]) & !is.na(raw))
    if (n_lost > 0L)
      warning("Trait '", tc, "': ", n_lost, " non-numeric value(s) -> NA.",
              call. = FALSE)
  }

  # -- Sample alignment ---------------------------------------------------------
  if (!is.null(sample_order)) {
    sample_order <- as.character(sample_order)

    missing_in_pheno <- setdiff(sample_order, dt[[sample_col]])
    if (length(missing_in_pheno))
      stop(length(missing_in_pheno),
           " genotype sample(s) absent from phenotype file.\n",
           "  First missing: ",
           paste(head(missing_in_pheno, 5L), collapse = ", "), "\n",
           "  Check that '", sample_col,
           "' values match genotype column names exactly.",
           call. = FALSE)

    extra <- setdiff(dt[[sample_col]], sample_order)
    if (length(extra))
      message("  [phenotype] ", length(extra),
              " sample(s) not in genotype -- dropped.")

    dt <- dt[match(sample_order, dt[[sample_col]]), ]
  }

  # -- Build phenotype return value ---------------------------------------------
  # Single trait  -> plain numeric vector   (backward compatible)
  # Multiple traits -> numeric matrix       (n_samples x n_traits)
  pheno_out <- if (length(trait_col) == 1L) {
    as.numeric(dt[[trait_col]])
  } else {
    m <- as.matrix(dt[, trait_col, with = FALSE])
    storage.mode(m) <- "numeric"
    rownames(m)     <- as.character(dt[[sample_col]])
    m
  }

  list(
    phenotype   = pheno_out,
    trait_names = trait_col,
    covariates  = if (!is.null(covar_cols)) as.data.frame(dt[, covar_cols, with = FALSE])
                  else NULL,
    sample_ids  = as.character(dt[[sample_col]])
  )
}


# -- Internal helpers -----------------------------------------------------------

#' Strip leading "chr" prefix from chromosome names (case-insensitive)
#' @keywords internal
#' @noRd
.normalise_chr <- function(x) sub("^chr", "", as.character(x), ignore.case = TRUE)


#' Hard-validate sample alignment after read_phenotype()
#'
#' Called immediately after read_phenotype() in run_sldp() as a second line
#' of defence.  Stops with an informative message if the vectors differ.
#' @keywords internal
#' @noRd
.harmonise_samples <- function(geno_ids, pheno_ids) {
  missing <- setdiff(geno_ids, pheno_ids)
  if (length(missing))
    stop(length(missing), " genotype sample(s) missing from phenotype: ",
         paste(head(missing, 5L), collapse = ", "), call. = FALSE)
  if (!identical(as.character(geno_ids), as.character(pheno_ids)))
    stop("Sample order mismatch after alignment. ",
         "read_phenotype() must be called with sample_order = g$sample_ids.",
         call. = FALSE)
  invisible(TRUE)
}


# -- Writers --------------------------------------------------------------------

#' Write numeric dosage genotype file
#'
#' @param snp_info A data.frame or data.table of SNP metadata with columns
#'   `SNP`, `CHR`, `POS`, `REF`, `ALT`.
#' @param geno_mat Numeric matrix of dosage values (0/1/2/NA), SNPs x samples.
#' @param file Output file path.
#' @param sep Field separator. Default `","`.
#' @return Invisibly returns `file`.
#' @export
write_numeric_genotype <- function(snp_info, geno_mat, file, sep = ",") {
  .assert_packages("data.table")
  .validate_genotype_object(snp_info, geno_mat)
  dt <- data.table::as.data.table(
    cbind(snp_info, as.data.frame(geno_mat, check.names = FALSE))
  )
  data.table::fwrite(dt, file = file, sep = sep, quote = FALSE, na = "NA")
  invisible(file)
}


#' Write HapMap genotype file
#'
#' Converts numeric dosage values (0/1/2/NA) to nucleotide calls
#' (`"AA"`, `"AT"`, `"TT"`, `"NN"`) and writes a tab-delimited HapMap file.
#'
#' @param snp_info A data.frame or data.table of SNP metadata with columns
#'   `SNP`, `CHR`, `POS`, `REF`, `ALT`.
#' @param geno_mat Numeric matrix of dosage values (0/1/2/NA), SNPs x samples.
#' @param file Output file path.
#' @return Invisibly returns `file`.
#' @export
write_hapmap_genotype <- function(snp_info, geno_mat, file) {
  .assert_packages("data.table")
  .validate_genotype_object(snp_info, geno_mat)

  geno_chr <- matrix(NA_character_, nrow = nrow(geno_mat), ncol = ncol(geno_mat),
                     dimnames = dimnames(geno_mat))
  for (i in seq_len(nrow(geno_mat))) {
    r <- snp_info$REF[i]; a <- snp_info$ALT[i]; gi <- geno_mat[i, ]
    geno_chr[i, gi == 0]   <- paste0(r, r)
    geno_chr[i, gi == 1]   <- paste0(r, a)
    geno_chr[i, gi == 2]   <- paste0(a, a)
    geno_chr[i, is.na(gi)] <- "NN"
  }
  out <- data.frame(
    check.names = FALSE,
    "rs#"       = snp_info$SNP,
    alleles     = paste0(snp_info$REF, "/", snp_info$ALT),
    chrom       = snp_info$CHR,
    pos         = snp_info$POS,
    strand      = "+",
    "assembly#" = NA_character_, center    = NA_character_,
    protLSID    = NA_character_, assayLSID = NA_character_,
    panelLSID   = NA_character_, QCcode    = NA_character_,
    as.data.frame(geno_chr, check.names = FALSE)
  )
  data.table::fwrite(out, file = file, sep = "\t", quote = FALSE, na = "NN")
  invisible(file)
}


#' Write pruning statistics report
#'
#' @param pruning_stats A data.frame or data.table of step-by-step pruning
#'   counts as returned in the `pruning_stats` element of `run_sldp()`.
#' @param stats_file Path for the CSV output file.
#' @param summary_file Optional path for a plain-text summary file. If `NULL`
#'   (default) no summary is written.
#' @return Invisibly returns a named list with elements `stats_file` and
#'   `summary_file`.
#' @export
write_pruning_report <- function(pruning_stats, stats_file,
                                  summary_file = NULL) {
  .assert_packages("data.table")
  data.table::fwrite(pruning_stats, stats_file, sep = ",", quote = FALSE)
  if (!is.null(summary_file)) {
    lines <- apply(as.data.frame(pruning_stats), 1L, function(x) {
      d <- if (nzchar(x[["details"]])) paste0(" | ", x[["details"]]) else ""
      sprintf("  - %-30s  %d -> %d  (removed %d)%s",
              x[["step"]], as.integer(x[["n_before"]]),
              as.integer(x[["n_after"]]), as.integer(x[["n_removed"]]), d)
    })
    writeLines(c("SLDP Pruning Summary", "====================",
                 sprintf("Steps recorded: %d", nrow(pruning_stats)), "", lines),
               con = summary_file)
  }
  invisible(list(stats_file = stats_file, summary_file = summary_file))
}
