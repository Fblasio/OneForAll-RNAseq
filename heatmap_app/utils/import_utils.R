# =============================================================================
# import_utils.R
# Utility functions for importing count data into the RNA-seq pipeline.
# Two branches are supported:
#   - Raw counts: CSV or Excel matrix (genes x samples)
#   - Salmon output: quant.sf files aggregated via tximport
# =============================================================================

# -----------------------------------------------------------------------------
# import_raw_counts()
# Reads a CSV or Excel counts matrix.
# Returns a named list:
#   $type         — "raw"
#   $counts       — data.frame, genes as rows, samples as columns
#   $sample_names — character vector of column names
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# read_counts_file()
# Reads a CSV or Excel counts table and returns the gene IDs SEPARATELY from the
# value columns:
#   $ids    — character vector, one per row, exactly as written in column 1
#   $values — data.frame of the remaining columns, with NO row names set
#   $notes  — human-readable observations for the caller to surface
#
# Keeping the IDs out of the row names is the whole point. A data.frame cannot
# hold duplicate row names, and neither `read.csv(row.names = 1)` nor
# `rownames(df) <- df[[1]]` will accept them -- both raise
# "duplicate 'row.names' are not allowed" while READING the file, before any
# validation runs. A repeated gene ID is ordinary in a per-transcript export,
# and the app's stated contract is to aggregate such rows, so the read step must
# tolerate what the aggregation step is designed to handle.
#
# Row order and row count are preserved exactly, so $ids stays aligned with
# $values through validate_count_matrix() (which only ever drops COLUMNS).
# -----------------------------------------------------------------------------
read_counts_file <- function(file_path, file_name) {
  ext <- tolower(tools::file_ext(file_name))

  raw <- if (ext == "csv") {
    utils::read.csv(file_path, check.names = FALSE, stringsAsFactors = FALSE)
  } else {
    as.data.frame(readxl::read_excel(file_path), stringsAsFactors = FALSE)
  }

  if (is.null(raw) || ncol(raw) < 2) {
    stop("The counts file needs a gene-ID column plus at least one sample column.")
  }

  ids    <- trimws(as.character(raw[[1]]))
  # An empty first-column cell arrives as NA from readxl, and so does the literal
  # text "NA" from either reader. NA cannot be a row name -- the next thing every
  # caller does is rownames(x) <- make.unique(ids), which stops with
  # "missing values in 'row.names' are not allowed". That error is not caught
  # anywhere, so in the app it took down mat_data() and with it the whole
  # sidebar, naming R internals rather than the offending row.
  # Mapped to "" so the blank-ID note below is TRUE: the row really is kept, and
  # collapse_counts_by_id() aggregates all such rows together under one empty ID.
  ids[is.na(ids)] <- ""
  values <- raw[, -1, drop = FALSE]
  rownames(values) <- NULL

  notes <- character(0)
  n_dup <- sum(duplicated(ids))
  if (n_dup > 0) {
    notes <- c(notes, sprintf(
      paste0("%d row(s) repeat a gene ID already seen (e.g. %s). Rows sharing an ID ",
             "are summed into one, which is the correct handling for a per-transcript ",
             "export but not for a file that was meant to have unique IDs -- check ",
             "which of the two you have."),
      n_dup, paste(utils::head(unique(ids[duplicated(ids)]), 3), collapse = ", ")))
  }
  n_blank <- sum(is.na(ids) | !nzchar(ids))
  if (n_blank > 0) {
    notes <- c(notes, sprintf(
      "%d row(s) have a blank gene ID. They are kept and aggregated together under an empty ID.",
      n_blank))
  }

  list(ids = ids, values = values, notes = notes)
}

# -----------------------------------------------------------------------------
# collapse_counts_by_id()
# Canonicalizes `ids` for the active species and sums rows that share a
# canonical ID. `values` must be the numeric matrix/data.frame whose rows are
# still aligned 1:1 with `ids`.
#
# rowsum() is the base R function documented for group-wise row aggregation and
# guarantees unique output row names, which is what makes the result a legal
# data.frame again.
# -----------------------------------------------------------------------------
collapse_counts_by_id <- function(values, ids) {
  canon <- canon_ids(ids)
  m <- as.matrix(values)
  if (anyDuplicated(canon) > 0) {
    out <- rowsum(m, group = canon, reorder = FALSE)
  } else {
    out <- m
    rownames(out) <- canon
  }
  as.data.frame(out, check.names = FALSE)
}

import_raw_counts <- function(file_path, file_name) {
  imported <- read_counts_file(file_path, file_name)
  ids      <- imported$ids
  counts   <- imported$values

  if (ncol(counts) == 0) {
    stop("The counts file has no sample columns after reading the gene-ID column.")
  }

  # Row names are needed for validate_count_matrix()'s own ID checks, but a
  # data.frame cannot carry duplicates -- so the ORIGINAL ids are kept in `ids`
  # and only a uniquified copy is attached here. The aggregation below uses
  # `ids`, never these placeholder names.
  rownames(counts) <- make.unique(ids, sep = "..dup")

  v <- validate_count_matrix(counts, require_counts = TRUE, source_label = "counts file")
  if (!v$ok) stop(paste(v$errors, collapse = " "))
  # validate_count_matrix() drops non-numeric COLUMNS but never rows, so `ids`
  # is still aligned with v$matrix row for row.
  counts <- as.data.frame(v$matrix, check.names = FALSE)

  # Canonicalize gene IDs at the earliest possible point, before the matrix is
  # consumed by any downstream reactive. Rows that collapse onto one canonical
  # ID -- whether they were duplicated in the file or made equal by dropping a
  # tomato version suffix -- are aggregated by summation.
  counts <- collapse_counts_by_id(counts, ids)

  attr(counts, "validation_warnings") <- c(imported$notes, v$warnings)

  list(
    type = "raw",
    counts = counts,
    sample_names = colnames(counts)
  )
}

# -----------------------------------------------------------------------------
# parse_salmon_uploads()
# Converts the data.frame produced by Shiny fileInput(multiple = TRUE) into
# a named character vector suitable for tximport:
#   names  = sample identifiers extracted from file names
#   values = temporary file paths on the server
#
# Name extraction rule:
#   "WT_1_quant.sf"  → "WT_1"
#   "mut_2_quant.sf" → "mut_2"
#   "sample_A.sf"    → "sample_A"
# The function removes the extension and, if the base name ends with
# "_quant", also removes that suffix.
# -----------------------------------------------------------------------------
parse_salmon_uploads <- function(file_input_df) {
  if (is.null(file_input_df) || nrow(file_input_df) == 0) {
    stop("No quant.sf files were provided.")
  }
  
  base_names <- tools::file_path_sans_ext(file_input_df$name)
  sample_ids <- sub("_quant$", "", base_names)
  
  if (any(duplicated(sample_ids))) {
    stop(paste0(
      "Duplicate sample names detected after parsing file names: ",
      paste(sample_ids[duplicated(sample_ids)], collapse = ", "),
      ". Rename the quant.sf files so each encodes a unique sample ID."
    ))
  }
  
  stats::setNames(file_input_df$datapath, sample_ids)
}

# -----------------------------------------------------------------------------
# import_salmon_counts()
# Aggregates transcript-level Salmon output to gene level using tximport.
# countsFromAbundance = "lengthScaledTPM" is used so that the resulting
# count estimates are on the scale of raw counts and suitable for
# DESeqDataSetFromTximport() without additional transformation.
# Reference: Love et al. (2016) F1000Research 4:1521
#            https://bioconductor.org/packages/release/bioc/vignettes/
#            tximport/inst/doc/tximport.html#DESeq2
#
# Arguments:
#   quant_paths  — named character vector: sample_id = path to quant.sf
#   tx2gene_path — path to a two-column CSV/TSV (transcript_id, gene_id)
#   quant_type   — "salmon" (default) or "alevin"
#   tx2gene_sep  — "," for CSV, "\t" for TSV (auto-detected by the module)
#
# Returns a named list:
#   $type         — "salmon"
#   $txi          — tximport result list (counts, abundance, length)
#   $sample_names — character vector of sample names (= names(quant_paths))
# -----------------------------------------------------------------------------
import_salmon_counts <- function(quant_paths,
                                 tx2gene_path,
                                 quant_type  = "salmon",
                                 tx2gene_sep = "\t") {
  
  if (!requireNamespace("tximport", quietly = TRUE)) {
    stop(
      "The tximport package is required for Salmon import. ",
      "Install it with: BiocManager::install('tximport')"
    )
  }
  
  missing_files <- quant_paths[!file.exists(quant_paths)]
  if (length(missing_files) > 0) {
    stop(paste0(
      "The following quant.sf files could not be found: ",
      paste(names(missing_files), collapse = ", ")
    ))
  }
  
  tx2gene <- tryCatch(
    read.table(tx2gene_path, header = FALSE, sep = tx2gene_sep,
               stringsAsFactors = FALSE, fill = TRUE,
               quote = "", comment.char = ""),
    error = function(e) {
      stop("Failed to read the tx2gene table: ", conditionMessage(e))
    }
  )
  
  if (ncol(tx2gene) < 2) {
    stop("The tx2gene table must have at least two columns: transcript ID and gene ID.")
  }
  
  # Keep only the first two columns and name them explicitly.
  tx2gene <- tx2gene[, 1:2, drop = FALSE]
  colnames(tx2gene) <- c("tx", "gene")
  
  # Canonicalize tomato GENE IDs BEFORE tximport so that tximport performs the
  # gene-level aggregation itself. tximport computes per-gene $length as an
  # abundance-weighted average (a valid DESeq2 offset); summing $length after
  # the fact — as the old code did — corrupts normalization.
  tx2gene$gene <- canon_ids(tx2gene$gene)
  
  txi <- tximport::tximport(
    quant_paths,
    type = quant_type,
    tx2gene = tx2gene,
    countsFromAbundance = "lengthScaledTPM",
    ignoreTxVersion = TRUE
  )
  # No post-hoc collapse needed: rownames(txi$counts) are already canonical
  # gene IDs and tximport aggregated counts/abundance/length correctly.
  
  list(
    type = "salmon",
    txi = txi,
    sample_names = names(quant_paths)
  )
}