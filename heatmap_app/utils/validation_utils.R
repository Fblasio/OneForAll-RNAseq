# =============================================================================
# validation_utils.R
# -----------------------------------------------------------------------------
# Input validation for uploaded count/expression matrices.
#
# One function inspects the matrix and returns a structured verdict with
# human-readable messages, so a malformed file is refused at the door with a
# specific reason rather than failing somewhere deep inside DESeq2, scale() or
# rowsum() with a message that names none of the user's own concepts.
#
# Returns a list:
#   $ok        logical  — TRUE if the matrix is usable (possibly after coercion)
#   $errors    character — blocking problems (ok = FALSE)
#   $warnings  character — non-blocking problems the user should see
#   $matrix    matrix    — cleaned numeric matrix (NULL if $ok is FALSE)
#   $dropped_cols character — non-numeric columns that were removed
#
# Design goal (the "trust gap"): every failure names WHAT is wrong and WHERE.
# =============================================================================

# Soft cap constants; keep here so both validation and memory checks agree.
# (MAX_GENES_MEMORY may already be defined in global.R; we only default it.)
if (!exists("MAX_GENES_MEMORY")) MAX_GENES_MEMORY <- 60000L
ONEFORALL_MAX_CELLS <- 40e6   # genes * samples above which we warn about size

validate_count_matrix <- function(df,
                                   require_counts = FALSE,
                                   source_label = "uploaded file") {
  errors   <- character(0)
  warnings <- character(0)
  dropped  <- character(0)

  # --- Shape --------------------------------------------------------------
  if (is.null(df) || nrow(df) == 0 || ncol(df) == 0) {
    return(list(ok = FALSE,
                errors = sprintf("The %s is empty or has no data columns.", source_label),
                warnings = warnings, matrix = NULL, dropped_cols = dropped))
  }
  if (ncol(df) < 2) {
    warnings <- c(warnings, sprintf(
      "The %s has only 1 sample column. Clustering and z-score need >= 2.",
      source_label))
  }
  if (nrow(df) < 2) {
    errors <- c(errors, sprintf("The %s has only 1 gene row; nothing to compare.",
                                source_label))
  }

  # --- Row names (gene IDs) ----------------------------------------------
  ids <- rownames(df)
  if (is.null(ids)) {
    errors <- c(errors, paste0(
      "No gene IDs found: the first column must contain gene identifiers ",
      "(used as row names)."))
  } else {
    blank <- is.na(ids) | !nzchar(trimws(ids))
    if (any(blank)) {
      warnings <- c(warnings, sprintf(
        "%d gene ID(s) are blank/NA and those rows will be unusable.", sum(blank)))
    }
    dup <- duplicated(ids)
    if (any(dup)) {
      warnings <- c(warnings, sprintf(
        "%d duplicate gene ID(s) detected (e.g. %s). They will be aggregated.",
        sum(dup), paste(utils::head(unique(ids[dup]), 3), collapse = ", ")))
    }
  }

  # --- Column types: drop non-numeric, coerce the rest --------------------
  is_num <- vapply(df, function(col) {
    is.numeric(col) || (is.character(col) &&
      suppressWarnings(!all(is.na(as.numeric(col[nzchar(col) & !is.na(col)])))))
  }, logical(1))

  # A column counts as non-numeric if coercion yields all-NA on non-empty values.
  hard_nonnum <- vapply(seq_along(df), function(j) {
    col <- df[[j]]
    if (is.numeric(col)) return(FALSE)
    vals <- suppressWarnings(as.numeric(as.character(col)))
    all(is.na(vals))
  }, logical(1))

  if (all(hard_nonnum)) {
    errors <- c(errors, paste0(
      "None of the data columns are numeric. Check that the first column is ",
      "gene IDs and the rest are expression values (no text headers inside data)."))
    return(list(ok = FALSE, errors = errors, warnings = warnings,
                matrix = NULL, dropped_cols = colnames(df)[hard_nonnum]))
  }
  if (any(hard_nonnum)) {
    dropped <- colnames(df)[hard_nonnum]
    warnings <- c(warnings, sprintf(
      "%d non-numeric column(s) removed: %s.",
      length(dropped), paste(dropped, collapse = ", ")))
    df <- df[, !hard_nonnum, drop = FALSE]
  }

  # --- Column names (sample IDs) ------------------------------------------
  # BLOCKING, unlike duplicate gene IDs. Group means are computed BY NAME
  # (rowMeans(expr_mat[, sel, drop = FALSE]) in app.r), and R resolves every
  # occurrence of a duplicated column name to the FIRST match. Two columns both
  # called "WT_1" therefore average the first one with itself and discard the
  # second replicate -- silently, from a file that otherwise validates cleanly.
  # These cannot be aggregated the way duplicate gene rows are: nothing here can
  # know whether they are technical replicates to merge or a labelling mistake,
  # so refuse and name them. Checked AFTER the non-numeric drop above, so text
  # columns that were removed anyway cannot trigger a false rejection. Compared
  # on trimmed names, so "WT_1" and "WT_1 " are caught as the collision they are.
  cn <- trimws(colnames(df))
  if (anyDuplicated(cn) > 0) {
    dup_cn <- unique(cn[duplicated(cn)])
    errors <- c(errors, sprintf(
      paste0("%d duplicate sample column name(s) in the %s: %s. Sample names are ",
             "used as column IDs, so a duplicate silently hides one replicate. ",
             "Please make them unique."),
      length(dup_cn), source_label, paste(dup_cn, collapse = ", ")))
  }

  mat <- as.matrix(df)
  storage.mode(mat) <- "numeric"

  # --- NA / negative / count sanity --------------------------------------
  n_na <- sum(is.na(mat))
  if (n_na > 0) {
    warnings <- c(warnings, sprintf(
      "%d missing value(s) (NA) present. They are treated as NA by rowMeans and may drop genes.",
      n_na))
  }
  if (any(mat < 0, na.rm = TRUE)) {
    warnings <- c(warnings, paste0(
      "Negative values present. For 'raw counts' these are invalid for ",
      "DESeq2/VST and will be clamped to 0 before building the model."))
  }
  if (require_counts) {
    frac_nonint <- mean(abs(mat - round(mat)) > 1e-8, na.rm = TRUE)
    if (isTRUE(frac_nonint > 0.01)) {
      warnings <- c(warnings, sprintf(
        "%.0f%% of values are non-integer but 'raw counts' is selected; they will be rounded.",
        100 * frac_nonint))
    }
    # Integer-ness does not separate counts from normalized data: RSEM
    # expected_count is non-integer, and CPM or normalized counts rounded
    # upstream are integer. Library size does separate them.
    ls <- colSums(mat, na.rm = TRUE)
    ls <- ls[is.finite(ls) & ls > 0]
    if (length(ls) >= 2) {
      ls_ratio <- max(ls) / min(ls)
      if (all(abs(ls - 1e6) / 1e6 < 0.01)) {
        errors <- c(errors, paste0(
          "Every column sums to ~1,000,000. This is CPM or TPM, not raw counts. ",
          "DESeq2 models counting noise and requires un-normalized counts; it ",
          "estimates size factors itself. Untick 'These are raw counts', or supply ",
          "the original count matrix."))
      } else if (ls_ratio < 1.05) {
        warnings <- c(warnings, sprintf(paste0(
          "Library sizes differ by only %.3fx across samples. Raw counts rarely look ",
          "like this; the matrix has probably already been normalized (DESeq2 ",
          "normalized counts, CPM, or depth-equalized). If so, the estimated size ",
          "factors will all be ~1 and the dispersion estimates will be wrong."),
          ls_ratio))
      }
    }
  }

  # --- Size guard ---------------------------------------------------------
  cells <- as.numeric(nrow(mat)) * as.numeric(ncol(mat))
  if (nrow(mat) > MAX_GENES_MEMORY) {
    warnings <- c(warnings, sprintf(
      "%d genes exceeds the recommended %d; consider pre-filtering for speed.",
      nrow(mat), MAX_GENES_MEMORY))
  }
  if (cells > ONEFORALL_MAX_CELLS) {
    warnings <- c(warnings, sprintf(
      "Large matrix (%s cells). Clustering may be slow or memory-heavy.",
      format(cells, big.mark = ",", scientific = FALSE)))
  }

  ok <- length(errors) == 0
  list(ok = ok,
       errors = errors,
       warnings = warnings,
       matrix = if (ok) mat else NULL,
       dropped_cols = dropped)
}

# Convenience: turn a validation result into Shiny notifications.
# Call inside an observer/reactive. Returns invisibly TRUE if ok.
notify_validation <- function(v, session = shiny::getDefaultReactiveDomain()) {
  for (w in v$warnings)
    shiny::showNotification(paste0("\u26a0\ufe0f ", w), type = "warning", duration = 9)
  for (e in v$errors)
    shiny::showNotification(paste0("\u274c ", e), type = "error", duration = 12)
  invisible(isTRUE(v$ok))
}
