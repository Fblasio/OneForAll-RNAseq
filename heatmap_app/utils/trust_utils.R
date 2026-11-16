# =============================================================================
# trust_utils.R
# -----------------------------------------------------------------------------
# A single convention for empty and failed results, so a blank panel always
# states WHY it is blank. Without it the user cannot distinguish "no signal in
# your data" from "the analysis never ran" -- two conclusions that call for
# opposite next steps.
#
# Use empty_result() to build the standard {table, object} shape that the GO/
# KEGG modules already expect, but with a *reason* the UI can render.
#
# Pattern in an enrichment function:
#   if (nothing_mapped)        return(empty_result("KEGG needs Entrez IDs; 0 of N mapped."))
#   if (nrow(res@result)==0)   return(empty_result("No pathways below the chosen cutoffs."))
#
# Pattern in a heatmap/reactive that returns NULL:
#   validate(need(nrow(sub) > 0, explain_empty("heatmap", n_input, n_matched)))
# =============================================================================

# Standard empty payload compatible with the GO/KEGG {table, object} contract.
empty_result <- function(reason) {
  list(table  = data.frame(Message = reason, stringsAsFactors = FALSE),
       object = NULL,
       reason = reason)
}

# Build a specific, numeric explanation for an empty gene-matching step.
# stage: short label ("GO heatmap", "cart heatmap", "KEGG").
explain_empty <- function(stage, n_input = NA, n_matched = NA, extra = NULL) {
  parts <- sprintf("No genes to show in the %s.", stage)
  if (!is.na(n_input) && !is.na(n_matched)) {
    parts <- paste0(parts, sprintf(
      " %d gene(s) selected, %d matched the current matrix.", n_input, n_matched))
  }
  if (!is.null(extra)) parts <- paste0(parts, " ", extra)
  parts
}

# Write a result table to CSV with that run's caveats as leading '#' lines.
#
# Every enrichment run can produce warnings that change how its table must be
# read -- simplify() adjusting p over MORE terms than are listed, KEGG dropping
# query genes that have no Entrez ID, a background that excluded part of the
# query. They were shown as transient notifications and then lost, so the
# exported CSV (the file that ends up in a supplement) carried a p.adjust column
# with nothing to say it was computed over a different set of rows.
#
# '#' lines are skipped by read.csv(comment.char = "#") and by most importers,
# and Excel shows them as leading rows. The table itself is written byte-for-
# byte as before, so anything already parsing these files keeps working.
#
# Embedded newlines are split so a multi-line warning cannot break out of the
# comment block and be parsed as data.
write_result_csv <- function(result, path, extra = character(0)) {
  notes <- c(if (!is.null(result$warnings)) as.character(result$warnings),
             if (!is.null(result$reason))   as.character(result$reason),
             as.character(extra))
  notes <- notes[!is.na(notes) & nzchar(notes)]
  con <- file(path, open = "wt")
  on.exit(close(con), add = TRUE)
  if (length(notes) > 0) {
    writeLines(c("# ---- Caveats recorded for this run ----",
                 paste0("# ", unlist(strsplit(notes, "\n", fixed = TRUE))),
                 "# ---------------------------------------"), con)
  }
  utils::write.csv(result$table, con, row.names = FALSE)
  invisible(length(notes))
}

# Diagnose an ID-mapping step and return a message if it failed, else NULL.
# Used before enrichment to distinguish "no enrichment" from "IDs didn't map".
# mapped: character vector of successfully mapped IDs (e.g. Entrez).
# input : the original gene vector.
mapping_diagnostic <- function(input, mapped, id_kind = "Entrez") {
  n_in  <- length(unique(canon_ids(input)))
  n_map <- length(unique(mapped[!is.na(mapped) & nzchar(mapped)]))
  if (n_map == 0) {
    return(sprintf(
      "None of the %d input gene IDs could be mapped to %s IDs for this species; enrichment skipped.",
      n_in, id_kind))
  }
  if (n_map < 0.25 * n_in) {
    return(sprintf(
      "Only %d of %d gene IDs mapped to %s; results may be unreliable.",
      n_map, n_in, id_kind))
  }
  NULL
}
