# =============================================================================
# consistency_utils.R
# -----------------------------------------------------------------------------
# Gene-ID consistency: one implementation each of the three operations every
# module needs, so no module has to reimplement them and none can disagree.
#   (a) canonicalize a gene ID              -> canon_ids()
#   (b) compare/match two ID vectors        -> ids_match() / ids_in()
#   (c) print a debug line                  -> dbg()   (silent unless enabled)
# =============================================================================

# --- Debug switch ------------------------------------------------------------
# Replaces the scattered message(">>> ... DEBUG ...") and message(dput(...))
# calls. Off by default; set options(CHESKI_DEBUG = TRUE) (or env var
# CHESKI_DEBUG=1) to see them. This is the "remove debug spam" bullet done in a
# way you can still turn back on when diagnosing.
.cheski_debug_on <- function() {
  isTRUE(getOption("CHESKI_DEBUG",
                   default = identical(Sys.getenv("CHESKI_DEBUG"), "1")))
}

dbg <- function(...) {
  if (.cheski_debug_on()) message("[One for All] ", ...)
  invisible(NULL)
}

# --- Single canonicalization entry point -------------------------------------
# The ONLY function any module should call to canonicalize IDs. It dispatches on
# the active species so callers never branch on species themselves.
#   tomato       -> canonicalize_tomato_ids() (Solyc##g######.<v>, case-fixed)
#   arabidopsis  -> normalize_gene_ids()      (toupper + trim)
#   anything else-> trimws()
canon_ids <- function(ids) {
  if (is.null(ids)) return(character(0))
  ids <- as.character(ids)
  if (identical(ACTIVE_SPECIES, "tomato")) {
    canonicalize_tomato_ids(ids)
  } else if (identical(ACTIVE_SPECIES, "arabidopsis")) {
    normalize_gene_ids(ids)
  } else {
    trimws(ids)
  }
}

# --- Single matching strategy ------------------------------------------------
# ids_in(x, table): TRUE/FALSE for each element of x, comparing on canonical
# form. Replaces every `toupper(a) %in% toupper(b)` and
# `normalize_gene_ids(a) %in% b` pattern in the modules.
ids_in <- function(x, table) {
  canon_ids(x) %in% canon_ids(table)
}

# ids_match(x, table): integer index of x within table on canonical form
# (NA where absent). Use when we need positions, e.g. to reorder a matrix.
ids_match <- function(x, table) {
  match(canon_ids(x), canon_ids(table))
}

# --- Small helper: parse a pasted textarea of gene IDs -----------------------
# One consistent splitter for gene_filter / deg_list / universe boxes.
# Splits on whitespace, commas, semicolons, tabs, newlines; drops blanks.
parse_id_textarea <- function(txt, canonicalize = TRUE) {
  if (is.null(txt) || !nzchar(trimws(txt))) return(character(0))
  ids <- unlist(strsplit(txt, "[,;\\s]+", perl = TRUE), use.names = FALSE)
  ids <- trimws(ids)
  ids <- ids[nzchar(ids)]
  ids <- unique(ids)
  if (canonicalize) ids <- canon_ids(ids)
  unique(ids)
}

# --- Clean invisible/odd whitespace from pasted text -------------------------
# Moved here from app.r so every module (not just app.r's own server function)
# can see it: app.r is sourced into its own environment, while this file (via
# global.R) loads into .GlobalEnv, so a function defined only in app.r is
# invisible to code in a separately-sourced module file.
#
# Copy-pasted text (from Excel, PDFs, web pages) often contains invisible
# Unicode space characters (non-breaking space U+00A0, zero-width space
# U+200B, etc.) that LOOK like a normal space but are not matched by R's
# whitespace regex (\s), even with perl = TRUE. Left uncleaned, a pasted
# list of several gene IDs silently collapses into one giant unmatched
# "ID" instead of being split apart. This converts them all to a plain
# ASCII space before any splitting happens.
#
# Quote characters are handled for the same reason. This audience pastes gene
# lists straight out of an R console or dput() output, where every ID arrives
# wrapped in double quotes: "AT1G01010". Those quotes are not whitespace, so
# parse_id_textarea() kept them inside the token, canon_ids() upper-cased the
# whole thing, and the ID silently matched nothing -- reported only as the
# generic "No genes matched. Check your IDs.". Straight and typographic quotes
# are both mapped to a space, which is safe: a quote can never be part of a
# gene ID, and the surrounding trimws()/split step removes the space again.
clean_pasted_ids <- function(x) {
  if (is.null(x)) return(x)
  gsub("[\u00A0\u2000-\u200B\u202F\u205F\u3000\"'\u2018\u2019\u201C\u201D]", " ", x, perl = TRUE)
}
