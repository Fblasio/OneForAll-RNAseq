# ── Pure ID-canonicalization helpers (top-level, dependency-free) ────────────
# Hoisted out of get_species_config()'s tomato branch: neither function's
# body touches ORG_DB or anything species-config-specific, so this is a pure
# code-move. It just means reaching them for testing no longer requires
# loading a Bioconductor annotation package as a side effect.
canonicalize_tomato_ids_impl  <- function(ids) {
  ids <- trimws(as.character(ids))
  parsed <- stringr::str_extract(
    ids,
    stringr::regex("^solyc[0-9]{2}g[0-9]+\\.[0-9]+", ignore_case = TRUE)
  )
  ifelse(is.na(parsed), ids, sub("^solyc", "Solyc", tolower(parsed)))
}

# Logic behind normalize_gene_ids(), parameterized so it's testable without
# a full species config object.
normalize_gene_ids_impl <- function(ids, uppercase) {
  if (uppercase) toupper(trimws(ids)) else trimws(ids)
}

# Columns an OrgDb may use to carry a human-readable gene name, in preference
# order. Defined once and shared by map_gene_symbols() and has_gene_symbols()
# below, so the UI can never offer symbol mapping for an annotation that
# map_gene_symbols() will then refuse.
SYMBOL_COLS <- c("SYMBOL", "GENENAME", "GENE_NAME")

# ── Side ID-map loader (STOPGAP) ────────────────────────────────────────────
# Loads one column of a side CSV that maps this species' gene IDs to an
# external ID space. Used for TWO things the tomato OrgDb cannot provide:
#   ENTREZID -> KEGG (KEGG's 'sly' keys on NCBI GeneID)
#   STRINGID -> STRING (tomato resolves on protein name, not Solyc gene ID)
# Without these, KEGG silently returns an empty result and STRING cannot map.
#
# NOTE ON DRIFT: this deliberately makes the app know mappings the OrgDb does
# not. That is a real inconsistency. The clean fix is to rebuild the OrgDb with
# these columns (AnnotationForge::makeOrgPackage); for ENTREZID the rebuild is
# picked up automatically by the columns(ORG_DB) branch in convert_to_entrez().
# Until then ID_MAP_SOURCE is recorded in the config so the reproducibility
# bundle can state where the mapping came from instead of hiding it.
#
# CSV format: header with a gene-ID column (named after the key type, e.g. GID,
# or simply first) plus any value columns:
#   GID,ENTREZID,STRINGID
#   Solyc01g005000.3,100147723,B1Q3F2
#   Solyc01g005010.4,101263636,
# Blank value cells are dropped for that column only, so one file can back both
# maps at their own (different) coverage levels.
#
# canon_fun is applied to the CSV's ID column so it is canonicalized EXACTLY as
# the app canonicalizes gene IDs — otherwise a prefix/version/case mismatch
# silently matches nothing.
load_id_map_impl <- function(path, key_type, value_col, canon_fun = identity) {
  if (!nzchar(path) || !file.exists(path)) return(NULL)
  df <- tryCatch(
    utils::read.csv(path, colClasses = "character", stringsAsFactors = FALSE),
    error = function(e) NULL
  )
  if (is.null(df) || nrow(df) == 0 || ncol(df) < 2) {
    warning("load_id_map_impl: could not read a 2+ column CSV from ", path)
    return(NULL)
  }
  colnames(df) <- toupper(trimws(colnames(df)))
  value_col <- toupper(value_col)
  if (!value_col %in% colnames(df)) return(NULL)   # column absent -> no map
  id_col <- if (toupper(key_type) %in% colnames(df)) toupper(key_type) else colnames(df)[1]

  out <- data.frame(
    ids = canon_fun(df[[id_col]]),
    val = trimws(df[[value_col]]),
    stringsAsFactors = FALSE
  )
  colnames(out) <- c(key_type, value_col)   # match bitr()'s output contract
  ok <- !is.na(out[[1]]) & nzchar(out[[1]]) &
        !is.na(out[[2]]) & nzchar(out[[2]])
  out <- unique(out[ok, , drop = FALSE])
  rownames(out) <- NULL
  if (nrow(out) == 0) return(NULL)
  message(sprintf("ID side map: %d %s->%s pairs from %s",
                  nrow(out), key_type, value_col, path))
  out
}


get_species_config <- function(species) {
  ITAG_VERSION <- NA_character_
  if (species == "arabidopsis") {
    library(org.At.tair.db)
    ORG_DB <- org.At.tair.db
    ORG_DB_PKG <- "org.At.tair.db"
    SPECIES_LABEL <- "Arabidopsis"
    SPECIES_COMMON <- "Arabidopsis thaliana"
    KEY_TYPE <- "TAIR"
    SYMBOL_COL <- "SYMBOL"
    ENTREZ_COL <- "ENTREZID"
    KEGG_CODE <- "ath"
    GENE_IDS_UPPERCASE <- TRUE
    UNIPROT_TAXID <- "3702"
    GENE_DB_LABEL <- "ThaleMine"
    GENE_DB_URL <- "https://bar.utoronto.ca/thalemine/portal.do?externalids="
    GENE_DETAIL_URL <- "https://www.arabidopsis.org/locus?name="
    GENE_ID_PLACEHOLDER <- "AT1G01010\nAT2G05510\nAT3G18780"
    canonicalize_tomato_ids <- NULL
    collapse_duplicate_rows <- NULL
  } else if (species == "tomato") {
    library(org.Slycopersicum.eg.db) 
    ORG_DB <- org.Slycopersicum.eg.db
    ORG_DB_PKG <- "org.Slycopersicum.eg.db"
    SPECIES_LABEL <- "Tomato"
    SPECIES_COMMON <- "Solanum lycopersicum"
    KEY_TYPE <- "GID"
    SYMBOL_COL <- NA_character_
    ENTREZ_COL <- NA_character_
    KEGG_CODE <- "sly"
    GENE_IDS_UPPERCASE <- FALSE
    UNIPROT_TAXID <- "4081"
    GENE_DB_LABEL <- "SGN"
    GENE_DB_URL <- "https://solgenomics.net/search/quick?term="
    GENE_DETAIL_URL <- "https://solgenomics.net/search/quick?term="
    # Defined once here so no module implements its own copy of the pattern.
    GENE_ID_PLACEHOLDER <- "Solyc01g005000.1\nSolyc02g010000.2\nSolyc03g005000.1"
    ITAG_VERSION <- "ITAG4.0"
    canonicalize_tomato_ids <- canonicalize_tomato_ids_impl
    # Aggregates matrix rows sharing the same canonical ID.
    # rowsum() is the base R function documented for group-wise row
    # aggregation and guarantees unique output row names, satisfying the
    # constraint that matrix/data.frame row names must be unique in R.
    # Reference: R base manual, ?rowsum.
    collapse_duplicate_rows <- function(mat, ids) {
      n_before <- nrow(mat)
      agg <- rowsum(mat, group = ids, reorder = FALSE)
      n_after <- nrow(agg)
      if (n_after < n_before) {
        message(sprintf(
          "collapse_duplicate_rows: %d duplicate tomato IDs aggregated into %d unique rows.",
          n_before - n_after, n_after
        ))
      }
      agg
    }
    
    # message(">>> [SPECIES DEBUG] ACTIVE_SPECIES: ", species)
    # message(">>> [SPECIES DEBUG] ORG_DB_PKG: ", ORG_DB_PKG)
    # message(">>> [SPECIES DEBUG] KEY_TYPE: ", KEY_TYPE)
    # message(">>> [SPECIES DEBUG] available keytypes (primeros): ",
    #         paste(head(AnnotationDbi::keytypes(ORG_DB), 20), collapse = " | "))
  } else {
    stop("Species must be 'arabidopsis' or 'tomato'.")
  }
  
  normalize_gene_ids <- function(ids) normalize_gene_ids_impl(ids, GENE_IDS_UPPERCASE)

  # ── Side ID maps (stopgap; see load_id_map_impl) ──────────────────────────
  # One CSV backs both maps. Path: CHESKI_ID_MAP env var, else
  # data/id_map_<species>.csv. Absent file / absent column -> NULL -> the
  # dependent feature keeps its previous behaviour exactly.
  ID_MAP_PATH <- Sys.getenv(
    "CHESKI_ID_MAP",
    unset = file.path("data", paste0("id_map_", species, ".csv"))
  )
  canon_for_map <- function(x) {
    x <- normalize_gene_ids(x)
    if (!is.null(canonicalize_tomato_ids)) x <- canonicalize_tomato_ids(x)
    x
  }
  ENTREZ_MAP <- load_id_map_impl(ID_MAP_PATH, KEY_TYPE, "ENTREZID", canon_for_map)
  STRING_MAP <- load_id_map_impl(ID_MAP_PATH, KEY_TYPE, "STRINGID", canon_for_map)
  ID_MAP_SOURCE <- if (is.null(ENTREZ_MAP) && is.null(STRING_MAP)) NA_character_ else ID_MAP_PATH

  # TRUE only when ORG_DB actually carries one of SYMBOL_COLS.
  # org.At.tair.db has SYMBOL and GENENAME; org.Slycopersicum.eg.db has
  # neither, so map_gene_symbols() returns NULL for tomato and the heatmap's
  # "Map gene IDs to gene symbols" checkbox was a silent no-op that still
  # named the package in its help text. Exposed so the UI can hide the control
  # instead of offering something the annotation cannot answer.
  #
  # Probed rather than hardcoded per species: rebuild the tomato OrgDb with a
  # SYMBOL column and the control comes back on its own, with no app change.
  has_gene_symbols <- function() {
    available_cols <- tryCatch(columns(ORG_DB), error = function(e) character(0))
    length(intersect(SYMBOL_COLS, available_cols)) > 0
  }

  map_gene_symbols <- function(gene_ids) {
    available_cols <- tryCatch(columns(ORG_DB), error = function(e) character(0))
    sym_col <- intersect(SYMBOL_COLS, available_cols)
    if (length(sym_col) == 0) return(NULL)
    tryCatch(
      mapIds(ORG_DB,
             keys = normalize_gene_ids(gene_ids),
             column = sym_col[1],
             keytype = KEY_TYPE,
             multiVals = "first"),
      error = function(e) NULL
    )
  }
  
  # Returns a data.frame whose FIRST column is the original gene IDs (named
  # after KEY_TYPE) and which has an ENTREZID column — the exact contract
  # bitr() produces, which run_kegg_enrichment() and .remap_rank_for_kegg()
  # both rely on. Returns NULL when no mapping is possible.
  convert_to_entrez <- function(gene_ids) {
    ids <- normalize_gene_ids(gene_ids)

    # 1. Preferred: the OrgDb itself, when it carries ENTREZID. Rebuilding the
    #    OrgDb with an ENTREZID column makes this branch win automatically and
    #    the side map below becomes dead code — no app change needed.
    available_cols <- tryCatch(columns(ORG_DB), error = function(e) character(0))
    if ("ENTREZID" %in% available_cols) {
      return(tryCatch(
        clusterProfiler::bitr(
          ids,
          fromType = KEY_TYPE,
          toType = "ENTREZID",
          OrgDb = ORG_DB
        ),
        error = function(e) NULL
      ))
    }

    # 2. Stopgap: the side CSV. Canonicalize the query IDs the same way the map
    #    was canonicalized at load time, so both sides are directly comparable.
    if (is.null(ENTREZ_MAP)) return(NULL)
    ids_c <- if (!is.null(canonicalize_tomato_ids)) canonicalize_tomato_ids(ids) else ids
    out <- ENTREZ_MAP[ENTREZ_MAP[[1]] %in% unique(ids_c), , drop = FALSE]
    if (nrow(out) == 0) return(NULL)
    rownames(out) <- NULL
    out
  }
  
  list(
    ORG_DB = ORG_DB,
    ORG_DB_PKG = ORG_DB_PKG,
    SPECIES_LABEL = SPECIES_LABEL,
    SPECIES_COMMON = SPECIES_COMMON,
    KEY_TYPE = KEY_TYPE,
    SYMBOL_COL = SYMBOL_COL,
    ENTREZ_COL = ENTREZ_COL,
    ENTREZ_MAP = ENTREZ_MAP,
    STRING_MAP = STRING_MAP,
    ID_MAP_PATH = ID_MAP_PATH,
    # NA when no side map was loaded; otherwise the CSV the mappings came from.
    # Record this in the reproducibility bundle so a run that used a side map is
    # never mistaken for one that used the OrgDb.
    ID_MAP_SOURCE = ID_MAP_SOURCE,
    KEGG_CODE = KEGG_CODE,
    GENE_IDS_UPPERCASE = GENE_IDS_UPPERCASE,
    UNIPROT_TAXID = UNIPROT_TAXID,
    GENE_DB_LABEL = GENE_DB_LABEL,
    GENE_DB_URL = GENE_DB_URL,
    GENE_DETAIL_URL = GENE_DETAIL_URL,
    GENE_ID_PLACEHOLDER = GENE_ID_PLACEHOLDER,
    ITAG_VERSION = ITAG_VERSION,
    normalize_gene_ids = normalize_gene_ids,
    map_gene_symbols = map_gene_symbols,
    has_gene_symbols = has_gene_symbols,
    convert_to_entrez = convert_to_entrez,
    canonicalize_tomato_ids = if (species == "tomato") canonicalize_tomato_ids else NULL,
    collapse_duplicate_rows = if (species == "tomato") collapse_duplicate_rows else NULL
  )
}

species_cfg <- get_species_config(Sys.getenv("CHESKI_SPECIES", unset = "arabidopsis"))

list2env(species_cfg, envir = .GlobalEnv)
ACTIVE_SPECIES <- Sys.getenv("CHESKI_SPECIES", unset = "arabidopsis")