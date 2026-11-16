# ── in-memory cache ────────────────────────────────────────
.enrich_cache <- new.env(hash = TRUE, parent = emptyenv())

.cache_key <- function(...) digest::digest(list(...))

clear_enrichment_cache <- function() {
  rm(list = ls(envir = .enrich_cache), envir = .enrich_cache)
  message("Enrichment cache cleared.")
}



# ── GO ─────────────────────────────────────────────────────
run_go_enrichment <- function(gene_ids,
                              process        = "BP",
                              padj_method    = "BH",
                              pval_cutoff    = 0.05,
                              qval_cutoff    = 0.05,
                              simplify_cutoff = 0.7,
                              use_cache      = TRUE,
                              universe       = NULL,
                              context = "default") {

  key <- .cache_key(sort(gene_ids), process, padj_method,
                    pval_cutoff, qval_cutoff, simplify_cutoff,
                    sort(universe), context, ACTIVE_SPECIES)

  if (use_cache && exists(key, envir = .enrich_cache)) {
    message("GO: returning cached result")
    return(get(key, envir = .enrich_cache))
  }

  geneids_upper <- canon_ids(gene_ids)
  go_warnings <- character(0)
  
  run_params <- list(
    gene_ids        = geneids_upper,
    universe        = if (!is.null(universe)) canon_ids(universe) else NULL,
    process         = process,
    padj_method     = padj_method,
    pval_cutoff     = pval_cutoff,
    qval_cutoff     = qval_cutoff,
    simplify_cutoff = simplify_cutoff,
    # Recorded rather than left implicit. rerun_GO.R emits these values from
    # here instead of hardcoding its own, so the generated script tests the same
    # set of terms as the app and cannot drift from it when either bound
    # changes.
    min_gs_size     = 3,
    max_gs_size     = 2000,
    org_db_pkg      = ORG_DB_PKG,
    key_type        = KEY_TYPE,
    species         = ACTIVE_SPECIES
  )
  
  # message("=== GO DEBUG ===")
  # message("ACTIVE_SPECIES: ", ACTIVE_SPECIES)
  # message("KEY_TYPE: ", KEY_TYPE)
  # message("ORG_DB_PKG: ", ORG_DB_PKG)
  # message("N input genes: ", length(gene_ids))
  # message("N normalized genes: ", length(geneids_upper))
  # message("Head input genes: ", paste(utils::head(gene_ids, 10), collapse = ", "))
  # message("Head normalized genes: ", paste(utils::head(geneids_upper, 10), collapse = ", "))
  # 
  # if (!is.null(universe)) {
  #   universe_norm <- normalize_gene_ids(universe)
  #   message("Universe length raw: ", length(universe))
  #   message("Universe length normalized: ", length(universe_norm))
  #   message("Head universe raw: ", paste(utils::head(universe, 10), collapse = ", "))
  #   message("Head universe normalized: ", paste(utils::head(universe_norm, 10), collapse = ", "))
  # } else {
  #   message("Universe: NULL")
  # }

  # ── Query / background reconciliation ──────────────────────────────────────
  # enrichGO() intersects the query with the universe internally and reports
  # nothing about it, so a background built from a different gene set than the
  # query is indistinguishable on screen from a genuine null result.
  # run_kegg_enrichment() below already reports exactly this loss; without the
  # same report here, the two panels silently disagree about how much of the
  # query was actually tested, and the GO panel's silence reads as "no problem".
  #
  # The query is deliberately NOT pre-filtered: enrichGO keeps doing its own
  # intersection and this block only measures and reports it, so no p-value
  # changes. Note the counts measure membership of the chosen BACKGROUND only;
  # genes further lost because the OrgDb has no annotation for them are a
  # separate matter and are not claimed here.
  #
  # The refusal floor is deliberately lower than KEGG's (which refuses below 5
  # surviving genes): enrichGO runs with minGSSize = 3, so a small query is
  # still a legitimate test. Only a background that excludes the query ENTIRELY
  # is refused.
  if (!is.null(universe) && length(universe) > 0) {
    uni_canon <- unique(canon_ids(universe))
    q_unique  <- unique(geneids_upper)
    n_before  <- length(q_unique)
    n_after   <- length(intersect(q_unique, uni_canon))
    n_lost    <- n_before - n_after
    pct_lost  <- if (n_before > 0) 100 * n_lost / n_before else 0

    run_params$n_query_in_universe <- n_after

    if (n_after == 0) {
      result <- empty_result(sprintf(paste0(
        "GO enrichment not run: none of the %d query gene(s) are present in the ",
        "selected background of %d gene(s). The background must contain every gene ",
        "that could have been selected, so a query lying entirely outside it means ",
        "the two lists were built from different gene sets, or in different ID formats."),
        n_before, length(uni_canon)))
      result$params   <- run_params
      result$warnings <- go_warnings
      if (use_cache) assign(key, result, envir = .enrich_cache)
      return(result)
    }

    if (length(uni_canon) < 100) {
      go_warnings <- c(go_warnings, sprintf(paste0(
        "The background contains only %d gene(s). The hypergeometric test is ",
        "unstable at this size: a single gene can move a term in or out of ",
        "significance. Treat these p-values as indicative only."),
        length(uni_canon)))
    }

    if (n_lost > 0) {
      go_warnings <- c(go_warnings, sprintf(paste0(
        "%d of %d query gene(s) (%.1f%%) are not in the selected background and are ",
        "excluded from the GO test. GeneRatio is computed over the %d remaining ",
        "gene(s), not the original list."),
        n_lost, n_before, pct_lost, n_after))
      if (pct_lost > 20) {
        go_warnings <- c(go_warnings, paste0(
          "More than 20% of the query was lost. This usually means the background and ",
          "the gene list were not built from the same experiment, rather than that the ",
          "gene list is wrong. Treat these p-values as provisional."))
      }
    }

    message(sprintf("GO: background = %d gene(s); query = %d of %d.",
                    length(uni_canon), n_after, n_before))
  }

  go_res <- enrichGO(
    gene     = geneids_upper,
    OrgDb    = ORG_DB,
    keyType  = KEY_TYPE,
    ont           = process,
    pAdjustMethod = padj_method,
    pvalueCutoff  = pval_cutoff,
    qvalueCutoff  = qval_cutoff,
    universe = if (!is.null(universe)) {
      canon_ids(universe)
    } else NULL,
    readable      = FALSE,
    minGSSize     = 3,
    maxGSSize     = 2000,
    pool          = ifelse(process == "ALL", TRUE, FALSE)
  )

  # ── Apply the p/q cutoffs explicitly, once, before anything reads the table ──
  # enrichGO() returns an object whose @result slot holds EVERY TESTED term.
  # In clusterProfiler >= 4 / DOSE >= 4 the pvalueCutoff/qvalueCutoff filter
  # lives in DOSE:::get_enriched(), which runs ONLY when the object is passed
  # through as.data.frame() / dim() / `[`. Reading @result directly skips it.
  #
  # Relying on clusterProfiler::simplify() to apply them is not sufficient: it
  # filters only as a side effect of its first statement,
  # `res <- as.data.frame(x)`, which makes significance filtering conditional on
  # a redundancy step that does not always run:
  #   * ontology = "ALL" skips simplify() (see the branch below) -> unfiltered;
  #   * the tryCatch below falls back to the raw `go_res` -> unfiltered.
  # Both returned every tested term, p.adjust up to 1, presented as results.
  #
  # Filtering here makes every downstream path inherit an already-filtered
  # table, and leaves simplify() doing only the job its name describes.
  # as.data.frame() is idempotent, so simplify() re-filtering afterwards is a
  # no-op -- no p-value and no term ordering changes on the happy path.
  n_terms_tested <- if (is.null(go_res)) 0L else nrow(go_res@result)
  if (!is.null(go_res)) go_res@result <- as.data.frame(go_res)

  if (is.null(go_res) || nrow(go_res@result) == 0) {
    result <- empty_result("No GO terms found below the chosen p/q-value cutoffs.")
    # Recorded on this path too: "4000 terms tested, none significant" and
    # "nothing was testable at all" are different outcomes and the bundle
    # should be able to tell them apart.
    run_params$n_terms_tested <- n_terms_tested
    result$params <- run_params
    # The warnings must travel on this path too. An empty table is the normal
    # way a background or ID mismatch surfaces, which is precisely when the
    # "N of M query genes are not in the background" note carries the
    # explanation the user needs.
    result$warnings <- go_warnings
    if (use_cache) assign(key, result, envir = .enrich_cache)
    return(result)
  }

  # Significant terms, i.e. what simplify() is about to receive. Named as
  # before so the fields below keep their meaning; n_terms_tested is the new
  # number and is what the p.adjust caveat must quote.
  n_terms_raw <- nrow(go_res@result)

  ego_final <- if (process != "ALL") {
    tryCatch(
      clusterProfiler::simplify(go_res, cutoff = simplify_cutoff,
                                by = "p.adjust", select_fun = min),
      error = function(e) {
        go_warnings <<- c(go_warnings, paste0(
          "simplify() failed (", conditionMessage(e), "): the table is not ",
          "redundancy-filtered and may contain overlapping parent/child terms."))
        go_res
      }
    )
  } else {
    # The old wording justified this with "semantic similarity is not
    # comparable across BP/MF/CC". That is not why it is skipped:
    # clusterProfiler's simplify() accepts ont = "ALL" (enrichGO relabels the
    # object to "GOALL") and dispatches simplify_ALL(), which splits by the
    # ONTOLOGY column and runs the similarity per ontology. It is skipped here
    # for cost -- one godata() build per ontology. The claim was removed rather
    # than restated, so the note says only what is true.
    go_warnings <- c(go_warnings, paste0(
      "ontology = ALL: simplify() is not applied, so redundant parent/child ",
      "terms remain. p.adjust is pooled across all three ontologies."))
    go_res
  }

  n_terms_simplified <- nrow(ego_final@result)

  if (n_terms_simplified < n_terms_raw) {
    # Both counts in this sentence are drawn from the SIGNIFICANT population, so
    # simplify() is credited only with the terms it actually merged. Taking
    # n_terms_raw from the unfiltered table instead would attribute every
    # non-significant term to redundancy as well, turning "18 of 34 redundant"
    # into "3966 of 4000". The p.adjust clause quotes n_terms_tested separately,
    # because that is the number BH adjusted over.
    #
    # write_result_csv() prepends this text to every exported GO CSV, so the
    # numbers travel into supplementary files and have to be right there too.
    go_warnings <- c(go_warnings, sprintf(paste0(
      "simplify(cutoff = %.2f) removed %d of %d significant terms as redundant ",
      "by semantic similarity. p.adjust was computed over all %d tested terms, ",
      "not the %d shown. Representatives were chosen by lowest p.adjust, which ",
      "favours the more general term of each redundant pair."),
      simplify_cutoff, n_terms_raw - n_terms_simplified, n_terms_raw,
      n_terms_tested, n_terms_simplified))
  }

  run_params$n_terms_tested          <- n_terms_tested
  run_params$n_terms_before_simplify <- n_terms_raw
  run_params$n_terms_after_simplify  <- n_terms_simplified

  ego_final@result$geneID <- vapply(
    ego_final@result$geneID,
    function(x) {
      ids   <- canon_ids(trimws(strsplit(x, "/")[[1]]))
      valid <- ids %in% geneids_upper
      paste(ids[valid], collapse = "/") 
    },
    character(1)
  )
  ego_final@result$Count <- vapply(
    ego_final@result$geneID,
    function(x) if (nzchar(x)) length(strsplit(x, "/")[[1]]) else 0L,
    integer(1)
  )
  ego_final@result <- ego_final@result[ego_final@result$Count > 0, ]

  # The `nrow == 0` guard above runs at the cutoff stage only. Two later steps
  # can ALSO empty the table -- simplify(), and the Count > 0 filter on the line
  # above -- and until now both fell through to the success return, which
  # carries no $reason. The GO Table tab then renders an empty DT with column
  # headers and the GO Plot tab takes its plotting branch (the `Description`
  # column IS present), so the user got two blank panels and no sentence.
  #
  # Reproducible with a degenerate background: 400 query genes against a
  # 40-gene custom universe makes qvalue() fail inside enricher_internal, DOSE
  # returns rows whose every column is NA, the cutoff guard sees nrow > 0, and
  # the geneID re-map then drops all of them.
  if (nrow(ego_final@result) == 0) {
    result <- empty_result(paste0(
      "No GO terms survived after redundancy filtering and gene re-mapping. ",
      "This usually means the background was too small or too different from ",
      "the query for enrichGO to compute usable statistics (q-values could not ",
      "be estimated), rather than that no biology is present. Try the default ",
      "whole-genome background, or a larger expression-derived one."))
    run_params$n_terms_after_gene_remap <- 0L
    result$params   <- run_params
    result$warnings <- go_warnings
    if (use_cache) assign(key, result, envir = .enrich_cache)
    return(result)
  }

  ego_final <- tryCatch({
    gr <- as.character(ego_final@result$GeneRatio)
    br <- as.character(ego_final@result$BgRatio)
    fe <- DOSE::parse_ratio(gr) / DOSE::parse_ratio(br)
    fe[!is.finite(fe)] <- NA
    ego_final@result$FoldEnrichment <- fe
    ego_final
  }, error = function(e) {
    message("FoldEnrichment computation failed: ", e$message)
    ego_final
  })

  result <- list(table = as.data.frame(ego_final@result), object = ego_final,
                 params = run_params, warnings = go_warnings)
  if (use_cache) assign(key, result, envir = .enrich_cache)
  result
}

# ── KEGG ───────────────────────────────────────────────────
run_kegg_enrichment <- function(gene_ids,
                                padj_method = "BH",
                                pval_cutoff = 0.05,
                                qval_cutoff = 0.05,
                                use_cache   = TRUE,
                                universe    = NULL) {
  # `universe` joins the cache key: the same gene list against a different
  # background is a different result, and omitting it here would serve a
  # stale genome-background table after the background changes.
  key <- .cache_key(sort(gene_ids), "kegg", padj_method, pval_cutoff, qval_cutoff,
                    sort(universe), ACTIVE_SPECIES)

  if (use_cache && exists(key, envir = .enrich_cache)) {
    message("KEGG: returning cached result")
    return(get(key, envir = .enrich_cache))
  }
  kegg_warnings <- character(0)

  # ── Run parameters ──────────────────────────────────────────────────────────
  # KEGG was the only enrichment path with no params block, so keggServer had
  # nothing to accumulate and the reproducibility bundle described KEGG runs by
  # re-reading live inputs at download time -- two of which pointed at widget
  # IDs that do not exist. Built here, before anything can return, and attached
  # to EVERY exit path so a refusal is as traceable as a success.
  # Mirrors run_go_enrichment()'s run_params field-for-field where the two
  # analyses share a concept, so both read the same way in manifest.json.
  run_params <- list(
    gene_ids      = canon_ids(gene_ids),
    universe      = if (!is.null(universe)) canon_ids(universe) else NULL,
    padj_method   = padj_method,
    pval_cutoff   = pval_cutoff,
    qval_cutoff   = qval_cutoff,
    min_gs_size   = 4,
    max_gs_size   = 1000,
    org_db_pkg    = ORG_DB_PKG,
    key_type      = KEY_TYPE,
    kegg_code     = KEGG_CODE,
    species       = ACTIVE_SPECIES,
    id_map_source = if (is.na(ID_MAP_SOURCE)) NA_character_ else ID_MAP_SOURCE,
    n_genes_input = length(unique(canon_ids(gene_ids)))
  )

  # Arabidopsis KEGG uses TAIR IDs directly; other species need Entrez/NCBI IDs.
  if (ACTIVE_SPECIES == "arabidopsis") {
    kegg_input <- normalize_gene_ids(gene_ids)
  } else {
    ids_canon <- canon_ids(gene_ids)
    entrez_df <- convert_to_entrez(ids_canon)
    if (is.null(entrez_df) || nrow(entrez_df) == 0) {
      result <- empty_result(paste0(
        "KEGG needs Entrez/NCBI IDs, but none of the ", length(unique(ids_canon)),
        " input gene ID(s) could be mapped for this species. KEGG enrichment skipped."))
      run_params$n_genes_mapped <- 0L
      result$params <- run_params
      if (use_cache) assign(key, result, envir = .enrich_cache)
      return(result)
    }
    kegg_input <- entrez_df$ENTREZID

    # ── Partial ID-mapping loss ───────────────────────────────────────────────
    # Above this point only TOTAL mapping failure is reported. Partial loss was
    # reported nowhere unless a universe was supplied, because the only counting
    # block in this function lives inside `if (!is.null(universe))` below — and
    # the default background for both KEGG panels is NULL (whole genome). So a
    # query where most genes have no Entrez ID ran silently on the remainder,
    # and GeneRatio's denominator quietly became the mapped subset while the
    # panel still reported the user's original selection.
    #
    # Counted on entrez_df[[1]], NOT on ENTREZID: column 1 is the KEY_TYPE
    # column in both branches of convert_to_entrez() (bitr's contract, and
    # load_id_map_impl() renames the side map's columns to match it), so this
    # counts GENES. Counting ENTREZID would count records instead, which is a
    # different and larger number wherever one gene carries several Entrez IDs.
    #
    # This is a warning, not a refusal: for species whose annotation genuinely
    # lacks Entrez IDs the shortfall is a permanent ceiling, not a fixable gap,
    # so refusing would block a legitimate analysis. Reporting is the only
    # remedy available. 20% mirrors the escalation point already used for
    # background loss below, so the two messages stay consistent.
    n_query      <- length(unique(ids_canon))
    n_mapped     <- length(unique(entrez_df[[1]]))
    n_unmapped   <- n_query - n_mapped
    pct_unmapped <- if (n_query > 0) 100 * n_unmapped / n_query else 0
    run_params$n_genes_mapped <- n_mapped

    if (n_unmapped > 0) {
      map_src_query <- if (is.na(ID_MAP_SOURCE)) "the OrgDb" else ID_MAP_SOURCE
      kegg_warnings <- c(kegg_warnings, sprintf(paste0(
        "%d of %d query gene(s) (%.1f%%) have no Entrez/NCBI ID for %s and were not ",
        "tested by KEGG. GeneRatio and Count are computed over the %d gene(s) that ",
        "mapped, not over the original list. This is a coverage limit of the ID ",
        "source (%s), not a filter you chose."),
        n_unmapped, n_query, pct_unmapped, SPECIES_COMMON, n_mapped, map_src_query))
      if (pct_unmapped > 20) {
        kegg_warnings <- c(kegg_warnings, paste0(
          "More than 20% of the query has no KEGG-compatible ID. Genes lacking an ",
          "NCBI GeneID are typically the less-characterised ones, so the tested set ",
          "is biased toward well-annotated genes. Quote the tested gene count, not ",
          "the selected one, when reporting these results."))
      }
      message(sprintf("KEGG: ID mapping = %d of %d query gene(s).", n_mapped, n_query))
    }
  }

  # Background must live in the SAME ID space as `kegg_input`, so it goes
  # through the identical conversion. universe = NULL means the CALLER asked for
  # no restriction -> enrichKEGG uses every gene in the TERM2GENE table. That is
  # a legitimate choice and is left untouched. What is not legitimate is a
  # requested background silently degrading into that.
  kegg_universe <- NULL
  if (!is.null(universe) && length(universe) > 0) {
    if (ACTIVE_SPECIES == "arabidopsis") {
      kegg_universe <- unique(normalize_gene_ids(universe))
    } else {
      uni_df <- convert_to_entrez(canon_ids(universe))
      if (!is.null(uni_df) && nrow(uni_df) > 0) {
        kegg_universe <- unique(uni_df$ENTREZID)
      }
    }

    map_src <- if (is.na(ID_MAP_SOURCE)) "the OrgDb" else ID_MAP_SOURCE

    # (a) Background maps to nothing: refuse rather than fall back to the whole
    # genome, which is the most permissive background available and the opposite
    # of what was requested.
    if (is.null(kegg_universe) || length(kegg_universe) == 0) {
      result <- empty_result(sprintf(paste0(
        "KEGG enrichment not run: none of the %d background gene(s) could be mapped ",
        "to the ID space KEGG uses for %s. Running anyway would silently replace the ",
        "chosen background with the whole genome, which inflates every p-value in the ",
        "table. Check the ID format of the background list, or the ID map (%s). To use ",
        "the whole genome, select that background explicitly."),
        length(unique(universe)), SPECIES_COMMON, map_src))
      # Carry the query-side mapping warning through: this is the only early
      # return in the function that dropped `kegg_warnings`, and it is reachable
      # with a query-mapping warning already pending.
      result$params   <- run_params
      result$warnings <- kegg_warnings
      if (use_cache) assign(key, result, envir = .enrich_cache)
      return(result)
    }

    # (b) Background maps but is too small for a stable hypergeometric test.
    if (length(kegg_universe) < 100) {
      kegg_warnings <- c(kegg_warnings, sprintf(paste0(
        "The background contains only %d gene(s) after KEGG ID mapping. The ",
        "hypergeometric test is unstable at this size: a single gene can move a ",
        "pathway in or out of significance. Treat these p-values as indicative only."),
        length(kegg_universe)))
    }

    # (c) Query genes absent from the background. The background defines the
    # population the test samples from, so a query gene not in it cannot be
    # drawn. Enlarging the background to accommodate it changes the test's
    # denominator and the p-values would describe a background that was never
    # chosen. Drop those genes from the query and report the loss.
    n_before   <- length(unique(kegg_input))
    kegg_input <- intersect(unique(kegg_input), kegg_universe)
    n_after    <- length(kegg_input)
    n_lost     <- n_before - n_after
    pct_lost   <- if (n_before > 0) 100 * n_lost / n_before else 0

    if (n_after < 5) {
      result <- empty_result(sprintf(paste0(
        "KEGG enrichment not run: only %d of the %d query gene(s) are present in the ",
        "selected background. The background must contain every gene that could have ",
        "been selected, so a query gene missing from it means the two lists were built ",
        "from different gene sets, or mapped to KEGG IDs at different coverage."),
        n_after, n_before))
      result$params   <- run_params
      result$warnings <- kegg_warnings
      if (use_cache) assign(key, result, envir = .enrich_cache)
      return(result)
    }

    if (n_lost > 0) {
      kegg_warnings <- c(kegg_warnings, sprintf(paste0(
        "%d of %d query genes (%.1f%%) are not in the selected background and were ",
        "excluded from the KEGG test. GeneRatio is computed over the %d remaining ",
        "genes, not the original list. Adding them to the background instead would ",
        "change the test's denominator and invalidate the p-values."),
        n_lost, n_before, pct_lost, n_after))
      if (pct_lost > 20) {
        kegg_warnings <- c(kegg_warnings, paste0(
          "More than 20% of the query was lost. This usually means the background and ",
          "the gene list were not built from the same experiment, rather than that the ",
          "gene list is wrong. Treat these p-values as provisional."))
      }
    }

    message(sprintf("KEGG: background = %d gene(s); query = %d of %d.",
                    length(kegg_universe), n_after, n_before))
  }

  kegg_res <- enrichKEGG(
    gene     = kegg_input,
    organism = KEGG_CODE,
    universe = kegg_universe,
    minGSSize         = 4,
    maxGSSize         = 1000,
    pAdjustMethod     = padj_method,
    pvalueCutoff      = pval_cutoff,
    qvalueCutoff      = qval_cutoff,
    use_internal_data = FALSE
  )

  # ── Apply the p/q cutoffs explicitly (same fix as run_go_enrichment) ────────
  # enrichKEGG() also returns an unfiltered @result, and KEGG never had a
  # simplify() step to filter it by accident -- so the table, the dotplot and
  # the exported CSV listed EVERY tested pathway, p.adjust up to 1, and the
  # cutoff boxes had no effect on what was displayed. make_enrich_dotplot()
  # ranks by GeneRatio rather than significance, so the plot showed the top 20
  # gene ratios regardless of p-value.
  if (!is.null(kegg_res)) {
    run_params$n_pathways_tested <- nrow(kegg_res@result)
    kegg_res@result <- as.data.frame(kegg_res)
  }

  if (is.null(kegg_res) || nrow(kegg_res@result) == 0) {
    result <- empty_result(sprintf(
      "No KEGG pathways found below the chosen cutoffs (%d gene ID(s) used as KEGG input).",
      length(unique(kegg_input))))
    result$params   <- run_params
    result$warnings <- kegg_warnings
    if (use_cache) assign(key, result, envir = .enrich_cache)
    return(result)
  }

  result <- list(table    = as.data.frame(kegg_res@result),
                 object   = kegg_res,
                 params   = run_params,
                 warnings = kegg_warnings)
  if (use_cache) assign(key, result, envir = .enrich_cache)
  result
}

# ── Custom dotplot ─────────────────────────────────────────
# pinned: character vector of Description values that must always appear in the plot.
# Remaining slots (top_n - n_pinned) are filled with top-ranked free terms by GeneRatio.
make_enrich_dotplot <- function(result_df, title = "", top_n = 20, pinned = character(0), font_size = 13) {
  req_cols <- c("Description", "GeneRatio", "p.adjust", "Count")
  missing  <- setdiff(req_cols, colnames(result_df))
  if (length(missing) > 0)
    stop("Missing columns: ", paste(missing, collapse = ", "))

  df <- result_df
  df$Description  <- as.character(df$Description)
  df$GeneRatio    <- as.character(df$GeneRatio)
  df$p.adjust     <- as.numeric(df$p.adjust)
  df$Count        <- as.numeric(df$Count)

  df$GeneRatioNum <- sapply(df$GeneRatio, function(x) {
    parts <- strsplit(x, "/")[[1]]
    if (length(parts) == 2) as.numeric(parts[1]) / as.numeric(parts[2])
    else NA_real_
  })

  df <- df[!is.na(df$GeneRatioNum) & !is.na(df$p.adjust), ]
  df <- df[order(df$GeneRatioNum, decreasing = TRUE), ]
  # Apply pinned logic: pinned terms always shown; remaining slots filled with top free terms
  if (length(pinned) > 0) {
    is_pinned <- df$Description %in% pinned
    df_pinned <- df[is_pinned,  ]
    df_free   <- df[!is_pinned, ]
    n_free    <- max(0L, top_n - nrow(df_pinned))
    df        <- rbind(df_free[seq_len(min(n_free, nrow(df_free))), ], df_pinned)
    df        <- df[order(df$GeneRatioNum, decreasing = TRUE), ]
  } else {
    df <- head(df, top_n)
  }
  df$Description <- factor(df$Description, levels = rev(df$Description))

  ggplot(df, aes(x = GeneRatioNum, y = Description,
                 size = Count, color = p.adjust)) +
    geom_point() +
    scale_color_gradient(low = "#C0392B", high = "#00008b", name = "p.adjust") +
    scale_size_continuous(name = "Gene count", range = c(3, 10)) +
    xlab("Gene Ratio") + ylab(NULL) +
    ggtitle(title) +
    theme_bw(base_size = font_size) +
    theme(axis.text.y  = element_text(size = font_size),
          axis.text.x  = element_text(size = font_size),
          axis.title.x = element_text(size = font_size + 2, face = "bold"),
          plot.title   = element_text(size = font_size + 4, face = "bold"),
          legend.text  = element_text(size = font_size - 1),
          legend.title = element_text(size = font_size,     face = "bold"))
}

# ── Flexible, publication-oriented enrichment plot ──────────────────────────
# make_enrichment_plot(): superset of make_enrich_dotplot()/.make_dotplot()
# (both kept as-is elsewhere for backward compatibility), with plot-type,
# x-axis metric, and color choices exposed as arguments instead of being
# hard-coded.
#
# Written against the {Description, GeneRatio, p.adjust, Count} columns that
# both run_go_enrichment() and run_kegg_enrichment() already produce, plus
# the optional FoldEnrichment column run_go_enrichment() adds. This makes it
# usable from GO today and from KEGG later without any changes here — see
# the FoldEnrichment fallback below.
#
# style        : "dot" (classic dotplot), "lollipop" (point + stem from 0),
#                or "bar" (horizontal bar chart).
# x_metric     : "GeneRatio" (k/n, always available) or "FoldEnrichment"
#                ("Enrichment Factor" in the UI; observed/expected ratio).
#                Falls back to GeneRatio with a message() if the column is
#                missing or entirely NA (e.g. a result table that predates
#                this feature, or a source that never computed it).
# color_mode   : "gradient" (color/fill mapped to a significance metric) or
#                "flat" (single user-chosen color, no color legend).
# color_metric : "p.adjust" or "neglog10_padj" (-log10(p.adjust); larger =
#                more significant, sometimes preferred so "more intense
#                color" reads as "more significant").
# size_by_count: TRUE maps point/lollipop-head size to Count; FALSE uses a
#                fixed size (ignored for style = "bar", whose length already
#                encodes the x metric).
# label_wrap   : character width to wrap long term Descriptions onto
#                multiple lines (0/NULL disables wrapping).
# theme_style  : "bw", "minimal", or "classic" — swaps the ggplot2 base theme.
make_enrichment_plot <- function(result_df,
                                  title         = "",
                                  top_n         = 20,
                                  pinned        = character(0),
                                  font_size     = 13,
                                  style         = c("dot", "lollipop", "bar"),
                                  x_metric      = c("GeneRatio", "FoldEnrichment"),
                                  color_mode    = c("gradient", "flat"),
                                  color_metric  = c("p.adjust", "neglog10_padj"),
                                  color_low     = "#3B7DD8",
                                  color_high    = "#C0392B",
                                  flat_color    = "#3B7DD8",
                                  size_by_count = TRUE,
                                  point_size    = 4,
                                  label_wrap    = 45,
                                  theme_style   = c("bw", "minimal", "classic"),
                                  x_lab         = NULL,
                                  sort_desc     = TRUE) {

  style        <- match.arg(style)
  x_metric     <- match.arg(x_metric)
  color_mode   <- match.arg(color_mode)
  color_metric <- match.arg(color_metric)
  theme_style  <- match.arg(theme_style)

  req_cols <- c("Description", "GeneRatio", "p.adjust", "Count")
  missing  <- setdiff(req_cols, colnames(result_df))
  if (length(missing) > 0)
    stop("Missing columns: ", paste(missing, collapse = ", "))

  df <- result_df
  df$Description <- as.character(df$Description)
  df$GeneRatio    <- as.character(df$GeneRatio)
  df$p.adjust     <- as.numeric(df$p.adjust)
  df$Count        <- as.numeric(df$Count)

  df$GeneRatioNum <- sapply(df$GeneRatio, function(x) {
    parts <- strsplit(x, "/")[[1]]
    if (length(parts) == 2) as.numeric(parts[1]) / as.numeric(parts[2]) else NA_real_
  })

  # FoldEnrichment ("Enrichment Factor"): use it only if present and usable;
  # otherwise transparently fall back to GeneRatio so callers never error
  # out just because a source table doesn't compute it (e.g. future KEGG use).
  has_fold <- "FoldEnrichment" %in% colnames(df) &&
    !all(is.na(suppressWarnings(as.numeric(df$FoldEnrichment))))
  if (identical(x_metric, "FoldEnrichment") && !has_fold) {
    message("make_enrichment_plot: FoldEnrichment not available in this ",
            "result table; falling back to GeneRatio.")
    x_metric <- "GeneRatio"
  }
  df$XMetric <- if (identical(x_metric, "FoldEnrichment")) {
    as.numeric(df$FoldEnrichment)
  } else {
    df$GeneRatioNum
  }
  x_axis_label <- x_lab %||% if (identical(x_metric, "FoldEnrichment")) {
    "Enrichment Factor (Fold Enrichment)"
  } else {
    "Gene Ratio"
  }

  df$ColorMetric <- if (identical(color_metric, "neglog10_padj")) {
    -log10(df$p.adjust)
  } else {
    df$p.adjust
  }
  color_legend_name <- if (identical(color_metric, "neglog10_padj")) {
    "-log10(p.adjust)"
  } else {
    "p.adjust"
  }

  df <- df[!is.na(df$XMetric) & !is.na(df$p.adjust), ]
  if (nrow(df) == 0) {
    return(
      ggplot() +
        annotate("text", x = 0.5, y = 0.5,
                 label = "No terms available to plot with the current settings.",
                 size = 6, color = "gray40") +
        theme_void()
    )
  }
  df <- df[order(df$XMetric, decreasing = sort_desc), ]

  # Pinned logic: pinned terms always shown; remaining slots filled with the
  # top-ranked free terms by the CHOSEN x metric (not hard-coded to GeneRatio).
  if (length(pinned) > 0) {
    is_pinned <- df$Description %in% pinned
    df_pinned <- df[is_pinned,  ]
    df_free   <- df[!is_pinned, ]
    n_free    <- max(0L, top_n - nrow(df_pinned))
    df        <- rbind(df_free[seq_len(min(n_free, nrow(df_free))), ], df_pinned)
    df        <- df[order(df$XMetric, decreasing = sort_desc), ]
  } else {
    df <- head(df, top_n)
  }

  if (!is.null(label_wrap) && label_wrap > 0) {
    df$Description <- vapply(df$Description,
                              function(s) paste(strwrap(s, width = label_wrap), collapse = "\n"),
                              character(1))
  }
  # Order retained top-to-bottom as sorted above (rev() so ggplot's bottom-up
  # discrete axis reads top-to-bottom the way the data frame is ordered).
  df$Description <- factor(df$Description, levels = rev(unique(df$Description)))

  p <- ggplot(df, aes(x = XMetric, y = Description))

  if (identical(style, "lollipop")) {
    p <- p + geom_segment(aes(x = 0, xend = XMetric, y = Description, yend = Description),
                           color = "grey60", linewidth = 0.6)
  }

  if (identical(style, "bar")) {
    if (identical(color_mode, "gradient")) {
      p <- p + geom_col(aes(fill = ColorMetric), width = 0.7) +
        scale_fill_gradient(low = color_low, high = color_high, name = color_legend_name)
    } else {
      p <- p + geom_col(fill = flat_color, width = 0.7)
    }
  } else {
    # Four explicit, unambiguous branches instead of trying to merge partial
    # aes()/NULL-size objects — safer than it looks, and easy to verify.
    if (identical(color_mode, "gradient") && isTRUE(size_by_count)) {
      p <- p + geom_point(aes(size = Count, color = ColorMetric)) +
        scale_color_gradient(low = color_low, high = color_high, name = color_legend_name) +
        scale_size_continuous(name = "Gene count", range = c(3, 10))
    } else if (identical(color_mode, "gradient") && !isTRUE(size_by_count)) {
      p <- p + geom_point(aes(color = ColorMetric), size = point_size) +
        scale_color_gradient(low = color_low, high = color_high, name = color_legend_name)
    } else if (!identical(color_mode, "gradient") && isTRUE(size_by_count)) {
      p <- p + geom_point(aes(size = Count), color = flat_color) +
        scale_size_continuous(name = "Gene count", range = c(3, 10))
    } else {
      p <- p + geom_point(color = flat_color, size = point_size)
    }
  }

  p <- p +
    xlab(x_axis_label) + ylab(NULL) +
    ggtitle(title)

  p <- p + switch(theme_style,
                   bw      = theme_bw(base_size = font_size),
                   minimal = theme_minimal(base_size = font_size),
                   classic = theme_classic(base_size = font_size))

  p + theme(axis.text.y  = element_text(size = font_size),
            axis.text.x  = element_text(size = font_size),
            axis.title.x = element_text(size = font_size + 2, face = "bold"),
            plot.title   = element_text(size = font_size + 4, face = "bold"),
            legend.text  = element_text(size = font_size - 1),
            legend.title = element_text(size = font_size,     face = "bold"))
}

# ── Shared "Plot style" control panel for the ORA plots ──────────────────────
# One definition of the controls that drive make_enrichment_plot(), used by
# both goPlotUI() and keggPlotUI(). Writing the block out inline in each module
# would produce two ~60-line duplicates that drift apart on the next change, so
# the two panels are guaranteed identical by construction rather than by
# discipline.
#
# The input IDs are exactly the argument names make_enrichment_plot() takes, so
# a module's plot reactive is a straight pass-through with %||% defaults, and
# the defaults here are the same ones that function already uses -- adding this
# panel to a module cannot change the figure it drew before.
#
# ns          : the calling module's NS() function. Needed both to namespace the
#               inputs and to build the conditionalPanel() conditions, which are
#               evaluated client-side against the fully-qualified input name.
# terms_label : "GO terms" / "pathways" -- the noun used in the help text, so
#               the X-axis tooltip reads correctly in either panel.
# open        : FALSE renders the panel folded (native <details>, closed by
#               default). These are set-once controls and both tabs are already
#               tall, so folded is the default.
#
# Native <details>/<summary> via htmltools::tag() rather than tags$details():
# same reasoning as gsea_module.R's own style panel -- no JS accordion
# dependency, and tag() is the generic constructor that exists in every
# htmltools version, unlike the tags$xxx() shorthands for newer HTML5 elements.
enrichPlotStyleUI <- function(ns, terms_label = "terms", open = FALSE) {
  x_axis_help <- paste0(
    "Gene Ratio = the share of YOUR list in the ", sub("s$", "", terms_label),
    ", so big ", terms_label, " score high. Enrichment Factor = how many times ",
    "more often than expected by chance (1 = none), which favours small ",
    "specific ", terms_label, ". Neither one shows significance.")

  details_args <- list(
    style = paste0("background:#f8f9fa; border:1px solid #dee2e6; border-radius:4px; ",
                   "padding:6px 10px 10px 10px; margin:8px 0;"),
    htmltools::tag("summary", list(
      style = "cursor:pointer; font-weight:bold; color:#1F3864; padding:4px 0;",
      "\U0001f3a8 Plot style"
    )),
    div(
      style = "margin-top:8px;",
      fluidRow(
        column(6, selectInput(ns("plot_style"), "Plot type",
                              choices  = c("Dot plot"  = "dot",
                                           "Lollipop"  = "lollipop",
                                           "Bar chart" = "bar"),
                              selected = "dot")),
        column(6, selectInput(ns("x_metric"), lbl("X-axis", x_axis_help),
                              choices  = c("Gene Ratio"                    = "GeneRatio",
                                           "Enrichment Factor (Fold Enr.)" = "FoldEnrichment"),
                              selected = "GeneRatio"))
      ),
      fluidRow(
        column(6, selectInput(ns("color_mode"), "Color",
                              choices  = c("By significance (gradient)" = "gradient",
                                           "Single color"               = "flat"),
                              selected = "gradient")),
        column(6, conditionalPanel(
          condition = paste0("input['", ns("color_mode"), "'] == 'gradient'"),
          selectInput(ns("color_metric"), "Color by",
                      choices  = c("p.adjust"         = "p.adjust",
                                   "-log10(p.adjust)" = "neglog10_padj"),
                      selected = "p.adjust")
        ))
      ),
      fluidRow(
        column(6, conditionalPanel(
          condition = paste0("input['", ns("color_mode"), "'] == 'gradient'"),
          fluidRow(
            column(6, colourpicker::colourInput(ns("color_low"),  "Low",  value = "#3B7DD8")),
            column(6, colourpicker::colourInput(ns("color_high"), "High", value = "#C0392B"))
          )
        )),
        column(6, conditionalPanel(
          condition = paste0("input['", ns("color_mode"), "'] == 'flat'"),
          colourpicker::colourInput(ns("flat_color"), "Point / bar color", value = "#3B7DD8")
        ))
      ),
      fluidRow(
        column(6, conditionalPanel(
          condition = paste0("input['", ns("plot_style"), "'] != 'bar'"),
          checkboxInput(ns("size_by_count"), "Size points by gene count", value = TRUE)
        )),
        column(6, selectInput(ns("theme_style"), "Theme",
                              choices  = c("Black & white" = "bw",
                                           "Minimal"       = "minimal",
                                           "Classic"       = "classic"),
                              selected = "bw"))
      ),
      fluidRow(
        column(12, sliderInput(ns("label_wrap"),
                               paste0("Wrap ", sub("s$", "", terms_label),
                                      " labels at (characters, 0 = off)"),
                               min = 0, max = 80, value = 45, step = 5, ticks = FALSE))
      )
    )
  )
  if (isTRUE(open)) details_args$open <- NA   # HTML boolean attribute
  htmltools::tag("details", details_args)
}

# ── Export sizing for the ORA plots ──────────────────────────────────────────
# The downloaded figure used to be a fixed 10x8in PDF and a
# png(1200, 900, res = 150) -- which is an 8x6in canvas, not a big one. ggplot
# font sizes are in POINTS, an absolute physical unit, so the same 13pt label
# that occupies 1.6% of the on-screen canvas occupies 3.0% of that PNG: the text
# and the points look ~1.85x larger in the download, and with many terms the
# y-axis labels collide. Neither download read the on-screen height slider, and
# the two disagreed with each other.
#
# These three helpers fix the sizing in one place so the GO and KEGG panels
# cannot drift apart, the same reason enrichPlotStyleUI() is shared.

# Height a dot/bar plot needs so its y-axis labels do not collide.
# Each term needs font_size * line_spacing * n_lines POINTS of vertical space;
# overhead_in is the title + axis + legend + margins, measured at a consistent
# 0.7in across the screen, PDF and PNG devices.
#
# n_lines: how many lines the LONGEST wrapped label takes. Callers pass the
# maximum over the whole result table rather than over the terms actually
# drawn -- deliberately conservative, because overestimating costs a little
# white space while underestimating puts the labels back on top of each other.
enrich_plot_height_in <- function(n_terms, font_size = 13, n_lines = 1,
                                  line_spacing = 1.2, overhead_in = 0.7,
                                  min_in = 3) {
  n_terms   <- suppressWarnings(as.numeric(n_terms))
  font_size <- suppressWarnings(as.numeric(font_size))
  n_lines   <- suppressWarnings(as.numeric(n_lines))
  if (!isTRUE(is.finite(n_terms)) || n_terms < 1) return(min_in)
  if (!isTRUE(is.finite(font_size)) || font_size <= 0) font_size <- 13
  if (!isTRUE(is.finite(n_lines))   || n_lines   <= 0) n_lines   <- 1
  max(min_in, n_terms * font_size * line_spacing * n_lines / 72 + overhead_in)
}

# Resolution for the PNG export.
#
# The bug being avoided: png(width, height, res) sets the PHYSICAL size to
# width/res inches. Fixing the pixels and choosing res separately silently
# changes the canvas size. Callers therefore pass INCHES and multiply by this
# dpi for both dimensions, so raising the dpi changes sharpness and nothing else.
#
# The cap exists because the height now grows with the number of terms: 100
# terms needs ~44in, which at 300dpi would be a 13000px-tall image. Reducing the
# dpi keeps the layout exactly as requested and only costs sharpness, which is
# the right thing to give up.
enrich_export_dpi <- function(w_in, h_in, target_dpi = 300, max_px = 10000) {
  w_in <- suppressWarnings(as.numeric(w_in)); h_in <- suppressWarnings(as.numeric(h_in))
  if (!isTRUE(is.finite(w_in)) || w_in <= 0) w_in <- 10
  if (!isTRUE(is.finite(h_in)) || h_in <= 0) h_in <- 8
  max(72, min(target_dpi, floor(max_px / max(w_in, h_in))))
}

# The export-size controls, shared by both plot tabs for the same reason
# enrichPlotStyleUI() is: two panels, one definition, no drift.
# prefix : lets one module carry more than one export block. The GO panel has
#          two (the dot plot and the term heatmap) and they share a namespace,
#          so their input ids must differ: "exp_w" vs "ht_exp_w".
# unit   : the noun in the help text and the note ("terms" / "genes").
enrichPlotExportUI <- function(ns, prefix = "exp", unit = "terms",
                               default_w = 10, default_h = 8) {
  id <- function(x) ns(paste0(prefix, "_", x))
  fluidRow(
    column(3, numericInput(
      id("w"),
      lbl("Export width (in)",
          paste0("The PDF/PNG is drawn on a canvas of exactly this size. Text is ",
                 "sized in points, so a narrower canvas makes the same font look ",
                 "bigger relative to the figure - this is what controls that.")),
      value = default_w, min = 3, max = 40, step = 0.5)),
    column(3, numericInput(
      id("h"),
      lbl("Export height (in)",
          paste0("Set automatically from the number of ", unit, " and the font ",
                 "size, because each ", sub("s$", "", unit), " needs a fixed amount ",
                 "of vertical space. Edit it freely; changing the ", unit,
                 " shown or the font size sets it again.")),
      value = default_h, min = 3, max = 80, step = 0.5)),
    column(6, uiOutput(id("note")))
  )
}

# Wires enrichPlotExportUI() up inside a module server: keeps the height field
# in step with the data, renders the note beside it, and returns the two
# download `content` functions so the GO and KEGG handlers are the same code.
#
# result_r : reactive returning the module's {table, ...} enrichment result.
# plot_r   : reactive returning the ggplot to write out.
# Returns list(pdf = function(file), png = function(file)).
# prefix     : which export block this is (see enrichPlotExportUI).
# needed_h_r : reactive giving the height in inches this figure needs. Supplied
#              by the caller for anything that is not the standard ORA dot plot
#              (the GO term heatmap sizes itself from GENES, not terms).
# draw_fun   : how to render onto the open device. print() for a ggplot,
#              ComplexHeatmap::draw() for a heatmap.
# n_items_r / font_r : only for the note's wording.
enrich_export_server <- function(input, output, session, result_r = NULL, plot_r = NULL,
                                 prefix = "exp", unit = "terms",
                                 needed_h_r = NULL, draw_fun = NULL,
                                 n_items_r = NULL, font_r = NULL) {
  id_w <- paste0(prefix, "_w"); id_h <- paste0(prefix, "_h")

  # Lines the longest label wraps onto, and hence the height the plot needs.
  # Used only when the caller does not supply its own needed_h_r.
  needed_h_default <- reactive({
    r   <- tryCatch(result_r(), error = function(e) NULL)
    tbl <- r$table
    if (is.null(tbl) || !"Description" %in% colnames(tbl) || nrow(tbl) == 0) return(NA_real_)
    n     <- min(as.integer(input$top_n %||% 20), nrow(tbl))
    wrap  <- suppressWarnings(as.numeric(input$label_wrap %||% 45))
    lines <- if (isTRUE(is.finite(wrap)) && wrap > 0) {
      max(vapply(as.character(tbl$Description),
                 function(s) length(strwrap(s, width = wrap)), integer(1)))
    } else 1L
    enrich_plot_height_in(n, input$plot_fs %||% 13, lines)
  })
  needed_h <- if (is.null(needed_h_r)) needed_h_default else needed_h_r

  # Follow the data. An observer rather than a default, because the right height
  # depends on Top N, the font size and the label wrap, all of which change
  # after the control is created.
  observeEvent(needed_h(), {
    h <- needed_h()
    if (isTRUE(is.finite(h))) updateNumericInput(session, id_h, value = round(h, 1))
  }, ignoreNULL = TRUE)

  output[[paste0(prefix, "_note")]] <- renderUI({
    h <- needed_h()
    if (!isTRUE(is.finite(h))) return(NULL)
    w   <- suppressWarnings(as.numeric(input[[id_w]] %||% 10))
    cur <- suppressWarnings(as.numeric(input[[id_h]] %||% h))
    dpi <- enrich_export_dpi(w, cur)
    short <- isTRUE(is.finite(cur)) && cur < h - 0.05
    n_it <- if (!is.null(n_items_r)) tryCatch(n_items_r(), error = function(e) NA) else
      tryCatch(min(as.integer(input$top_n %||% 20), nrow(result_r()$table)),
               error = function(e) NA)
    fs   <- if (!is.null(font_r)) tryCatch(font_r(), error = function(e) NA) else
      (input$plot_fs %||% 13)
    tags$div(
      style = "font-size:0.85em; margin-top:26px;",
      tags$span(style = if (short) "color:#a94442;" else "color:#3c763d;",
                sprintf("%s %.1f in of height for %s %s at %s pt.",
                        if (short) "⚠ Needs" else "✓ Fits:", h,
                        if (is.na(n_it)) "these" else format(n_it), unit, format(fs))),
      if (short) tags$span(" Labels will overlap at the current height."),
      tags$br(),
      tags$span(style = "color:#777;",
                sprintf("PNG will be %.0f x %.0f px at %d dpi.", w * dpi, cur * dpi, dpi))
    )
  })

  size <- function() {
    w <- suppressWarnings(as.numeric(input[[id_w]] %||% 10))
    h <- suppressWarnings(as.numeric(input[[id_h]] %||% NA_real_))
    if (!isTRUE(is.finite(w)) || w <= 0) w <- 10
    # Falls back to the COMPUTED height, not to a constant. The observer above
    # normally has the field already set, but that is a round-trip through the
    # browser: if the user downloads before it lands, or clears the box, the
    # export must still be tall enough for the terms rather than silently
    # reverting to a fixed 8in and putting the labels back on top of each other.
    if (!isTRUE(is.finite(h)) || h <= 0) {
      h <- tryCatch(needed_h(), error = function(e) NA_real_)
      if (!isTRUE(is.finite(h)) || h <= 0) h <- 8
    }
    list(w = w, h = h)
  }

  # print() for a ggplot; ComplexHeatmap::draw() for a heatmap.
  render <- if (is.null(draw_fun)) function() print(plot_r()) else draw_fun

  list(
    pdf = function(file) {
      s <- size()
      grDevices::pdf(file, width = s$w, height = s$h)
      on.exit(grDevices::dev.off(), add = TRUE)
      render()
    },
    png = function(file) {
      s   <- size()
      dpi <- enrich_export_dpi(s$w, s$h)
      # Pixels AND res scale together, so the canvas stays exactly s$w x s$h
      # inches and the dpi only decides sharpness.
      grDevices::png(file, width = s$w * dpi, height = s$h * dpi, res = dpi)
      on.exit(grDevices::dev.off(), add = TRUE)
      render()
    }
  )
}

# Reads the enrichPlotStyleUI() controls out of a module's `input` and returns
# them as the named list make_enrichment_plot() expects. Every field falls back
# to that function's own default, so a module can call this before the panel has
# ever been rendered (Shiny returns NULL for an input that does not exist yet)
# and still get the figure it would have drawn without the panel.
enrich_plot_style_args <- function(input) {
  list(
    style         = input$plot_style    %||% "dot",
    x_metric      = input$x_metric      %||% "GeneRatio",
    color_mode    = input$color_mode    %||% "gradient",
    color_metric  = input$color_metric  %||% "p.adjust",
    color_low     = input$color_low     %||% "#3B7DD8",
    color_high    = input$color_high    %||% "#C0392B",
    flat_color    = input$flat_color    %||% "#3B7DD8",
    size_by_count = if (is.null(input$size_by_count)) TRUE else isTRUE(input$size_by_count),
    theme_style   = input$theme_style   %||% "bw",
    label_wrap    = input$label_wrap    %||% 45
  )
}

# ── Freeze the figure's appearance AT RUN TIME ───────────────────────────────
# The style controls are LIVE: the plot re-renders as you drag a slider, so
# there is no "the figure was built now" moment the way there is for the
# heatmap. Reading them when the reproducibility bundle is downloaded would
# therefore record whatever the panel happened to say at download time, which
# is precisely the drift heatmap_history was restructured to eliminate.
#
# Instead this is called INSIDE goServer()/keggServer()'s eventReactive(input$run)
# body, so the values recorded are the ones that were on screen when Run was
# pressed. eventReactive isolates its body from everything but its event
# expression, so reading the style inputs here creates no extra dependency and
# cannot re-fire the run.
#
# Returns the full make_enrichment_plot() argument set (style + framing), which
# is what write_rerun_GO()/write_rerun_KEGG() need to emit a runnable figure.
# `title` is passed in because only the caller knows it (GO includes the
# ontology, KEGG is fixed).
enrich_plot_snapshot <- function(input, title) {
  c(enrich_plot_style_args(input),
    list(
      title     = title,
      top_n     = input$top_n     %||% 20,
      pinned    = as.character(input$pin_terms %||% character(0)),
      font_size = input$plot_fs   %||% 13,
      sort_desc = if (is.null(input$sort_desc)) TRUE else isTRUE(input$sort_desc)
    ))
}

# Attach a run-time style snapshot to an enrichment result's params block.
# No-op when the result carries no params (nothing to attach it to), so a
# malformed result can never turn into an error at the call site.
#
# The result object is COPIED before it is modified, which matters on a cache
# hit: run_go_enrichment() may return the very object stored in .enrich_cache,
# and writing the current run's style onto it would retro-fit that style to
# every earlier run that shares the cache key.
attach_plot_snapshot <- function(res, input, title) {
  if (!is.list(res) || !is.list(res$params)) return(res)
  res$params$plot <- enrich_plot_snapshot(input, title)
  res
}
