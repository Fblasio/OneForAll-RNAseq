# =============================================================================
# tair_module.R
# Shiny module: selected genes table with ThaleMine, UniProt and Perplexity links
# Receives: selected_genes  (reactiveVal with gene ID vector)
#           gene_labels_rv  (reactiveVal with named vector: names=ID, value=symbol)
# =============================================================================

# ── Helper: batch UniProt lookup via REST API ─────────────────────────────────
# Queries UniProt using gene_exact terms restricted to the active species'
# organism_id (UNIPROT_TAXID). IMPORTANT: `gene_ids` must already be in the
# form UniProt/SGN index — for tomato that means LOCUS ONLY, no ITAG version
# suffix (Solyc04g072930, not Solyc04g072930.1); the caller passes
# link_gene_id()-derived IDs for exactly this reason.
# Matching against UniProt results is done case-insensitively (UniProt's
# orderedLocusNames often come back uppercased), but results are keyed by the
# ORIGINAL input casing so downstream lookups stay consistent with the app.
# Returns a named list keyed by gene_ids: list(protein_name, accession)
.fetch_uniprot_batch <- function(gene_ids) {
  if (length(gene_ids) == 0) return(list())
  if (!requireNamespace("httr",      quietly = TRUE)) return(list())
  if (!requireNamespace("jsonlite",  quietly = TRUE)) return(list())
  
  gene_ids_q  <- trimws(as.character(gene_ids))
  upper_to_in <- stats::setNames(gene_ids_q, toupper(gene_ids_q))  # UPPER -> original casing

  # ── Tomato needs a different query AND a different matching key ────────────
  # Verified against the live REST API: for organism_id 4081, `gene_exact:` on a
  # Solyc locus returns ZERO hits in every form (with version, without version,
  # locus only) -- UniProt does not index Solyc IDs under the gene-name field.
  # The gene names it does hold are things like LOC778315 / TIP1;1, so the
  # orderedLocusNames matching below never fires either. Net effect before this
  # fix: every tomato row fell back to a bare "Search UniProt" link.
  #
  # What does work is a free-text query, with the Solyc ID recovered from the
  # entry's EnsemblPlants cross-reference (e.g. "Solyc06g074820.3.1").
  # Free text also returns unrelated entries, so the xref is not a nicety --
  # it is what makes the match trustworthy.
  is_tomato <- identical(ACTIVE_SPECIES, "tomato")

  # Gene-level key for tomato: the bare locus. UniProt's EnsemblPlants xref is
  # transcript-level ("Solyc06g074820.3.1") and the count matrix may carry a
  # different ITAG version (".2" vs ".3"), so both sides are reduced to the
  # locus -- no dot, no version, no transcript.
  .locus <- function(x) toupper(sub("^(Solyc[0-9]{2}g[0-9]+).*$", "\\1", trimws(x),
                                    ignore.case = TRUE))

  # One REST call. Returns the parsed `results` list, or NULL.
  .fetch <- function(query, extra_fields = "") {
    url <- paste0(
      "https://rest.uniprot.org/uniprotkb/search?",
      "query=",  utils::URLencode(query, reserved = TRUE),
      "&fields=accession,gene_names,protein_name", extra_fields,
      "&format=json&size=500"
    )
    tryCatch({
      resp <- httr::GET(url, httr::timeout(15))
      if (httr::http_error(resp)) return(NULL)
      jsonlite::fromJSON(httr::content(resp, "text", encoding = "UTF-8"),
                         simplifyDataFrame = FALSE)$results
    }, error = function(e) {
      message("UniProt batch fetch error: ", e$message); NULL
    })
  }

  .xref <- function(entry, db) {
    tryCatch(unlist(lapply(entry$uniProtKBCrossReferences, function(x) {
      if (identical(x$database, db)) x$id else NULL
    })), error = function(e) character(0))
  }

  result <- list()
  # First match wins. That is what gives the Entrez pass priority over the
  # free-text pass below, and it makes the answer independent of how many
  # entries UniProt happens to return in one response.
  .add <- function(gid, pname, acc) {
    if (is.null(result[[gid]]))
      result[[gid]] <<- list(
        protein_name = if (!is.na(pname)) pname else "Unknown protein",
        accession    = acc)
  }

  .pname <- function(entry) {
    p <- tryCatch(entry$proteinDescription$recommendedName$fullName$value,
                  error = function(e) NULL)
    # TrEMBL entries have no recommendedName -- they carry submissionNames
    # instead. Most tomato hits are TrEMBL, so without this fallback a large
    # share of rows showed "Unknown protein" while UniProt actually holds a
    # usable description (e.g. A0A3Q7H0N6 -> "Similar to aquaporin").
    if (is.null(p) || length(p) == 0)
      p <- tryCatch(entry$proteinDescription$submissionNames[[1]]$fullName$value,
                    error = function(e) NULL)
    if (is.null(p) || length(p) == 0) NA_character_ else p
  }

  # ══ TOMATO: Entrez cross-reference first, Solyc free text as the fallback ══
  # Pass 1 uses data the app already holds -- data/id_map_tomato.csv, loaded as
  # ENTREZ_MAP for KEGG. `xref:geneid-<id>` is a FIELD-RESTRICTED query, so it
  # cannot pick up the free-text false positives pass 2 has to filter out, and
  # it tends to surface the better-curated entry of a locus.
  # Pass 2 covers the ~25% of genes with no Entrez ID in the map (measured on a
  # real selection: 16 of 20 had one), where the only usable handle is the
  # Solyc ID -- which UniProt files under EnsemblPlants, not as a gene name.
  # Together they resolved 19 of 20 real genes; either alone resolved fewer.
  if (is_tomato) {
    ent_for <- rep(NA_character_, length(gene_ids_q))
    if (exists("ENTREZ_MAP") && !is.null(ENTREZ_MAP) && nrow(ENTREZ_MAP) > 0) {
      ent_for <- ENTREZ_MAP$ENTREZID[match(.locus(gene_ids_q), .locus(ENTREZ_MAP[[1]]))]
    }
    has_ent <- !is.na(ent_for) & nzchar(ent_for)

    if (any(has_ent)) {
      q1 <- paste0("(", paste(paste0("xref:geneid-", unique(ent_for[has_ent])),
                              collapse = " OR "), ")",
                   " AND organism_id:", UNIPROT_TAXID)
      for (entry in .fetch(q1, ",xref_geneid") %||% list()) {
        gids <- .xref(entry, "GeneID")
        pn   <- .pname(entry)
        for (k in which(has_ent & ent_for %in% gids))
          .add(gene_ids_q[k], pn, entry$primaryAccession)
      }
    }

    todo <- setdiff(gene_ids_q, names(result))
    if (length(todo) > 0) {
      q2 <- paste0("(", paste(unique(todo), collapse = " OR "), ")",
                   " AND organism_id:", UNIPROT_TAXID)
      by_locus <- split(todo, .locus(todo))
      for (entry in .fetch(q2, ",xref_ensemblplants") %||% list()) {
        pn <- .pname(entry)
        for (loc in intersect(names(by_locus), .locus(.xref(entry, "EnsemblPlants"))))
          for (gid in by_locus[[loc]]) .add(gid, pn, entry$primaryAccession)
      }
    }
    return(result)
  }

  # ══ ARABIDOPSIS: unchanged -- gene_exact: + orderedLocusNames ══════════════
  query <- paste0("(", paste(paste0("gene_exact:", gene_ids_q), collapse = " OR "), ")",
                  " AND organism_id:", UNIPROT_TAXID)
  url <- paste0(
    "https://rest.uniprot.org/uniprotkb/search?",
    "query=",  utils::URLencode(query, reserved = TRUE),
    "&fields=accession,gene_names,protein_name",
    "&format=json&size=500"
  )

  tryCatch({
    resp <- httr::GET(url, httr::timeout(15))
    if (httr::http_error(resp)) return(list())

    parsed  <- jsonlite::fromJSON(
      httr::content(resp, "text", encoding = "UTF-8"),
      simplifyDataFrame = FALSE
    )
    entries <- parsed$results
    if (is.null(entries) || length(entries) == 0) return(list())

    result <- list()
    for (entry in entries) {
      acc   <- entry$primaryAccession
      pname <- .pname(entry)

      {
        # Match genes via orderedLocusNames (locus IDs stored here in UniProt)
        oll_ids <- tryCatch({
          all_oll <- lapply(entry$genes, function(g) {
            if (!is.null(g$orderedLocusNames))
              vapply(g$orderedLocusNames, function(x) x$value, character(1))
            else
              character(0)
          })
          toupper(unlist(all_oll))
        }, error = function(e) character(0))

        for (gid_upper in intersect(names(upper_to_in), oll_ids)) {
          gid <- upper_to_in[[gid_upper]]
          result[[gid]] <- list(
            protein_name = if (!is.na(pname)) pname else "Unknown protein",
            accession    = acc
          )
        }
      }
    }
    result
  }, error = function(e) {
    message("UniProt batch fetch error: ", e$message)
    list()
  })
}

# ── Helper: build a Perplexity search URL for a given gene ID ─────────────────
# The prompt is fixed except for the gene ID, which replaces [GENE_NAME].
# thalemine_url and uniprot_url are appended as context to help Perplexity find the gene.
.make_perplexity_url <- function(gene_id, thalemine_url, uniprot_url) {
  prompt <- paste0(
    "You are an expert molecular biologist and scientific literature analyst.\n\n",
    "Your task is to compile a structured, evidence-based and DETAILED report on the following gene: ",
    gene_id, ".\n\n",
    "\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\n",
    "\u26a0\ufe0f STRICT CONSTRAINT (non-negotiable):\n",
    "You MUST only report information explicitly documented in\n",
    "peer-reviewed scientific articles. Do NOT generate, infer,\n",
    "assume, or extrapolate any data beyond what is directly\n",
    "stated in the published literature. If information for a\n",
    "specific section is unavailable in the scientific record,\n",
    "write exactly: \"No peer-reviewed data available.\"\n",
    "\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\n\n",
    "Organize your response using the following sections:\n\n",
    "## 1. Official Gene Name\n",
    "Provide the official gene symbol and full name (e.g., HGNC-approved nomenclature or equivalent database).\n\n",
    "## 2. Alternative Names & Aliases\n",
    "List all known aliases found in the literature, including:\n",
    "- Alternative gene symbols or IDs\n",
    "- Encoded protein name(s) and isoforms\n",
    "- Any other designations used across publications\n\n",
    "## 3. Biological Functions & Pathways\n",
    "Provide a DETAILED description of biological roles, molecular functions, and cellular or signaling pathways ",
    "this gene/protein participates in. Only include findings explicitly reported in peer-reviewed ",
    "sources. Cite the article for each claim.\n\n",
    "## 4. Supporting Scientific References\n",
    "For every piece of information above, provide the full citation in this format:\n",
    "> [Author(s), Year]. Title. *Journal Name*. DOI: xxxxx\n\n",
    "\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\n",
    "REMINDER: Every single claim in sections 3 and 4 must map to at least one citation. ",
    "If you cannot cite it, do not include it.\n",
    "\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\n",
    # Context links appended after the prompt to help Perplexity locate the gene
    "\nThis information can help you:\n",
    paste0("- ", GENE_DB_LABEL, ": "), thalemine_url, "\n"
  )
  paste0("https://www.perplexity.ai/search?q=",
         utils::URLencode(prompt, reserved = TRUE))
}

# ── UI ───────────────────────────────────────────────────────────────────────
tairGenesUI <- function(id) {
  ns <- NS(id)
  tagList(
    br(),
    fluidRow(
      column(12,
        helpText(HTML(paste0("🔗 Click on a Gene ID to open its ", GENE_DB_LABEL, " page. UniProt column fetched live — may take a moment."))),
        downloadButton(ns("dl_csv"), "\u2b07 Download CSV",
                       class = "btn-sm btn-success"),
        br(), br(),
        DTOutput(ns("table"))
      )
    )
  )
}

# ── Server ───────────────────────────────────────────────────────────────────
tairGenesServer <- function(id, selected_genes, gene_labels_rv) {
  moduleServer(id, function(input, output, session) {

    # Parse name to build the table if the active species = TOMATO.
    # For tomato this also doubles as the UniProt/SGN QUERY ID: both
    # services index the locus only, without the ITAG version suffix.
    link_gene_id <- function(x) {
      x <- trimws(as.character(x))
      if (identical(ACTIVE_SPECIES, "tomato")) {
        # ignore.case = TRUE: gene IDs may arrive uppercase (e.g.
        # SOLYC04G005160.3.1) depending on the source matrix. With
        # ignore.case = FALSE this regex silently failed to match and returned
        # the ID unchanged, so the UniProt/SGN query carried the full versioned
        # transcript ID and matched nothing — every row fell back to a bare
        # "Search UniProt" link. Callers also canon_ids() first now, so this is
        # belt-and-braces rather than the only defence.
        sub("^((Solyc[0-9]{2}g[0-9]+)).*$", "\\1", x, ignore.case = TRUE)
      } else {
        x
      }
    }
    
    uniprot_rv <- reactive({
      genes <- canon_ids(selected_genes())
      req(genes)
      query_ids <- vapply(genes, link_gene_id, character(1))
      withProgress(message = "Fetching UniProt protein names...", value = 0.5, {
        .fetch_uniprot_batch(query_ids)
      })
    })

    # Build the full dataframe with HTML links for the DT and plain text for CSV
    genes_df <- reactive({
      req(selected_genes())
      # canon_ids() so the displayed IDs, the SGN/ThaleMine links and the
      # UniProt lookups all use the same canonical form the rest of the app
      # uses (heatmap subsetting, GO/KEGG/STRING). Without it the raw matrix
      # IDs were shown verbatim — e.g. SOLYC04G005160.3.1 instead of
      # Solyc04g005160.3 — and every external query was built from a string
      # those services do not index.
      genes    <- canon_ids(selected_genes())
      link_ids <- vapply(genes, link_gene_id, character(1))
      gl       <- gene_labels_rv()
      uni_info <- tryCatch(uniprot_rv(), error = function(e) list())

      # Symbols come straight from the org.db so they show regardless of the
      # heatmap "Map gene IDs to gene symbols" toggle.
      # Fallbacks: heatmap labels, then the raw ID.
      symbols <- genes
      direct  <- tryCatch(map_gene_symbols(genes), error = function(e) NULL)
      if (!is.null(direct)) {
        sym     <- as.character(direct[normalize_gene_ids(genes)])
        symbols <- ifelse(is.na(sym) | !nzchar(sym), genes, sym)
      } else if (!is.null(gl) && length(gl) > 0) {
        idx     <- match(normalize_gene_ids(genes), normalize_gene_ids(names(gl)))
        symbols <- ifelse(is.na(idx), genes, as.character(gl[idx]))
      }

      # ThaleMine links (replaces TAIR)
      thalemine_url   <- paste0(GENE_DB_URL, link_ids)
      thalemine_links <- paste0('<a href="', thalemine_url, '" target="_blank">', genes, '</a>')

      # UniProt protein name + direct entry link (or fallback search link)
      prot_html <- vapply(seq_along(genes), function(i) {
        
        gid    <- genes[i]
        gidkey <- link_ids[i]
        info   <- uni_info[[gidkey]]
        
        if (!is.null(info)) {
          entry_url <- paste0("https://www.uniprot.org/uniprotkb/", info$accession, "/entry")
          pname <- gsub('"', "&quot;", info$protein_name)
          paste0('<a href="', entry_url, '" target="_blank">', pname, '</a>')
        } else {
          search_url <- if (identical(ACTIVE_SPECIES, "tomato")) {
            paste0(
              "https://www.uniprot.org/uniprotkb?query=",
              utils::URLencode(gidkey, reserved = TRUE)
            )
          } else {
            paste0(
              "https://www.uniprot.org/uniprotkb?query=gene_exact%3A",
              utils::URLencode(gid, reserved = TRUE),
              "%20AND%20organism_id%3A",
              UNIPROT_TAXID
            )
          }
          paste0('<a href="', search_url, '" target="_blank">Search UniProt</a>')
        }
      }, character(1))

      # Plain-text protein name for CSV export
      prot_plain <- vapply(seq_along(genes), function(i) {
        info <- uni_info[[link_ids[i]]]
        if (!is.null(info)) info$protein_name else NA_character_
      }, character(1))
      
      prot_accession <- vapply(seq_along(genes), function(i) {
        info <- uni_info[[link_ids[i]]]
        if (!is.null(info)) info$accession else NA_character_
      }, character(1))

      # Perplexity links — one per gene, with ThaleMine and UniProt URLs as context
      perplexity_html <- vapply(seq_along(genes), function(i) {
        gid    <- genes[i]
        gidkey <- link_ids[i]
        tm_url <- thalemine_url[i]
        info   <- uni_info[[gidkey]]
        uurl <- if (!is.null(info)) {
          paste0("https://www.uniprot.org/uniprotkb/", info$accession, "/entry")
        } else if (identical(ACTIVE_SPECIES, "tomato")) {
          paste0(
            "https://www.uniprot.org/uniprotkb?query=",
            utils::URLencode(link_ids[i], reserved = TRUE)
          )
        } else {
          paste0(
            "https://www.uniprot.org/uniprotkb?query=gene_exact%3A",
            utils::URLencode(gidkey, reserved = TRUE),
            "%20AND%20organism_id%3A",
            UNIPROT_TAXID
          )
        }
        url <- .make_perplexity_url(gid, tm_url, uurl)
        paste0('<a href="', url, '" target="_blank">\U0001f50d Ask Perplexity</a>')
      }, character(1))

      # Store both HTML version (for DT) and plain version (for CSV) in same df
      data.frame(
        `Gene ID`           = thalemine_links,
        `Symbol`            = symbols,
        `Protein (UniProt)` = prot_html,
        `Ask Perplexity`    = perplexity_html,
        # Plain-text columns used only for CSV (not shown in DT)
        `_gene_id_plain`    = genes,
        `_thalemine_url`    = thalemine_url,
        `_prot_plain`       = prot_plain,
        `_prot_accession`   = prot_accession,
        stringsAsFactors    = FALSE,
        check.names         = FALSE
      )
    })

    # DT table — renders HTML columns, hides the plain-text helper columns
    output$table <- renderDT({
      req(genes_df())
      df      <- genes_df()
      display <- df[, c("Gene ID", "Symbol", "Protein (UniProt)", "Ask Perplexity")]
      datatable(
        display,
        escape    = FALSE,
        rownames  = FALSE,
        selection = "none",
        options   = list(
          pageLength  = 25,
          scrollX     = TRUE,
          columnDefs  = list(
            list(className = "dt-left", targets = "_all")
          )
        )
      )
    })

    # CSV download — exports plain-text columns, no HTML
    output$dl_csv <- downloadHandler(
      filename = function() {
        paste0("selected_genes_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
      },
      content = function(f) {
        df <- genes_df()
        export_df <- data.frame(
          `Gene ID`           = df[["_gene_id_plain"]],
          `Symbol`            = df[["Symbol"]],
          `Protein name`      = df[["_prot_plain"]],
          `UniProt accession` = df[["_prot_accession"]],
          `Gene DB URL`   = df[["_thalemine_url"]],
          stringsAsFactors    = FALSE,
          check.names         = FALSE
        )
        write.csv(export_df, f, row.names = FALSE)
      }
    )

  })
}
