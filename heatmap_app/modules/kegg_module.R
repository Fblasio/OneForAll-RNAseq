keggControlsUI <- function(id) {
  ns <- NS(id)
  tagList(
    h5("KEGG Settings"),
    # Same wording as goControlsUI()'s p.adjust help on purpose: it is the same
    # control with the same trap, and two different explanations of one concept
    # would be worse than none.
    selectInput(ns("padj_method"),
                lbl("p.adjust method", "BH controls the false discovery rate - the default, and the right choice here. Bonferroni / holm / hochberg / hommel are far stricter. 'none' means NO correction at all; never report those p-values as significant."),
                choices  = c("BH","BY","bonferroni","holm","hochberg","hommel","none"),
                selected = "BH"),
    numericInput(ns("pval"), "p-value cutoff",
                 value = 0.05, min = 0, max = 1, step = 0.01),
    numericInput(ns("qval"), "q-value cutoff",
                 value = 0.05, min = 0, max = 1, step = 0.01),
    checkboxInput(ns("use_cache"), "Use cache", value = TRUE),
    actionButton(ns("run"), "\u25b6 Run KEGG",
                 class = "btn-warning btn-sm", width = "100%")
  )
}

keggTableUI <- function(id) {
  ns <- NS(id)
  tagList(
    br(),
    downloadButton(ns("dl_table"), "\u2b07 Download CSV", class = "btn-sm btn-success"),
    br(), br(),
    helpText("\U0001f4a1 Click a row \u2192 generate pathway map."),
    DTOutput(ns("table")),
    uiOutput(ns("pathview_panel"))
  )
}

keggPlotUI <- function(id) {
  ns <- NS(id)
  tagList(
    br(),
    fluidRow(
      column(6, sliderInput(ns("plot_h"),  "Plot height (px)", 300, 1600, 800, 50, ticks=FALSE)),
      column(6, sliderInput(ns("plot_fs"), "Font size",        6,   28,  13,  1,  ticks=FALSE))
    ),
    fluidRow(
      column(6, sliderInput(ns("top_n"), "Top N KEGG pathways", 5, 100, 20, 5, ticks=FALSE)),
      column(6, checkboxInput(ns("sort_desc"), "Sort highest first", value = TRUE))
    ),
    fluidRow(
      column(12, selectizeInput(ns("pin_terms"), "Always show these pathways (pinned):",
                                choices = NULL, multiple = TRUE,
                                options = list(placeholder = "Select KEGG pathways to pin...",
                                               plugins = list("remove_button"))))
    ),
    # ── Plot style / customization panel ──────────────────────────────────
    # The same controls the GO Plot tab offers, from the one shared definition
    # in utils/enrichment_functions.R, so the two panels cannot diverge.
    enrichPlotStyleUI(ns, terms_label = "pathways"),
    # Same shared export-size block as the GO Plot tab -- one definition, so the
    # two tabs cannot drift apart. See utils/enrichment_functions.R.
    enrichPlotExportUI(ns),
    fluidRow(
      column(4, downloadButton(ns("dl_pdf"), "\u2b07 PDF", class = "btn-sm btn-success")),
      column(4, downloadButton(ns("dl_png"), "\u2b07 PNG", class = "btn-sm btn-info"))
    ),
    br(),
    uiOutput(ns("plot_container"))
  )
}

keggServer <- function(id, selected_genes, scaled_matrix_rv = NULL,
                       deg_genes_rv = NULL, universe = NULL) {
  moduleServer(id, function(input, output, session) {
    
    pathview_png_path <- reactiveVal(NULL)
    
    result <- eventReactive(input$run, {
      req(selected_genes())
      # Same run-time style freeze as goServer() -- see attach_plot_snapshot()
      # in utils/enrichment_functions.R for why it happens here and not at
      # bundle-download time.
      attach_plot_snapshot(withProgress(message = "Running KEGG enrichment...", {
        withCallingHandlers(
          tryCatch(
            run_kegg_enrichment(
              selected_genes(),
              padj_method = input$padj_method,
              pval_cutoff = input$pval,
              qval_cutoff = input$qval,
              use_cache   = input$use_cache,
              # Mirrors the goServer() wiring: NULL stays whole-genome, so any
              # caller that does not pass a universe behaves exactly as before.
              universe    = if (is.null(universe)) NULL else universe()
            ),
            error = function(e) {
              # Same reasoning as go_module.r's error branch: without params the
              # history observer drops the run, so a failed attempt left no
              # trace anywhere the user could see.
              log_error(paste0("KEGG run failed: ", conditionMessage(e)))
              list(table  = data.frame(Error = conditionMessage(e)),
                   object = NULL,
                   params = list(species       = ACTIVE_SPECIES,
                                 kegg_code     = KEGG_CODE,
                                 padj_method   = input$padj_method,
                                 pval_cutoff   = input$pval,
                                 qval_cutoff   = input$qval,
                                 n_genes_input = tryCatch(length(unique(canon_ids(selected_genes()))),
                                                          error = function(e2) NA_integer_),
                                 failed        = TRUE,
                                 error         = conditionMessage(e)))
            }
          ),
          warning = function(w) {
            message("KEGG WARNING: ", w$message)
            invokeRestart("muffleWarning")
          }
        )
      }), input, "KEGG Pathway Enrichment")
    })
    
    observeEvent(result(), {
      for (w in result()$warnings %||% character(0)) {
        showNotification(paste0("\u26a0\ufe0f ", w), type = "warning", duration = 14)
        log_warning(w)
      }
      if (!is.null(result()$reason)) {
        showNotification(paste0("\u2139\ufe0f ", result()$reason),
                         type = "message", duration = 14)
      }
    })

    output$table <- renderDT({
      req(result())
      tbl <- result()$table
      if ("ID" %in% colnames(tbl)) {
        tbl$ID <- paste0(
          '<a href="https://www.kegg.jp/pathway/', tbl$ID,
          '" target="_blank">', tbl$ID, '</a>'
        )
      }
      datatable(tbl,
                escape    = FALSE,
                selection = "single",
                options   = list(pageLength = 10, scrollX = TRUE),
                rownames  = FALSE)
    })
    
    output$pathview_panel <- renderUI({
      sel <- input$table_rows_selected
      req(sel)
      tbl    <- result()$table
      raw_id <- gsub("<[^>]+>", "", tbl[sel, "ID"])
      
      mat <- if (!is.null(scaled_matrix_rv))
        tryCatch(scaled_matrix_rv(), error = function(e) NULL) else NULL
      grp_choices <- if (!is.null(mat)) colnames(mat) else character(0)
      
      tagList(
        hr(),
        fluidRow(
          column(5,
                 tags$p(style = "margin-top: 6px;",
                        tags$b("Selected: "),
                        tags$span(tbl[sel, "Description"], style = "color: #1F3864;"),
                        tags$span(paste0(" (", raw_id, ")"), style = "color:gray; font-size:0.85em;"))
          ),
          column(4,
                 selectizeInput(session$ns("pathview_groups"),
                                label    = "Groups to display (1 or 2):",
                                choices  = grp_choices,
                                selected = grp_choices[1],
                                multiple = TRUE,
                                options  = list(maxItems = 2,
                                                placeholder = "Select 1 or 2 groups..."))
          ),
          column(3,
                 br(),
                 if (grepl("Global and overview", tbl[sel, "Description"], ignore.case = TRUE) ||
                     grepl("01200|01230|01240|01232", raw_id)) {
                   tags$p(style = "color: #e67e22; font-size: 0.85em; margin-top: 6px;",
                          icon("exclamation-triangle"),
                          "Global overview maps cannot be rendered by pathview.")
                 } else {
                   actionButton(session$ns("run_pathview"),
                                "\U0001f5fa\ufe0f Generate Pathway Map",
                                class = "btn-info btn-sm", width = "100%")
                 }
          )
        )
      )
    })
    
    observeEvent(input$run_pathview, {
      req(selected_genes(), result())
      sel <- input$table_rows_selected
      req(sel)
      
      tbl      <- result()$table
      raw_id   <- gsub("<[^>]+>", "", tbl[sel, "ID"])
      path_num <- gsub("^[a-z]+", "", raw_id)
      species  <- gsub("[0-9]+$", "", raw_id)
      if (!nzchar(species)) species <- KEGG_CODE
      if (!nzchar(path_num)) path_num <- raw_id
      
      # Helper: wait for PNG to be fully written, then embed as base64
      show_pathview_modal <- function(png_path) {
        # Poll until file has content (max 10 seconds)
        png_ready <- FALSE
        for (i in 1:20) {
          Sys.sleep(0.5)
          if (file.exists(png_path) && file.info(png_path)$size > 0) {
            png_ready <- TRUE
            break
          }
        }
        message(">>> PNG ready: ", png_ready, " | size: ", file.info(png_path)$size)
        
        if (!png_ready) {
          showNotification(
            paste0("\u26a0\ufe0f PNG not ready for ", raw_id, ". Try again."),
            type = "warning", duration = 8)
          return()
        }
        
        raw_bytes  <- readBin(png_path, "raw", file.info(png_path)$size)
        img_base64 <- jsonlite::base64_enc(raw_bytes)
        img_src    <- paste0("data:image/png;base64,", img_base64)
        pathview_png_path(png_path)
        showModal(modalDialog(
          title     = paste0("\U0001f5fa\ufe0f ", tbl[sel, "Description"],
                             " (", raw_id, ")"),
          size      = "l",
          easyClose = TRUE,
          footer    = tagList(
            downloadButton(session$ns("dl_pathview"),
                           "\u2b07 Download PNG", class = "btn-sm btn-success"),
            modalButton("Close")
          ),
          tags$img(src = img_src, style = "width:100%; height:auto;")
        ))
      }
      
      withProgress(message = "Generating pathway map...", {
        tryCatch({
          tmp_dir <- file.path(tempdir(), paste0("pv_", as.integer(Sys.time())))
          dir.create(tmp_dir, showWarnings = FALSE)
          old_wd <- setwd(tmp_dir)          # pathview writes output PNG to getwd()
          on.exit(setwd(old_wd), add = TRUE) # kegg.dir only caches templates
          
          mat <- if (!is.null(scaled_matrix_rv))
            tryCatch(scaled_matrix_rv(), error = function(e) NULL) else NULL
          
          # Extract genes directly from KEGG enrichment result for this pathway
          kegg_tbl  <- result()$table
          row_match <- which(gsub("<[^>]+>", "", kegg_tbl$ID) == raw_id)
          
          genes_raw <- if (length(row_match) > 0 && "geneID" %in% colnames(kegg_tbl)) {
            trimws(unlist(strsplit(as.character(kegg_tbl[row_match[1], "geneID"]), "/")))
          } else {
            as.character(selected_genes())
          }
          
          # ── Reconcile two different ID spaces ────────────────────────────
          # For non-Arabidopsis, run_kegg_enrichment() sends Entrez IDs and
          # keeps readable = FALSE, so kegg_tbl$geneID holds ENTREZ IDs. The
          # heatmap matrix (scaled_matrix_rv) is keyed on Solyc gene IDs. Those
          # never match, so the matrix lookup below found nothing and gene_data
          # silently stayed at rep(1) — a pathway map with every box flat.
          # So: translate Entrez -> Solyc here for the matrix lookup, and keep
          # gene_to_kegg to translate back before pathview(), which expects the
          # Entrez space implied by gene.idtype = "KEGG".
          # No ENTREZ_MAP (or Arabidopsis) -> unchanged previous behaviour.
          gene_to_kegg <- NULL
          if (!identical(ACTIVE_SPECIES, "arabidopsis") &&
              exists("ENTREZ_MAP") && !is.null(ENTREZ_MAP)) {
            hit <- match(genes_raw, as.character(ENTREZ_MAP$ENTREZID))
            solyc <- canon_ids(ifelse(is.na(hit), genes_raw,
                                      as.character(ENTREZ_MAP[[1]][hit])))
            genes_to_use <- solyc
            ok <- !is.na(hit)
            if (any(ok)) {
              # Solyc -> Entrez, first Entrez per gene (matches the
              # multiVals = "first" convention used elsewhere in the app).
              gene_to_kegg <- stats::setNames(genes_raw[ok], solyc[ok])
              gene_to_kegg <- gene_to_kegg[!duplicated(names(gene_to_kegg))]
            }
          } else {
            genes_to_use <- canon_ids(genes_raw)
          }
          
          if (length(genes_to_use) == 0) {
            showNotification("\u26a0\ufe0f No genes available for this pathway.",
                             type = "warning", duration = 6)
            return()
          }
          
          selected_grps <- input$pathview_groups
          # An explicit check rather than a bare req(), which aborts
          # silently. In the DEG panel keggServer() is wired without
          # scaled_matrix_rv, so the group selector has no choices and
          # input$pathview_groups stays NULL. A req() here would fire on every
          # click and the button would do nothing at all -- no error, no
          # message, no spinner. State the reason instead.
          if (length(selected_grps) < 1) {
            showNotification(
              paste0("ℹ️ Pathway maps colour genes by expression, so they need an ",
                     "expression matrix. This panel was opened from the DE results, which ",
                     "carry gene IDs only — open KEGG from the Heatmap panel to render ",
                     "pathway maps."),
              type = "warning", duration = 10)
            return()
          }
          
          # Build gene_data from scaled matrix; default to rep(1) if no match
          gene_data <- setNames(rep(1, length(genes_to_use)), genes_to_use)
          
          if (!is.null(mat) && length(genes_to_use) > 0) {
            genes_in <- ids_in(rownames(mat), genes_to_use)
            sub      <- mat[genes_in, , drop = FALSE]
            if (nrow(sub) > 0) {
              cols_ok <- intersect(selected_grps, colnames(sub))
              if (length(cols_ok) == 1) {
                gene_data <- setNames(sub[, cols_ok, drop = TRUE], canon_ids(rownames(sub)))
              } else if (length(cols_ok) >= 2) {
                gd <- sub[, cols_ok, drop = FALSE]
                rownames(gd) <- canon_ids(rownames(sub))
                gene_data <- gd
              }
            }
          }
          
          if (length(gene_data) == 0) {
            showNotification("\u26a0\ufe0f No genes matched the pathway.",
                             type = "warning", duration = 6)
            return()
          }
          
          # ── Back to pathview's ID space ──────────────────────────────────
          # gene_data is keyed on Solyc at this point (either from the matrix
          # rownames or the rep(1) fallback), but gene.idtype = "KEGG" means
          # pathview matches on Entrez. Translate the names back and drop any
          # gene we cannot translate — pathview would silently ignore those
          # anyway, and dropping them keeps the colour scale honest by not
          # letting untranslatable entries sit in the data unused.
          if (!is.null(gene_to_kegg)) {
            nm   <- if (is.matrix(gene_data)) rownames(gene_data) else names(gene_data)
            keep <- nm %in% names(gene_to_kegg)
            if (!any(keep)) {
              showNotification(
                "\u26a0\ufe0f None of the pathway's genes could be mapped back to KEGG IDs.",
                type = "warning", duration = 6)
              return()
            }
            new_nm <- unname(gene_to_kegg[nm[keep]])
            if (is.matrix(gene_data)) {
              gene_data <- gene_data[keep, , drop = FALSE]
              rownames(gene_data) <- new_nm
              gene_data <- gene_data[!duplicated(rownames(gene_data)), , drop = FALSE]
            } else {
              gene_data <- gene_data[keep]
              names(gene_data) <- new_nm
              gene_data <- gene_data[!duplicated(names(gene_data))]
            }
          }
          
          gene_idtype <- if (identical(ACTIVE_SPECIES, "arabidopsis")) "TAIR" else "KEGG"
          pv_suffix     <- paste0("cheski_", format(Sys.time(), "%H%M%S"))
          if (is.data.frame(gene_data)) gene_data <- as.matrix(gene_data)
          storage.mode(gene_data) <- "numeric"
          is_two_groups <- is.matrix(gene_data) && ncol(gene_data) == 2
          
          # Run pathview silently — PNG may be written even if an internal error occurs
          try(silent = TRUE, {
            if (is_two_groups) {
              suppressWarnings(
                pathview::pathview(
                  gene.data   = gene_data,
                  kegg.dir = tmp_dir,
                  pathway.id  = path_num,
                  species     = species,
                  gene.idtype = gene_idtype,
                  out.suffix  = pv_suffix,
                  kegg.native = TRUE,
                  same.layer  = FALSE,
                  limit       = list(gene = 2, cpd = 1)
                )
              )
            } else {
              gene_range <- max(abs(gene_data), na.rm = TRUE)
              gene_range <- if (is.finite(gene_range) && gene_range > 0) gene_range else 1
              pathview::pathview(
                gene.data   = gene_data,
                pathway.id  = path_num,
                kegg.dir = tmp_dir,
                species     = species,
                gene.idtype = gene_idtype,
                out.suffix  = pv_suffix,
                kegg.native = TRUE,
                low         = list(gene = "#00008B", cpd = "blue"),
                mid         = list(gene = "#f5f5f5", cpd = "gray"),
                high        = list(gene = "#FFD700", cpd = "yellow"),
                limit       = list(gene = gene_range, cpd = 1)
              )
            }
          })
          
          try(silent = TRUE, while (!is.null(dev.list())) dev.off())
          
          # Resolve PNG path (multi for two groups, standard otherwise)
          png_file <- file.path(tmp_dir, paste0(raw_id, ".", pv_suffix, ".png"))
          if (is_two_groups) {
            png_multi <- file.path(tmp_dir, paste0(raw_id, ".", pv_suffix, ".multi.png"))
            if (file.exists(png_multi)) png_file <- png_multi
          }
          
          message(">>> PNG exists: ", file.exists(png_file))
          show_pathview_modal(png_file)
          
        }, error = function(e) {
          showNotification(paste("\u274c Pathview error:", e$message),
                           type = "error", duration = 8)
        })
      })
    })
    
    # Download handler for pathview PNG
    output$dl_pathview <- downloadHandler(
      filename = function() {
        sel <- input$table_rows_selected
        pid <- if (!is.null(sel))
          gsub("<[^>]+>", "", result()$table[sel, "ID"]) else "pathway"
        paste0(pid, "_pathview.png")
      },
      content = function(f) {
        req(pathview_png_path())
        file.copy(pathview_png_path(), f)
      }
    )
    
    # Update pin_terms choices when enrichment result is available
    observe({
      req(result())
      tbl <- result()$table
      req(!is.null(tbl), "Description" %in% colnames(tbl))
      updateSelectizeInput(session, "pin_terms",
                           choices  = tbl$Description,
                           selected = character(0),
                           server   = TRUE)
    })
    
    plot_obj <- reactive({
      req(result())
      tbl <- result()$table
      if (is.null(tbl) || !"Description" %in% colnames(tbl)) {
        return(
          ggplot() +
            annotate("text", x = 0.5, y = 0.5,
                     label = "No significant KEGG pathways found.",
                     size = 6, color = "gray40") +
            theme_void()
        )
      }
      # make_enrichment_plot() needs the {Description, GeneRatio, p.adjust,
      # Count} columns run_kegg_enrichment() already produces. enrichKEGG's
      # result additionally carries FoldEnrichment natively, so the "Enrichment
      # Factor" x-axis is a real column here and not a silent fallback to
      # GeneRatio.
      do.call(make_enrichment_plot, c(
        list(
          result_df = tbl,
          title     = "KEGG Pathway Enrichment",
          top_n     = input$top_n     %||% 20,
          pinned    = input$pin_terms %||% character(0),
          font_size = input$plot_fs   %||% 13,
          sort_desc = if (is.null(input$sort_desc)) TRUE else isTRUE(input$sort_desc)
        ),
        enrich_plot_style_args(input)
      ))
    })
    
    output$plot_container <- renderUI({
      h <- if (is.null(input$plot_h) || !is.numeric(input$plot_h)) 800 else input$plot_h
      plotOutput(session$ns("plot"), height = paste0(h, "px"))
    })
    output$plot <- renderPlot({ plot_obj() })
    
    # Downloads
    output$dl_table <- downloadHandler(
      filename = function() paste0("KEGG_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
      # Carries the ID-mapping loss note, which for tomato is usually the
      # single most important thing to read alongside this table.
      content  = function(file) { req(result()); write_result_csv(result(), file) }
    )
    # Same shared export sizing as the GO Plot tab: both formats read the two
    # inch fields, and the PNG's pixel count and dpi move together so the canvas
    # is exactly the size asked for.
    .exp <- enrich_export_server(input, output, session,
                                 result_r = result, plot_r = plot_obj)
    output$dl_pdf <- downloadHandler(
      filename = function() paste0("KEGG_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf"),
      content  = .exp$pdf
    )
    output$dl_png <- downloadHandler(
      filename = function() paste0("KEGG_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
      content  = .exp$png
    )

    # ── Accumulate EVERY run this session (for full-history reproducibility) ──
    # keggServer returned nothing at all, so KEGG was the one enrichment path
    # with no per-run record: the bundle described it by re-reading live inputs
    # at download time, which reported the last widget state rather than what
    # each run used, and named two widgets that do not exist.
    #
    # NO ignoreInit here, deliberately -- same reasoning as go_module.r:620-626.
    # result() is an eventReactive that ABORTS at init, so ignoreInit's skip
    # would be consumed by the user's FIRST Run KEGG click instead of by init,
    # and run 1 would be permanently missing. The is.null(r$params) guard below
    # already makes an init-time fire a no-op.
    kegg_history <- reactiveVal(list())
    observeEvent(result(), {
      r <- result()
      if (is.null(r) || is.null(r$params)) return()
      p <- r$params
      # Same provenance stamp as go_module.R. keggServer() has no `context`
      # argument, so the module's own `id` is used instead -- "kegg" for the
      # heatmap panel, "kegg_de" for the DEG panel. Taking it from `id` rather
      # than a new argument means the label cannot drift from the module it
      # actually belongs to. run_at is the sort key app.r merges the two
      # panels' histories with.
      p$context <- id
      p$run_at  <- format(Sys.time(), "%Y-%m-%d %H:%M:%OS3")
      kegg_history(c(kegg_history(), list(p)))
    })

    return(list(result = result, history = kegg_history))
  })
}