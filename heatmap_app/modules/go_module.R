goControlsUI <- function(id) {
  ns <- NS(id)
  tagList(
    h5("GO Settings"),
    selectInput(ns("ontology"),
                lbl("GO Ontology", "BP = biological process (the usual choice), MF = molecular function, CC = cellular component. ALL pools the three into ONE p-value correction and skips the redundancy filter, so expect overlapping terms."),
                choices = c("BP","MF","CC","ALL"), selected = "BP"),
    selectInput(ns("padj_method"),
                lbl("p.adjust method", "BH controls the false discovery rate - the default, and the right choice for enrichment. Bonferroni / holm / hochberg / hommel are far stricter. 'none' means NO correction at all; never report those p-values as significant."),
                choices = c("BH","BY","bonferroni","holm","hochberg","hommel","none"),
                selected = "BH"),
    numericInput(ns("pval"), "p-value cutoff", value = 0.05, min = 0, max = 1, step = 0.01),
    numericInput(ns("qval"), "q-value cutoff", value = 0.05, min = 0, max = 1, step = 0.01),
    numericInput(ns("simplify_cutoff"),
                 lbl("Simplify cutoff", "Collapses redundant parent/child GO terms by semantic similarity. Lower = more aggressive; 1 = off. The term kept is the one with the lowest p.adjust, which favours the more general of each pair. Ignored when Ontology = ALL."),
                 value = 0.7, min = 0, max = 1, step = 0.05),
    checkboxInput(ns("use_cache"), "Use cache", value = TRUE),
    actionButton(ns("run"), "\u25b6 Run GO", class = "btn-info btn-sm", width = "100%")
  )
}

goTableUI <- function(id) {
  ns <- NS(id)
  tagList(
    br(),
    downloadButton(ns("dl_table"), "\u2b07 Download CSV", class = "btn-sm btn-success"),
    br(), br(),
    helpText("\U0001f4a1 Click a row \u2192 see those genes in GO Heatmap tab."),
    DTOutput(ns("table"))
  )
}

goPlotUI <- function(id) {
  ns <- NS(id)
  tagList(
    br(),
    fluidRow(
      column(6, sliderInput(ns("plot_h"),  "Plot height (px)", 300, 1600, 800, 50, ticks=FALSE)),
      column(6, sliderInput(ns("plot_fs"), "Font size",        6,   28,  13,  1,  ticks=FALSE))
    ),
    # Slider to control how many GO terms appear in the plot
    fluidRow(
      column(6, sliderInput(ns("top_n"), "Top N GO terms", 5, 100, 20, 5, ticks=FALSE)),
      column(6, checkboxInput(ns("sort_desc"), "Sort highest first", value = TRUE))
    ),
    # Selectize to pin specific GO terms — they are always shown regardless of rank
    fluidRow(
      column(12, selectizeInput(ns("pin_terms"), "Always show these GO terms (pinned):",
        choices = NULL, multiple = TRUE,
        options = list(placeholder = "Select GO terms to pin...",
                       plugins = list("remove_button"))))
    ),
    # ── Plot style / customization panel ──────────────────────────────────
    # Defined in enrichPlotStyleUI() (utils/enrichment_functions.R) and shared
    # with the KEGG Plot tab, so the two panels offer the same controls with the
    # same defaults. The input IDs are the argument names of
    # make_enrichment_plot(), which is what plot_obj() below relies on.
    # Folded by default: these are set-once controls and this tab is already
    # tall.
    enrichPlotStyleUI(ns, terms_label = "GO terms"),
    # Export canvas size. The downloads used to be a fixed 10x8in PDF and an
    # 8x6in PNG that ignored every size control on this tab; see
    # enrichPlotExportUI() in utils/enrichment_functions.R.
    enrichPlotExportUI(ns),
    fluidRow(
      column(4, downloadButton(ns("dl_pdf"), "\u2b07 PDF", class = "btn-sm btn-success")),
      column(4, downloadButton(ns("dl_png"), "\u2b07 PNG", class = "btn-sm btn-info"))
    ),
    br(),
    uiOutput(ns("plot_container"))
  )
}

goHeatmapUI <- function(id) {
  ns <- NS(id)
  tagList(
    br(),
    helpText("Select a row in the GO Table tab first."),
    uiOutput(ns("term_info")),
    br(),
    fluidRow(
      column(6, sliderInput(ns("ht_h"),  "Heatmap height (px)", 300, 1600, 600, 50, ticks=FALSE)),
      column(6, sliderInput(ns("ht_fs"), "Font size",           6,   28,  11,  1,  ticks=FALSE))
    ),
    # Its own export block: this tab lives in the same module namespace as the
    # dot plot, so the ids are prefixed "ht_exp" to avoid colliding with it.
    # Sized from GENES rather than terms -- one row per gene in the GO term.
    enrichPlotExportUI(ns, prefix = "ht_exp", unit = "genes",
                       default_w = 8, default_h = 10),
    fluidRow(
      column(4, downloadButton(ns("dl_ht_pdf"), "\u2b07 PDF", class = "btn-sm btn-success")),
      column(4, downloadButton(ns("dl_ht_png"), "\u2b07 PNG", class = "btn-sm btn-info"))
    ),
    br(),
    uiOutput(ns("go_heatmap_container")),
    hr(),
    h5("\U0001f9ec Gene Picker", style = "color:#1F3864; font-weight:bold;"),
    helpText("Pick genes from this heatmap and accumulate them into a custom list."),
    fluidRow(
      column(8, selectizeInput(ns("picker_genes"), label = NULL, choices = NULL,
               multiple = TRUE, options = list(placeholder = "Genes in current GO heatmap...",
                                               plugins = list("remove_button")))),
      column(4,
        actionButton(ns("add_to_cart"),     "\u2795 Add selected",
                     class = "btn-info btn-sm",    width = "100%"), br(), br(),
        actionButton(ns("add_all_to_cart"), "\u2795 Add ALL",
                     class = "btn-default btn-sm", width = "100%"))
    ),
    hr(),
    h5("\U0001f6d2 Custom gene list", style = "color:#1F3864; font-weight:bold;"),
    uiOutput(ns("cart_summary")),         # 1. conteo + lista
    br(),
    uiOutput(ns("cart_remove_ui")),       # 2. selector para quitar genes
    br(),
    fluidRow(
      column(6, actionButton(ns("new_heatmap"), "\U0001f525 New Heatmap from list",
                             class = "btn-primary btn-sm", width = "100%")),
      column(6, actionButton(ns("clear_cart"), "\U0001f5d1\ufe0f Clear list",
                             class = "btn-danger btn-sm", width = "100%"))
    )
  )
}

# ── helper: GO dotplot ───────────────────────────────────────────────────
# .make_dotplot() now just forwards to the shared, customizable
# make_enrichment_plot() (see enrichment_functions.R) with style = "dot" and
# the historical default colors, so any call site that predates the new
# plot-customization panel keeps producing pixel-equivalent output. New code
# (see plot_obj() below) calls make_enrichment_plot() directly to access the
# full set of style/x-axis/color options.
.make_dotplot <- function(tbl, title = "", top_n = 20, pinned = character(0), font_size = 13) {
  make_enrichment_plot(
    tbl, title = title, top_n = top_n, pinned = pinned, font_size = font_size,
    style = "dot", x_metric = "GeneRatio",
    color_mode = "gradient", color_metric = "p.adjust",
    color_low = "#C0392B", color_high = "#3B7DD8",
    size_by_count = TRUE, label_wrap = 0, theme_style = "bw"
  )
}

goServer <- function(id, selected_genes,
                     scaled_matrix_rv = NULL,
                     # Both optional and NULL-defaulted, so existing callers
                     # that don't pass them (e.g. "go_de", which has no VST
                     # twin) keep working unchanged.
                     cluster_matrix_rv = NULL,
                     linkage_rv       = NULL,
                     gene_labels_rv   = NULL,
                     universe         = NULL,
                     on_new_heatmap   = NULL,
                     col_fun_rv       = NULL,
                     col_order_rv     = NULL,
                     deg_genes_rv     = NULL,
                     scale_label_rv   = NULL,
                     context = "default") {
  moduleServer(id, function(input, output, session) {

    gene_cart <- reactiveVal(character(0))

    # ── enrichment ─────────────────────────────────────────
    result <- eventReactive(input$run, {
      genes <- selected_genes()
      if (is.null(genes) || length(genes) == 0) {
        showNotification(
          "\u26a0\ufe0f No genes to run GO on yet. Load/select a gene list first, then click Run GO.",
          type = "warning", duration = 8)
      }
      req(genes, cancelOutput = TRUE)
      uni <- if (is.reactive(universe) || is.function(universe)) {
        tryCatch(universe(), error = function(e) NULL)
      } else {
        universe
      }
      
      # attach_plot_snapshot(): freezes the Plot-style panel as it is RIGHT NOW,
      # so the reproducibility bundle can regenerate this run's figure instead
      # of describing whatever the (live) style controls say at download time.
      # Applied to both the success and the failure branch, so a failed run is
      # still recorded with the settings it was attempted under.
      # with_progress_messages() rather than withProgress(): same panel, but it
      # also relays enrichGO()'s and simplify()'s own messages ("preparing gene
      # to GO mapping data...") as they arrive, so a slow GO term/gene mapping
      # step no longer looks like a hang. See log_utils.R.
      attach_plot_snapshot(with_progress_messages({
        tryCatch(
          run_go_enrichment(
            selected_genes(),
            process         = input$ontology,
            padj_method     = input$padj_method,
            pval_cutoff     = input$pval,
            qval_cutoff     = input$qval,
            simplify_cutoff = input$simplify_cutoff,
            use_cache       = input$use_cache,
            universe        = uni,
            context         = context
          ),
          error = function(e) {
            # A params block is attached even on failure, because the history
            # observer below skips anything without one -- so a crashed run used
            # to vanish entirely and a bundle showing 3 GO runs could come from
            # a session that attempted 6. `failed` lets rerun_GO.R note the
            # attempt without trying to reproduce it.
            log_error(paste0("GO run failed: ", conditionMessage(e)))
            list(table  = data.frame(Error = conditionMessage(e)),
                 object = NULL,
                 params = list(context       = context,
                               species       = ACTIVE_SPECIES,
                               process       = input$ontology,
                               padj_method   = input$padj_method,
                               pval_cutoff   = input$pval,
                               qval_cutoff   = input$qval,
                               gene_ids      = tryCatch(canon_ids(selected_genes()),
                                                        error = function(e2) character(0)),
                               failed        = TRUE,
                               error         = conditionMessage(e)))
          }
        )
      }, "Running GO enrichment"), input, paste("GO", input$ontology, "Enrichment"))
    })

    observeEvent(result(), {
      for (w in result()$warnings %||% character(0)) {
        showNotification(paste0("\u26a0\ufe0f ", w), type = "warning", duration = 14)
        log_warning(w)
      }
    })

    # ── table ──────────────────────────────────────────────
    output$table <- renderDT({
      req(result())
      datatable(result()$table,
                selection = "single",
                options   = list(pageLength = 10, scrollX = TRUE),
                rownames  = FALSE)
    })

    # ── update pin_terms choices whenever enrichment result changes ──────
    observe({
      req(result())
      tbl <- result()$table
      req(!is.null(tbl), "Description" %in% colnames(tbl))
      updateSelectizeInput(session, "pin_terms",
                           choices  = tbl$Description,
                           selected = character(0),
                           server   = TRUE)
    })

    # ── plot ───────────────────────────────────────────────
    plot_obj <- reactive({
      req(result())
      tbl <- result()$table

      if (is.null(tbl) || !"Description" %in% colnames(tbl)) {
        return(
          ggplot() +
            annotate("text", x = 0.5, y = 0.5,
                     label = "No significant GO terms found.",
                     size = 6, color = "gray40") +
            theme_void()
        )
      }
      # The ten style arguments come from enrich_plot_style_args(), the reader
      # that pairs with enrichPlotStyleUI(). The KEGG panel composes its call
      # the same way, so neither can drift from the other.
      do.call(make_enrichment_plot, c(
        list(
          result_df = tbl,
          title     = paste("GO", input$ontology, "Enrichment"),
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

    # ── selected row genes ─────────────────────────────────
    go_genes_selected <- reactive({
      req(result())
      sel <- input$table_rows_selected
      req(sel)
      tbl <- result()$table
      req(nrow(tbl) >= sel)
      gene_col <- tbl[sel, "geneID"]
      req(!is.null(gene_col), !is.na(gene_col), nzchar(as.character(gene_col)))
      genes <- trimws(strsplit(as.character(gene_col), "/")[[1]])
      genes <- canon_ids(genes)
      genes
    })

    output$term_info <- renderUI({
      req(result())
      sel <- input$table_rows_selected
      req(sel)
      tbl <- result()$table
      req(nrow(tbl) >= sel)
      genes <- tryCatch(go_genes_selected(), error = function(e) NULL)
      n <- if (!is.null(genes)) length(genes) else "?"
      tags$p(tags$b(paste0("GO term: ", tbl[sel, "Description"], " (", n, " genes)")))
    })

    # ── heatmap ────────────────────────────────────────────
    go_heatmap_obj <- reactive({
      mat   <- scaled_matrix_rv()
      gl    <- gene_labels_rv()
      genes <- go_genes_selected()
      
      req(!is.null(mat), length(genes) > 0)
      
      # rownames(mat), selected_genes() and genes are already canonical for
      # tomato (matrix canonicalized at import; GO term genes canonicalized in
      # go_genes_selected()). Only trimming is required.
      sel_idx     <- which(ids_in(rownames(mat), selected_genes()))
      mat_sub     <- mat[sel_idx, , drop = FALSE]
      mat_sub_ids <- rownames(mat_sub)
      
      # message("=== GO heatmap DEBUG ===")
      # message("N genes from GO term: ", length(go_ids_for_match))
      # message("Head genes from GO term: ", paste(utils::head(go_ids_for_match, 10), collapse = ", "))
      # message("N rows in scaled matrix: ", nrow(mat))
      # message("Head matrix rownames: ", paste(utils::head(rownames(mat), 10), collapse = ", "))
      
      message("N rows after selected_genes() subset: ", nrow(mat_sub))
      message("Head rownames after selected_genes() subset: ", paste(utils::head(rownames(mat_sub), 10), collapse = ", "))
      
      if (nrow(mat_sub) == 0) {
        msg <- explain_empty("GO heatmap", length(selected_genes()), 0,
                             "None of the currently selected genes are present in the expression matrix.")
        showNotification(paste0("\u26a0\ufe0f ", msg), type = "warning", duration = 6)
        log_warning(msg)
        return(NULL)
      }
      
      go_idx <- ids_match(genes, mat_sub_ids)
      go_idx <- go_idx[!is.na(go_idx)]
      
      message("N GO genes matched in current matrix: ", length(go_idx))
      message("Head matched GO genes: ", paste(utils::head(mat_sub_ids[go_idx], 10), collapse = ", "))
      message("================================")
      
      if (length(go_idx) == 0) {
        msg <- explain_empty("GO heatmap", length(genes), 0,
                             "None of this GO term's genes are in your current sub-heatmap selection.")
        showNotification(paste0("\u26a0\ufe0f ", msg), type = "warning", duration = 6)
        log_warning(msg)
        return(NULL)
      }
      
      sub_mat <- mat_sub[go_idx, , drop = FALSE]
      row_labels <- mat_sub_ids[go_idx]
      
      degs_go <- if (!is.null(deg_genes_rv)) {
        tryCatch(canon_ids(deg_genes_rv()), error = function(e) character(0))
      } else character(0)
      
      is_deg_go <- row_labels %in% degs_go
      go_col <- ifelse(is_deg_go, "#FF00FF", "black")
      go_face <- ifelse(is_deg_go, "bold", "plain")
      
      if (!is.null(col_order_rv)) {
        ord <- tryCatch(col_order_rv(), error = function(e) NULL)
        if (!is.null(ord)) {
          ord_valid <- intersect(ord, colnames(sub_mat))
          if (length(ord_valid) == ncol(sub_mat)) {
            sub_mat <- sub_mat[, ord_valid, drop = FALSE]
          }
        }
      }
      
      if (!is.null(gl) && length(gl) > 0) {
        gl_names <- names(gl)
        gl_vals  <- as.character(gl)
        idx <- match(row_labels, gl_names)
        valid <- !is.na(idx)
        if (any(valid)) {
          row_labels[valid] <- gl_vals[idx[valid]]
        }
      }
      
      fs <- if (is.null(input$ht_fs)) 11 else input$ht_fs
      
      col_fun_go <- if (!is.null(col_fun_rv)) {
        tryCatch(
          col_fun_rv(),
          error = function(e) colorRamp2(c(-2, 0, 2), c("#3B5998", "white", "#C0392B"))
        )
      } else {
        colorRamp2(c(-2, 0, 2), c("#3B5998", "white", "#C0392B"))
      }
      
      # ── Row clustering ───────────────────────────────────────────────────
      # cluster_rows was never passed below, so this heatmap silently used
      # ComplexHeatmap's default: TRUE, euclidean, on sub_mat — the DISPLAY
      # matrix. That contradicts the main heatmap (Pearson, on the VST
      # substrate), so the same genes could come out in a different row order
      # here than in the heatmap they were brushed from, with no warning.
      #
      # Now: cluster on cluster_matrix_rv when it's available, through the
      # same cluster_rows_or_reason() contract, with the same linkage.
      # Column order is irrelevant to row-row correlation (permuting columns
      # permutes the paired observations identically for every row), so the
      # column mismatch between sub_mat and the substrate is harmless.
      # When no substrate was supplied ("go_de" has no VST twin), fall back to
      # sub_mat as before — but under the guards, and stated explicitly rather
      # than inherited from a library default.
      clust_src_go <- if (!is.null(cluster_matrix_rv)) {
        tryCatch(cluster_matrix_rv(), error = function(e) NULL)
      } else NULL
      sub_clust_go <- NULL
      if (!is.null(clust_src_go)) {
        idx_go <- match(rownames(sub_mat), rownames(clust_src_go))
        # The substrate is sample-level, so its columns do not match sub_mat's
        # group columns. Column identity is irrelevant to row-row correlation
        # and every substrate sample is a valid observation, so all are used.
        # Only row alignment is checked.
        if (!anyNA(idx_go)) {
          sub_clust_go <- clust_src_go[idx_go, , drop = FALSE]
        }
      }
      if (is.null(sub_clust_go)) sub_clust_go <- sub_mat
      
      linkage_go <- if (!is.null(linkage_rv)) {
        tryCatch(linkage_rv(), error = function(e) "complete")
      } else "complete"
      cluster_rows_go <- cluster_rows_or_reason(sub_clust_go,
                                                linkage = linkage_go)$hc
      
      Heatmap(
        sub_mat,
        name = "zscore",
        cluster_rows = cluster_rows_go,
        cluster_columns = FALSE,
        show_row_names = TRUE,
        row_labels = row_labels,
        row_names_gp = gpar(fontsize = fs, col = go_col, fontface = go_face),
        column_names_gp = gpar(fontsize = fs, fontface = "bold"),
        col = col_fun_go,
        heatmap_legend_param = list(
          title     = if (is.function(scale_label_rv)) scale_label_rv() else "zscore",
          title_gp = gpar(fontsize = fs, fontface = "bold"),
          labels_gp = gpar(fontsize = fs - 1)
        )
      )
    })

    output$go_heatmap_container <- renderUI({
      h <- if (is.null(input$ht_h) || !is.numeric(input$ht_h)) 600 else input$ht_h
      plotOutput(session$ns("go_heatmap"), height = paste0(h, "px"))
    })
    output$go_heatmap <- renderPlot({
      ht <- go_heatmap_obj()
      req(!is.null(ht))
      grid::grid.newpage()
      draw(ht)
    })

    # ── gene picker ────────────────────────────────────────
    observe({
      mat   <- scaled_matrix_rv()
      genes <- tryCatch(go_genes_selected(), error = function(e) NULL)
      req(!is.null(mat), !is.null(genes), length(genes) > 0)
      in_sel  <- ids_in(rownames(mat), selected_genes())
      sub     <- mat[in_sel, , drop = FALSE]
      present <- rownames(sub)[ids_in(rownames(sub), genes)]
      gl      <- gene_labels_rv()
      if (!is.null(gl) && length(gl) > 0) {
        idx   <- ids_match(present, names(gl))
        valid <- !is.na(idx)
        lbls  <- present; lbls[valid] <- as.character(gl[idx[valid]])
        choices <- setNames(present, lbls)
      } else { choices <- present }
      updateSelectizeInput(session, "picker_genes",
                           choices = choices, selected = character(0))
    })

    observeEvent(input$add_to_cart, {
      req(length(input$picker_genes) > 0)
      combined <- unique(c(gene_cart(), canon_ids(input$picker_genes)))
      gene_cart(combined)
      showNotification(paste0("\u2795 ", length(input$picker_genes),
        " added. Total: ", length(combined)), type = "message", duration = 3)
    })

    observeEvent(input$add_all_to_cart, {
      mat   <- scaled_matrix_rv()
      genes <- tryCatch(go_genes_selected(), error = function(e) NULL)
      req(!is.null(mat), !is.null(genes), length(genes) > 0)
      in_sel  <- ids_in(rownames(mat), selected_genes())
      sub     <- mat[in_sel, , drop = FALSE]
      present <- rownames(sub)[ids_in(rownames(sub), genes)]
      combined <- unique(c(gene_cart(), present))
      gene_cart(combined)
      showNotification(paste0("\u2795 All ", length(present),
        " added. Total: ", length(combined)), type = "message", duration = 3)
    })

    observeEvent(input$clear_cart, {
      gene_cart(character(0))
      showNotification("\U0001f5d1 List cleared.", type = "warning", duration = 3)
    })
    
    # Dynamic selectize showing current cart genes — for individual removal
    output$cart_remove_ui <- renderUI({
      cart <- gene_cart()
      if (length(cart) == 0) return(NULL)
      tagList(
        fluidRow(
          column(8,
                 selectizeInput(session$ns("remove_genes"),
                                label = "Remove genes from list:",
                                choices  = cart,
                                selected = NULL,
                                multiple = TRUE,
                                options  = list(
                                  placeholder = "Select genes to remove...",
                                  plugins     = list("remove_button")
                                ))
          ),
          column(4,
                 br(),
                 actionButton(session$ns("remove_from_cart"), "\u2796 Remove selected",
                              class = "btn-warning btn-sm", width = "100%")
          )
        )
      )
    })
    
    observeEvent(input$remove_from_cart, {
      req(length(input$remove_genes) > 0)
      remaining <- setdiff(gene_cart(), canon_ids(input$remove_genes))
      gene_cart(remaining)
      showNotification(
        paste0("\u2796 ", length(input$remove_genes), " removed. Total: ", length(remaining)),
        type = "warning", duration = 3
      )
    })

    output$cart_summary <- renderUI({
      cart <- gene_cart()
      if (length(cart) == 0)
        tags$p(style = "color:gray; font-style:italic;", "No genes in list yet.")
      else
        tagList(
          tags$p(style = "color:#1F3864;",
                 tags$b(paste0(length(cart), " gene(s):"))),
          tags$ul(
            style = "max-height: 180px; overflow-y: auto; padding-left: 18px;
                 font-size: 0.85em; color: #1F3864; border: 1px solid #dee2e6;
                 border-radius: 4px; padding: 6px 6px 6px 22px; background: #f8f9fa;",
            lapply(cart, function(g) tags$li(g))
          )
        )
      })
    
    observeEvent(input$new_heatmap, {
      cart <- gene_cart()
      message(">>> new_heatmap clicked, cart length: ", length(cart))
      req(length(cart) > 0)
      if (!is.null(on_new_heatmap)) {
        message(">>> calling on_new_heatmap")
        on_new_heatmap(list(genes = cart, ts = Sys.time()))
      } else {
        message(">>> on_new_heatmap is NULL!")
      }
    })

    # ── downloads ──────────────────────────────────────────
    ts <- function() format(Sys.time(), "%Y%m%d_%H%M%S")

    output$dl_table <- downloadHandler(
      filename = function() paste0("GO_", input$ontology, "_", ts(), ".csv"),
      # write_result_csv() prepends this run's caveats as '#' lines -- most
      # importantly that simplify() adjusted p over MORE terms than are listed.
      content  = function(f) { req(result()); write_result_csv(result(), f) }
    )
    # Export size is driven by enrichPlotExportUI()'s two inch fields, which the
    # shared server keeps in step with Top N / font size / label wrap. Both
    # formats now read the SAME two numbers, so the PDF and the PNG can no
    # longer come out at different scales.
    .exp <- enrich_export_server(input, output, session,
                                 result_r = result, plot_r = plot_obj)
    output$dl_pdf <- downloadHandler(
      filename = function() paste0("GO_", input$ontology, "_", ts(), ".pdf"),
      content  = .exp$pdf
    )
    output$dl_png <- downloadHandler(
      filename = function() paste0("GO_", input$ontology, "_", ts(), ".png"),
      content  = .exp$png
    )
    # ── GO heatmap export ─────────────────────────────────────────────────────
    # Same three fixes as the dot plot above, but the row count is the number of
    # GENES in the selected GO term that are present in the current matrix, and
    # each row's label is a single line, so n_lines = 1.
    #
    # overhead_in is 1.2 rather than the dot plot's 0.7: a ComplexHeatmap adds
    # column names and padding below the body, and the heatmap title above it.
    # The colour legend sits on the RIGHT by default, so it costs width, not
    # height.
    ht_n_genes <- reactive({
      ht <- tryCatch(go_heatmap_obj(), error = function(e) NULL)
      if (is.null(ht)) return(NA_integer_)
      nrow(ht@matrix)
    })
    ht_needed_h <- reactive({
      n <- ht_n_genes()
      if (!isTRUE(is.finite(n))) return(NA_real_)
      enrich_plot_height_in(n, input$ht_fs %||% 11, n_lines = 1, overhead_in = 1.2)
    })
    .exp_ht <- enrich_export_server(
      input, output, session,
      prefix     = "ht_exp",
      unit       = "genes",
      needed_h_r = ht_needed_h,
      n_items_r  = ht_n_genes,
      font_r     = reactive(input$ht_fs %||% 11),
      draw_fun   = function() {
        ht <- go_heatmap_obj()
        req(!is.null(ht))
        ComplexHeatmap::draw(ht)
      })
    output$dl_ht_pdf <- downloadHandler(
      filename = function() paste0("GO_heatmap_", ts(), ".pdf"),
      content  = .exp_ht$pdf
    )
    output$dl_ht_png <- downloadHandler(
      filename = function() paste0("GO_heatmap_", ts(), ".png"),
      content  = .exp_ht$png
    )
    
    # ── Accumulate EVERY run this session (for full-history reproducibility) ──
    # Mirrors de_comparison_log.R: each successful GO run is appended so the
    # reproducibility bundle can regenerate all runs, not just the current one.
    # `result` (last run) is kept unchanged so existing consumers don't break.
    go_history <- reactiveVal(list())
    # NO ignoreInit here, deliberately. Unlike de_comparison_log.R:51 -- which
    # observes reactive(input$run_de) and whose init pass returns a value --
    # result() is an eventReactive that ABORTS at init (input$run == 0 is a null
    # event), so ignoreInit's skip was consumed by the user's FIRST Run GO click
    # instead of by init: run 1 was never recorded and the repro bundle was
    # permanently one GO run short. The is.null(r$params) guard below already
    # makes an init-time fire a no-op.
    observeEvent(result(), {
      r <- result()
      if (is.null(r) || is.null(r$params)) return()
      p <- r$params
      # ── Provenance, stamped here rather than inside run_go_enrichment() ─────
      # `context` is already an argument of this module ("heatmap" / "rnaseq" /
      # "go_section") and the FAILURE branch above already records it -- only the
      # success path dropped it. Without it, the bundle concatenates runs from
      # two panels with nothing saying which is which, and write_rerun_GO()'s
      # `label <- p$context %||% sprintf("GO_run_%d", idx)` always takes its
      # fallback. run_at is the sort key that lets app.r put the merged history
      # back into the order the runs actually happened in; sub-second precision
      # because two clicks can land in the same second.
      p$context <- context
      p$run_at  <- format(Sys.time(), "%Y-%m-%d %H:%M:%OS3")
      go_history(c(go_history(), list(p)))
    })

    # Expose the last-run reactive (table + object + the exact gene_ids/
    # universe/cutoffs used) so the reproducibility bundle can freeze the
    # real run instead of re-reading possibly-changed live inputs.
    # `history` is the accumulating list of every run's params this session.
    return(list(result = result, history = go_history))
  })
}
