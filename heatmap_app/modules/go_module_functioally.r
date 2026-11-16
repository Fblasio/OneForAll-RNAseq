goControlsUI <- function(id) {
  ns <- NS(id)
  tagList(
    h5("GO Settings"),
    selectInput(ns("ontology"), "GO Ontology",
                choices = c("BP","MF","CC","ALL"), selected = "BP"),
    selectInput(ns("padj_method"), "p.adjust method",
                choices = c("BH","BY","bonferroni","holm","hochberg","hommel","none"),
                selected = "BH"),
    numericInput(ns("pval"), "p-value cutoff", value = 0.05, min = 0, max = 1, step = 0.01),
    numericInput(ns("qval"), "q-value cutoff", value = 0.05, min = 0, max = 1, step = 0.01),
    numericInput(ns("simplify_cutoff"), "Simplify cutoff",
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
    # Slider to control how many GO terms appear in the dotplot
    fluidRow(
      column(6, sliderInput(ns("top_n"), "Top N GO terms", 5, 100, 20, 5, ticks=FALSE))
    ),
    # Selectize to pin specific GO terms — they are always shown regardless of rank
    fluidRow(
      column(12, selectizeInput(ns("pin_terms"), "Always show these GO terms (pinned):",
        choices = NULL, multiple = TRUE,
        options = list(placeholder = "Select GO terms to pin...",
                       plugins = list("remove_button"))))
    ),
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

# ── helper: pure ggplot2 dotplot, no enrichplot ─────────────
# pinned: character vector of Description values that must always appear in the plot.
# The remaining slots (top_n - length(pinned)) are filled with the top free terms by GeneRatio.
.make_dotplot <- function(tbl, title = "", top_n = 20, pinned = character(0), font_size = 13) {
  df <- tbl
  df$Description <- as.character(df$Description)
  df$GeneRatio   <- as.character(df$GeneRatio)
  df$p.adjust    <- as.numeric(as.character(df$p.adjust))
  df$Count       <- as.numeric(as.character(df$Count))

  df$GeneRatioNum <- sapply(df$GeneRatio, function(x) {
    p <- strsplit(x, "/")[[1]]
    if (length(p) == 2) as.numeric(p[1]) / as.numeric(p[2]) else NA_real_
  })

  df <- df[!is.na(df$GeneRatioNum) & !is.na(df$p.adjust), ]
  df <- df[order(df$GeneRatioNum, decreasing = TRUE), ]
  # Apply pinned logic: pinned terms are always shown;
  # remaining free slots (top_n - n_pinned) are filled with top-ranked free terms.
  if (length(pinned) > 0) {
    is_pinned  <- df$Description %in% pinned
    df_pinned  <- df[is_pinned,  ]
    df_free    <- df[!is_pinned, ]
    n_free     <- max(0L, top_n - nrow(df_pinned))
    df         <- rbind(df_free[seq_len(min(n_free, nrow(df_free))), ], df_pinned)
    df         <- df[order(df$GeneRatioNum, decreasing = TRUE), ]
  } else {
    df <- head(df, top_n)
  }
  df$Description <- factor(df$Description, levels = rev(df$Description))

  ggplot(df, aes(x = GeneRatioNum, y = Description,
                 size = Count, color = p.adjust)) +
    geom_point() +
    scale_color_gradient(low = "#C0392B", high = "#3B7DD8", name = "p.adjust") +
    scale_size_continuous(name = "Gene count", range = c(3, 10)) +
    xlab("Gene Ratio") + ylab(NULL) +
    ggtitle(title) +
    theme_bw(base_size = font_size) +
    theme(axis.text.y  = element_text(size = font_size),
          axis.text.x  = element_text(size = font_size),
          axis.title.x = element_text(size = font_size + 2, face = "bold"),
          plot.title   = element_text(size = font_size + 4, face = "bold"),
          legend.text  = element_text(size = font_size - 1),
          legend.title = element_text(size = font_size, face = "bold"))
}

goServer <- function(id, selected_genes,
                     scaled_matrix_rv = NULL,
                     gene_labels_rv   = NULL,
                     universe         = NULL,
                     on_new_heatmap   = NULL,
                     col_fun_rv       = NULL,
                     col_order_rv     = NULL,
                     deg_genes_rv     = NULL,
                     context = "default") {
  moduleServer(id, function(input, output, session) {

    gene_cart <- reactiveVal(character(0))

    # ── enrichment ─────────────────────────────────────────
    result <- eventReactive(input$run, {
      req(selected_genes())
      uni <- if (is.reactive(universe) || is.function(universe)) {
        tryCatch(universe(), error = function(e) NULL)
      } else {
        universe
      }
      
      withProgress(message = "Running GO enrichment...", {
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
            message("GO ERROR: ", e$message)
            list(table = data.frame(Error = e$message), object = NULL)
          }
        )
      })
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
      message(">>> USANDO NUEVO plot_obj <<<")
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
      .make_dotplot(tbl,
               title     = paste("GO", input$ontology, "Enrichment"),
               top_n     = if (is.null(input$top_n))     20            else input$top_n,
               pinned    = if (is.null(input$pin_terms)) character(0)  else input$pin_terms,
               font_size = if (is.null(input$plot_fs))   13            else input$plot_fs)
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
      
      Heatmap(
        sub_mat,
        name = "zscore",
        cluster_columns = FALSE,
        show_row_names = TRUE,
        row_labels = row_labels,
        row_names_gp = gpar(fontsize = fs, col = go_col, fontface = go_face),
        column_names_gp = gpar(fontsize = fs, fontface = "bold"),
        col = col_fun_go,
        heatmap_legend_param = list(
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
      content  = function(f) { req(result()); write.csv(result()$table, f, row.names = FALSE) }
    )
    output$dl_pdf <- downloadHandler(
      filename = function() paste0("GO_", input$ontology, "_", ts(), ".pdf"),
      content  = function(f) { pdf(f, width=10, height=8); print(plot_obj()); dev.off() }
    )
    output$dl_png <- downloadHandler(
      filename = function() paste0("GO_", input$ontology, "_", ts(), ".png"),
      content  = function(f) { png(f, width=1200, height=900, res=150); print(plot_obj()); dev.off() }
    )
    output$dl_ht_pdf <- downloadHandler(
      filename = function() paste0("GO_heatmap_", ts(), ".pdf"),
      content  = function(f) { ht <- go_heatmap_obj(); req(!is.null(ht)); pdf(f, width=8, height=10); draw(ht); dev.off() }
    )
    output$dl_ht_png <- downloadHandler(
      filename = function() paste0("GO_heatmap_", ts(), ".png"),
      content  = function(f) { ht <- go_heatmap_obj(); req(!is.null(ht)); png(f, width=900, height=1200, res=150); draw(ht); dev.off() }
    )
    
    # Expose the last-run reactive (table + object + the exact gene_ids/
    # universe/cutoffs used) so the reproducibility bundle can freeze the
    # real run instead of re-reading possibly-changed live inputs.
    return(list(result = result))
  })
}
