# =============================================================================
# gene_detail_module.R
# Módulo Shiny: seleccionar un gen específico y ver su heatmap + z-scores
# Recibe:  selected_genes   (reactiveVal — genes del brush actual)
#          scaled_matrix_rv (reactiveVal — matriz escalada completa)
#          gene_labels_rv   (reactiveVal — named vector ID->symbol)
# =============================================================================

# ── UI ───────────────────────────────────────────────────────────────────────
geneDetailUI <- function(id) {
  ns <- NS(id)
  tagList(
    br(),
    fluidRow(
      column(5,
        selectizeInput(ns("gene_select"),
                       label    = "🔍 Selecciona o busca un gen:",
                       choices  = NULL,
                       selected = NULL,
                       options  = list(
                         placeholder = "Escribe un ID o símbolo...",
                         maxOptions  = 500
                       ))
      ),
      column(7,
        uiOutput(ns("gene_badge"))
      )
    ),
    hr(),
    fluidRow(
      column(7,
        h5("🎨 Z-score Heatmap"),
        fluidRow(
          column(6, sliderInput(ns("ht_h"),  "Height (px)", 150, 800, 200, 50, ticks=FALSE)),
          column(6, sliderInput(ns("ht_fs"), "Font size",   6,   24,  11,  1,  ticks=FALSE))
        ),
        uiOutput(ns("heatmap_container"))
      ),
      column(5,
        h5("📋 Z-score por grupo"),
        DTOutput(ns("table"))
      )
    )
  )
}

# ── Server ───────────────────────────────────────────────────────────────────
geneDetailServer <- function(id, selected_genes, scaled_matrix_rv, gene_labels_rv,
                             col_fun_rv   = NULL,
                             col_order_rv = NULL,
                             deg_genes_rv = NULL) {
  moduleServer(id, function(input, output, session) {

    # ── Actualizar choices cuando cambia la selección del brush ──────────
    observe({
      genes <- selected_genes()
      gl    <- gene_labels_rv()

      if (is.null(genes) || length(genes) == 0) {
        updateSelectizeInput(session, "gene_select",
                             choices = character(0), server = TRUE)
        return()
      }

      # Etiquetas del selector: "Symbol (AT_ID)" o solo "AT_ID"
      if (!is.null(gl) && length(gl) > 0) {
        idx <- ids_match(genes, names(gl))
        symbols <- ifelse(is.na(idx), genes, as.character(gl[idx]))
        labels  <- ifelse(symbols == genes,
                          genes,
                          paste0(symbols, "  (", genes, ")"))
      } else {
        labels <- genes
      }

      choices <- setNames(genes, labels)
      updateSelectizeInput(session, "gene_select",
                           choices = choices, server = TRUE)
    })

    # ── Fila de z-scores del gen seleccionado ─────────────────────────────
    gene_data <- reactive({
      req(input$gene_select, scaled_matrix_rv())
      mat  <- scaled_matrix_rv()
      gene <- input$gene_select
      idx  <- which(ids_in(rownames(mat), gene))
      if (length(idx) == 0) return(NULL)
      mat[idx[1], , drop = FALSE]
    })

    # ── Badge informativo del gen ─────────────────────────────────────────
    output$gene_badge <- renderUI({
      req(input$gene_select)
      gene <- input$gene_select
      gl   <- gene_labels_rv()
      idx <- ids_match(gene, names(gl))
      sym  <- if (!is.null(gl) && !is.na(idx)) gl[idx] else gene

      tair_url <- paste0(GENE_DETAIL_URL, gene)
      tags$p(
        style = "margin-top: 28px;",
        tags$strong(sym),
        if (sym != gene) tags$span(paste0(" · ", gene),
                                    style = "color:grey; font-size:0.9em;"),
        tags$a(paste0("🔗 ", GENE_DB_LABEL), href = tair_url, target = "_blank",
               style = "margin-left: 12px; font-size: 0.9em;")
      )
    })

    # ── Heatmap de 1 fila ─────────────────────────────────────────────────
    output$heatmap_container <- renderUI({
      h <- if (is.null(input$ht_h) || !is.numeric(input$ht_h)) 200 else input$ht_h
      plotOutput(session$ns("heatmap"), height = paste0(h, "px"))
    })
    output$heatmap <- renderPlot({
      req(gene_data())
      mat  <- gene_data()
      gl   <- gene_labels_rv()
      # Inherit column order from main heatmap
      if (!is.null(col_order_rv)) {
        ord <- tryCatch(col_order_rv(), error = function(e) NULL)
        if (!is.null(ord)) {
          ord_valid <- intersect(ord, colnames(mat))
          if (length(ord_valid) == ncol(mat))
            mat <- mat[, ord_valid, drop = FALSE]
        }
      }
      gene <- rownames(mat)
      idx <- ids_match(gene, names(gl))
      sym  <- if (!is.null(gl) && length(gl) > 0 && !is.na(idx)) gl[idx] else gene
      row_label <- if (sym != gene) paste0(sym, "\n(", gene, ")") else gene

      # Usar los mismos límites que el heatmap principal
      zmax     <- max(abs(mat), 2, na.rm = TRUE)
      
      col_fun <- if (!is.null(col_fun_rv)) {
        tryCatch(col_fun_rv(), error = function(e) {
          zmax <- max(abs(mat), 2, na.rm = TRUE)
          colorRamp2(c(-zmax, 0, zmax), c("#3B5998", "white", "#C0392B"))
        })
      } else {
        zmax <- max(abs(mat), 2, na.rm = TRUE)
        colorRamp2(c(-zmax, 0, zmax), c("#3B5998", "white", "#C0392B"))
      }
      
      degs_det  <- if (!is.null(deg_genes_rv)) tryCatch(deg_genes_rv(), error = function(e) character(0)) else character(0)
      is_deg_det <- ids_in(rownames(mat), degs_det)
      det_col    <- if (any(is_deg_det)) "#FF00FF" else "black"
      det_face   <- if (any(is_deg_det)) "bold"    else "plain"

      ht <- Heatmap(
        mat,
        name             = "zscore",
        col              = col_fun,
        cluster_rows     = FALSE,
        cluster_columns  = FALSE,
        show_row_names   = TRUE,
        row_names_side   = "left",
        row_labels       = row_label,
        show_column_names = TRUE,
        column_names_gp  = gpar(fontsize = if (is.null(input$ht_fs)) 12 else input$ht_fs + 1),
        row_names_gp = gpar(fontsize = if (is.null(input$ht_fs)) 11 else input$ht_fs,
                            fontface = det_face, col = det_col),
        rect_gp          = gpar(col = "black", lwd = 1.5),
        heatmap_legend_param = list(
          title_gp  = gpar(fontsize = if (is.null(input$ht_fs)) 11 else input$ht_fs,
                           fontface = "bold"),
          labels_gp = gpar(fontsize = if (is.null(input$ht_fs)) 10 else input$ht_fs - 1),
          grid_width  = unit(4, "mm"),
          grid_height = unit(4, "mm")
        ),
        # Mostrar el valor z-score dentro de cada celda
        cell_fun = function(j, i, x, y, width, height, fill) {
          val   <- mat[i, j]
          color <- if (abs(val) > 0.5) "#FF00FF" else "black"
          grid.text(sprintf("%.2f", val), x, y,
                    gp = gpar(fontsize = if (is.null(input$ht_fs)) 11 else input$ht_fs, fontface = "bold", col = color))
        }
      )
      draw(ht, padding = unit(c(8, 8, 8, 30), "mm"))
    })

    # ── Tabla de z-scores ─────────────────────────────────────────────────
    output$table <- renderDT({
      req(gene_data())
      mat <- gene_data()
      df  <- data.frame(
        Group    = colnames(mat),
        `Z-score` = round(as.numeric(mat[1, ]), 4),
        stringsAsFactors = FALSE,
        check.names      = FALSE
      )
      # Color de fondo proporcional al z-score
      datatable(
        df,
        rownames  = FALSE,
        selection = "none",
        options   = list(pageLength = 20, dom = "t", scrollX = TRUE)
      ) |>
        formatStyle(
          "Z-score",
          background = styleInterval(
            c(-2, -1, 0, 1, 2),
            c("#00008B", "#7B9BD5", "#f5f5f5", "#FFF176", "#FFD700", "#B8860B")
          ),
          color = styleInterval(
            c(-1.5, 1.5),
            c("white", "black", "black")
          )
        )
    })

  })
}
