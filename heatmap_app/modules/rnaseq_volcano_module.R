# =============================================================================
# rnaseq_volcano_module.R
# -----------------------------------------------------------------------------
# Shiny module for volcano plot visualization from differential expression
# results. The module receives a reactive DE results table and opens a modal
# dialog with interactive plot controls and download handlers.
# =============================================================================

rnaseqVolcanoUI <- function(id) {
  ns <- NS(id)
  
  tagList(
    uiOutput(ns("open_button_ui"))
  )
}

# ── Pure classification logic (extracted so it's testable in isolation) ──────
# Given padj/log2FoldChange vectors and the two cutoffs, returns the same
# 4-way label the app draws on the volcano plot and uses for point coloring.
classify_volcano <- function(padj, log2FoldChange, padj_cutoff = 0.05, lfc_cutoff = 0.5) {
  dplyr::case_when(
    padj < padj_cutoff & log2FoldChange >  lfc_cutoff ~ "Up",
    padj < padj_cutoff & log2FoldChange < -lfc_cutoff ~ "Down",
    padj < padj_cutoff                                ~ "Sig (low FC)",
    TRUE                                               ~ "NS"
  )
}

rnaseqVolcanoServer <- function(id, de_results_r) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    volcano_data <- reactive({
      req(de_results_r())
      df <- de_results_r()
      
      req(all(c("gene_id", "log2FoldChange", "padj") %in% colnames(df)))
      
      df <- df[!is.na(df$padj) & !is.na(df$log2FoldChange), , drop = FALSE]
      
      padj_cutoff <- input$volcano_padj_cutoff %||% 0.05
      lfc_cutoff  <- input$volcano_lfc_cutoff %||% 0.5
      
      df$significant_class <- classify_volcano(df$padj, df$log2FoldChange,
                                               padj_cutoff, lfc_cutoff)
      
      df$neglog10padj <- -log10(df$padj)
      
      df
    })
    
    volcano_top_labels <- reactive({
      req(volcano_data())
      df <- volcano_data()
      
      n_labels <- input$volcano_n_labels %||% 20
      n_labels <- max(0, as.integer(n_labels))
      
      if (n_labels == 0) {
        return(df[0, , drop = FALSE])
      }
      
      df |>
        dplyr::filter(significant_class %in% c("Up", "Down")) |>
        dplyr::arrange(padj) |>
        dplyr::slice_head(n = n_labels)
    })
    
    volcano_plot_r <- reactive({
      req(volcano_data())
      
      df <- volcano_data()
      top_genes <- volcano_top_labels()
      
      n_up <- sum(df$significant_class == "Up", na.rm = TRUE)
      n_down <- sum(df$significant_class == "Down", na.rm = TRUE)
      
      max_fc <- max(abs(df$log2FoldChange), na.rm = TRUE)
      if (!is.finite(max_fc) || max_fc <= 0) {
        max_fc <- 1
      }
      
      padj_cutoff <- input$volcano_padj_cutoff %||% 0.05
      lfc_cutoff  <- input$volcano_lfc_cutoff %||% 0.5
      
      point_size <- input$volcano_point_size %||% 2.2
      label_size <- input$volcano_label_size %||% 4.5
      axis_size  <- input$volcano_axis_size %||% 14
      title_size <- input$volcano_title_size %||% 18
      subtitle_size <- input$volcano_subtitle_size %||% 14
      
      col_up <- input$volcano_col_up %||% "#cc0000"
      col_down <- input$volcano_col_down %||% "#0055cc"
      col_sig_low <- input$volcano_col_sig_low %||% "#ff9900"
      col_ns <- input$volcano_col_ns %||% "grey75"
      
      plot_title <- input$volcano_title %||% "Volcano plot"
      
      ggplot2::ggplot(
        df[order(df$significant_class == "NS", decreasing = TRUE), , drop = FALSE],
        ggplot2::aes(x = log2FoldChange, y = neglog10padj, color = significant_class)
      ) +
        ggplot2::annotate(
          "rect",
          xmin = lfc_cutoff, xmax = Inf,
          ymin = -log10(padj_cutoff), ymax = Inf,
          fill = scales::alpha(col_up, 0.15),
          color = NA
        ) +
        ggplot2::annotate(
          "rect",
          xmin = -Inf, xmax = -lfc_cutoff,
          ymin = -log10(padj_cutoff), ymax = Inf,
          fill = scales::alpha(col_down, 0.15),
          color = NA
        ) +
        ggplot2::geom_point(
          data = df[df$significant_class == "NS", , drop = FALSE],
          color = col_ns,
          alpha = 0.35,
          size = point_size
        ) +
        ggplot2::geom_point(
          data = df[df$significant_class != "NS", , drop = FALSE],
          alpha = 0.9,
          size = point_size
        ) +
        ggplot2::geom_vline(
          xintercept = c(-lfc_cutoff, lfc_cutoff),
          linetype = "dashed",
          color = "grey30",
          linewidth = 0.5
        ) +
        ggplot2::geom_hline(
          yintercept = -log10(padj_cutoff),
          linetype = "dashed",
          color = "grey30",
          linewidth = 0.5
        ) +
        ggplot2::annotate(
          "text",
          x = max_fc * 0.95,
          y = -log10(padj_cutoff) + 0.15,
          label = paste0("padj = ", format(padj_cutoff, trim = TRUE)),
          hjust = 1,
          size = label_size * 0.8,
          color = "grey40"
        ) +
        ggrepel::geom_text_repel(
          data = top_genes,
          ggplot2::aes(label = gene_id),
          size = label_size,
          max.overlaps = Inf,
          box.padding = 0.4,
          segment.color = "grey50",
          show.legend = FALSE
        ) +
        ggplot2::scale_color_manual(
          values = c(
            "Up" = col_up,
            "Down" = col_down,
            "Sig (low FC)" = col_sig_low,
            "NS" = col_ns
          )
        ) +
        ggplot2::labs(
          title = plot_title,
          subtitle = paste0("↑ Up: ", n_up, "  |  ↓ Down: ", n_down),
          x = "log2 Fold Change",
          y = "-log10(padj)",
          color = NULL
        ) +
        ggplot2::theme_bw(base_size = 12) +
        ggplot2::theme(
          legend.position = "bottom",
          plot.title = ggplot2::element_text(face = "bold", size = title_size),
          plot.subtitle = ggplot2::element_text(color = "grey40", size = subtitle_size),
          axis.text.x = ggplot2::element_text(size = axis_size),
          axis.text.y = ggplot2::element_text(size = axis_size),
          axis.title.x = ggplot2::element_text(size = axis_size + 1),
          axis.title.y = ggplot2::element_text(size = axis_size + 1),
          legend.text = ggplot2::element_text(size = axis_size - 1)
        )
    })
    
    output$open_button_ui <- renderUI({
      req(de_results_r())
      
      if (is.null(de_results_r()) || nrow(de_results_r()) == 0) {
        return(NULL)
      }
      
      actionButton(
        ns("open_volcano"),
        "Open volcano plot",
        class = "btn-primary btn-sm",
        width = "100%"
      )
    })
    
    observeEvent(input$open_volcano, {
      showModal(
        modalDialog(
          title = "Volcano plot",
          size = "l",
          easyClose = TRUE,
          footer = tagList(
            downloadButton(ns("download_volcano_png"), "Save PNG"),
            downloadButton(ns("download_volcano_pdf"), "Save PDF"),
            modalButton("Close")
          ),

          # Same width as the Venn modal. The rule is global (.modal-dialog), and
          # venn_module_v2.r injects an identical one, so the two are consistent
          # by construction rather than by two numbers that can drift apart.
          tags$style(HTML("
      .modal-dialog {
        width: 85vw !important;
        max-width: 85vw !important;
        margin: 1.00rem auto !important;
      }
      .modal-content {
        width: 100% !important;
      }
    ")),

          div(
            style = paste(
              "width: 80vw;",
              "max-width: 80vw;",
              "max-height: 88vh;",
              "overflow-y: auto;",
              "overflow-x: hidden;",
              "padding-right: 8px;"
            ),

            # ---- Always visible: the settings that change WHICH points are
            # Up / Down / NS, and which genes get named. These are not styling:
            # classify_volcano() reads the two cutoffs, so they decide the
            # classification of every point and the subtitle's up/down counts.
            fluidRow(
              column(4, numericInput(ns("volcano_padj_cutoff"), "padj cutoff",
                                     value = 0.05, min = 0, max = 1, step = 0.01)),
              column(4, numericInput(ns("volcano_lfc_cutoff"), "log2FC cutoff",
                                     value = 0.5, min = 0, step = 0.1)),
              column(4, numericInput(ns("volcano_n_labels"), "Number of gene labels",
                                     value = 20, min = 0, step = 1))
            ),

            # ---- Everything cosmetic, folded and grouped ----------------------
            # Same controls, same inputIds, same defaults -- regrouped from one
            # 13-deep stack in a narrow 4-wide column into four themed rows, and
            # put inside the app's standard collapsible_panel() (global.R), the
            # same fold the GO/KEGG "Plot style" blocks and the Venn modal use.
            collapsible_panel(
              "🎨 Appearance & plot size",

              # title and the two element sizes that are not text
              fluidRow(
                column(6, textInput(ns("volcano_title"), "Title", value = "Volcano plot")),
                column(3, sliderInput(ns("volcano_point_size"), "Point size",
                                      min = 0.5, max = 6, value = 2.2, step = 0.1)),
                column(3, sliderInput(ns("volcano_label_size"), "Gene label size",
                                      min = 2, max = 10, value = 4.5, step = 0.1))
              ),

              # the four point classes, in the order the legend shows them
              fluidRow(
                column(3, colourpicker::colourInput(ns("volcano_col_up"),
                                                    "Up", value = "#cc0000")),
                column(3, colourpicker::colourInput(ns("volcano_col_down"),
                                                    "Down", value = "#0055cc")),
                column(3, colourpicker::colourInput(ns("volcano_col_sig_low"),
                                                    "Sig (low FC)", value = "#ff9900")),
                column(3, colourpicker::colourInput(ns("volcano_col_ns"),
                                                    "NS", value = "grey75"))
              ),

              # every text size together, largest to smallest
              fluidRow(
                column(4, sliderInput(ns("volcano_title_size"), "Title size",
                                      min = 10, max = 30, value = 18, step = 1)),
                column(4, sliderInput(ns("volcano_subtitle_size"), "Subtitle size",
                                      min = 8, max = 24, value = 14, step = 1)),
                column(4, sliderInput(ns("volcano_axis_size"), "Axis text size",
                                      min = 8, max = 24, value = 14, step = 1))
              ),

              # canvas size, which is also what the PNG/PDF exports use
              fluidRow(
                column(6, sliderInput(ns("plot_width"), "Plot width (px)",
                                      min = 600, max = 1800, value = 750, step = 50,
                                      width = "100%")),
                column(6, sliderInput(ns("plot_height"), "Plot height (px)",
                                      min = 400, max = 1400, value = 700, step = 50,
                                      width = "100%"))
              ),
              open = FALSE
            ),

            uiOutput(ns("volcano_plot_ui"))
          )
        )
      )
    })
    
    output$volcano_plot_ui <- renderUI({
      req(input$plot_width, input$plot_height)
      
      div(
        style = paste0(
          "width: 100%; overflow-x: auto;"
        ),
        div(
          style = paste0(
            "width: ", input$plot_width, "px;"
          ),
          plotOutput(
            ns("volcano_plot"),
            width = input$plot_width,
            height = paste0(input$plot_height, "px")
          )
        )
      )
    })
    
    output$volcano_plot <- renderPlot({
      req(volcano_plot_r())
      volcano_plot_r()
    }, res = 120)
    
    output$download_volcano_png <- downloadHandler(
      filename = function() {
        paste0("volcano_plot_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png")
      },
      content = function(file) {
        ggplot2::ggsave(
          filename = file,
          plot = volcano_plot_r(),
          width = (if (!is.null(input$plot_width)) input$plot_width else 1000) / 100,
          height = (if (!is.null(input$plot_height)) input$plot_height else 700) / 100,
          dpi = 300,
          bg = "white"
        )
      }
    )
    
    output$download_volcano_pdf <- downloadHandler(
      filename = function() {
        paste0("volcano_plot_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf")
      },
      content = function(file) {
        ggplot2::ggsave(
          filename = file,
          plot = volcano_plot_r(),
          width = (if (!is.null(input$plot_width)) input$plot_width else 1000) / 100,
          height = (if (!is.null(input$plot_height)) input$plot_height else 700) / 100,
          device = cairo_pdf,
          bg = "white"
        )
      }
    )
  })
}