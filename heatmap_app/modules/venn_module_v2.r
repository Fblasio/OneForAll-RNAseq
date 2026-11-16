# =============================================================================
# venn_module.R
# -----------------------------------------------------------------------------
# Interactive Venn diagram from DE gene-list selections.
# Rendered inside a modal dialog triggered by an external button.
#
# Input  (server): de_history_rv — the $history reactiveVal from
#                  selectionHistoryServer("de_hist").
#                  Each list entry: list(name, genes, n, time)
#
# Usage in app.R:
#   UI     : vennModalUI("venn")          # declares the modal shell + trigger
#   Server : vennServer("venn", de_history_rv = de_hist_out$history)
#   Trigger: actionButton("open_venn", ...) anywhere in the UI
#            + observeEvent in server (see app.R instructions)
#
# Packages: ggvenn (plot), ggVennDiagram (region genes), ggplot2, DT
# =============================================================================

# -----------------------------------------------------------------------------
# UI — only the modal definition; the trigger button lives in app.R
# -----------------------------------------------------------------------------
vennModalUI <- function(id) {
  ns <- NS(id)
  # Modal is injected into the DOM at startup; shown/hidden by showModal()
  tags$div(
    # Hidden UI container — content rendered server-side via uiOutput
    uiOutput(ns("venn_modal_content"))
  )
}

# -----------------------------------------------------------------------------
# SERVER
# -----------------------------------------------------------------------------
vennServer <- function(id, de_history_rv) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    .venn_placeholder_plot <- function() {
      ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0.5, y = 0.5,
                          label = "No Venn diagram drawn yet.\nSelect 2\u20135 lists and click \"Draw Venn\".",
                          size = 5, color = "grey40") +
        ggplot2::theme_void()
    }

    # -----------------------------------------------------------------------
    # Show modal when the external trigger fires.
    # The trigger is input$open_venn (defined in app.R outside the module),
    # so it is accessed via session$rootScope via the parent session —
    # simplest: expose a callable that app.R's observeEvent can call.
    # -----------------------------------------------------------------------
    show_modal <- function() {
      showModal(modalDialog(
        title = tagList(
          tags$span("🔵 Venn Diagram — DE gene lists"),
          tags$small(
            style = "color:#888; margin-left:10px; font-weight:normal;"
          )
        ),
        size = "l",
        easyClose = TRUE,
        footer = modalButton("Close"),
        
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

        # ---- Selector ----
        uiOutput(ns("list_selector_ui")),
        br(),

        # ---- Action row ----
        fluidRow(
          column(
            4,
            actionButton(
              ns("draw_venn"),
              "🔵 Draw Venn",
              class = "btn-primary btn-sm",
              width = "100%"
            )
          ),
          column(
            4,
            downloadButton(
              ns("dl_venn_png"),
              "⬇ Save PNG",
              class = "btn-default btn-sm",
              width = "100%"
            )
          )
        ),
        br(),
        
        # ---- Personalize venn (folded) ----
        # Purely a layout change: every control below is unchanged -- same
        # inputIds, same defaults, same order -- they are only moved inside a
        # <details> disclosure that starts CLOSED, so this modal opens on the
        # three things the task actually needs (tick the lists, Draw, Save PNG)
        # instead of on ten sliders. collapsible_panel() is the app's existing
        # fold from global.R, the same construction the GO and KEGG "Plot style"
        # blocks use, so this introduces no new mechanism and no new package.
        collapsible_panel(
          "🎨 Appearance & plot size",
          h5("🎨 Set colours", style = "margin-bottom: 6px;"),
          uiOutput(ns("venn_set_colour_ui")),
          sliderInput(ns("venn_alpha"), "Fill opacity", min = 0.05, max = 1, value = 0.75, step = 0.05),

          fluidRow(
            column(
              3,
              sliderInput(ns("venn_set_size"), "Set text", min = 3, max = 10, value = 4, step = 0.5)
            ),
            column(
              3,
              sliderInput(ns("venn_label_size"), "Region text", min = 1, max = 8, value = 2, step = 0.5)
            ),
            column(
              3,
              sliderInput(ns("venn_edge_size"), "Edge width", min = 0.3, max = 2.5, value = 0.6, step = 0.1)
            )
          ),
          fluidRow(
            column(
              4,
              selectInput(
                ns("venn_label_mode"),
                "Region labels",
                choices = c("Count" = "count", "Percent" = "percent", "Both" = "both"),
                selected = "both"
              )
            ),
            column(
              4,
              sliderInput(ns("venn_percent_digits"), "Percent digits", min = 0, max = 2, value = 1, step = 1)
            ),
            column(
              4,
              colourpicker::colourInput(
                ns("venn_label_color"),
                "Region text color",
                value = "#111111"
              )
            )
          ),

          # ---- Plot size ----
          fluidRow(
            column(
              6,
              sliderInput(
                ns("venn_plot_width"),
                "Plot width (px)",
                min = 700,
                max = 1600,
                value = 875,
                step = 25
              )
            ),
            column(
              6,
              sliderInput(
                ns("venn_plot_height"),
                "Plot height (px)",
                min = 500,
                max = 1200,
                value = 725,
                step = 25
              )
            )
          ),
          open = FALSE
        ),

        # ---- Region feedback ----
        uiOutput(ns("venn_plot_ui")),
        uiOutput(ns("region_selector_ui")),
        uiOutput(ns("region_header_ui")),
        DT::DTOutput(ns("region_table")),
        uiOutput(ns("region_dl_ui"))
        )
      ))
    }


    # -----------------------------------------------------------------------
    # Flat named list: name → character vector of genes
    # -----------------------------------------------------------------------
    flat_sets <- reactive({
      h <- de_history_rv()
      if (length(h) == 0) return(list())

      out    <- list()
      for (entry in h) {
        nm  <- trimws(entry$name  %||% "Unnamed")
        gns <- as.character(entry$genes)
        if (length(gns) == 0) next
        key    <- nm
        suffix <- 1L
        while (key %in% names(out)) {
          key    <- paste0(nm, "_", suffix)
          suffix <- suffix + 1L
        }
        out[[key]] <- gns
      }
      out
    })

    # -----------------------------------------------------------------------
    # Checkbox selector
    # -----------------------------------------------------------------------
    output$list_selector_ui <- renderUI({
      fl <- flat_sets()
      if (length(fl) == 0) {
        return(helpText("⚠️ No DE gene lists saved yet. Run DE analysis and save selections first."))
      }
      labels <- paste0(names(fl), " (", sapply(fl, length), " genes)")
      vals <- names(fl)
      
      checkboxGroupInput(
        inputId = ns("chosen_sets"),
        label = "DE lists to compare:",
        choiceNames = as.list(labels),
        choiceValues = as.list(vals),
        selected = vals[seq_len(min(length(vals), 3L))]
      )
    })

    # -----------------------------------------------------------------------
    # Validated sets (2–5) for venn diagramm
    # -----------------------------------------------------------------------
      valid_sets <- reactive({
      fl <- flat_sets()
      req(length(fl) > 0)
      req(length(input$chosen_sets) > 0)
      
      sel <- intersect(input$chosen_sets, names(fl))
      
      shiny::validate(
        shiny::need(length(sel) >= 2, "⚠️ Select at least 2 lists."),
        shiny::need(length(sel) <= 5, "⚠️ Maximum 5 sets for a Venn diagram.")
      )
      
      fl[sel]
    })
      # dimension plot
      plot_dims <- reactive({
        list(
          width_px = input$venn_plot_width %||% 875,
          height_px = input$venn_plot_height %||% 725,
          res = 200
        )
      })
       # ---- style to personalize 
      style_inputs <- reactive({
        sets <- tryCatch(valid_sets(), error = function(e) NULL)
        n <- if (is.null(sets)) 0L else length(sets)
        
        default_cols <- c(
          "#1F77B4", "#D62728", "#2CA02C",
          "#FF7F0E", "#9467BD", "#8C564B"
        )
        
        set_cols <- if (n == 0) {
          character(0)
        } else {
          stats::setNames(vapply(seq_len(n), function(i) {
            nm_safe <- make.names(names(sets)[i])
            col_i   <- input[[paste0("venn_set_fill__", nm_safe)]]
            if (is.null(col_i) || !nzchar(col_i)) {
              default_cols[((i - 1) %% length(default_cols)) + 1]
            } else col_i
          }, character(1)), names(sets))
        }
        
        list(
          set_cols = set_cols,
          alpha = input$venn_alpha %||% 0.75,
          edge_lwd = input$venn_edge_size %||% 0.8,
          set_textsize = input$venn_set_size %||% 5,
          region_textsize = input$venn_label_size %||% 4.5,
          label_mode = input$venn_label_mode %||% "both",
          pct_digits = input$venn_percent_digits %||% 1,
          label_col = input$venn_label_color %||% "#111111"
        )
      })
      
      output$venn_set_colour_ui <- renderUI({
        sets <- tryCatch(valid_sets(), error = function(e) NULL)
        if (is.null(sets) || length(sets) == 0) return(NULL)
        
        default_cols <- c(
          "#1F77B4", "#D62728", "#2CA02C",
          "#FF7F0E", "#9467BD", "#8C564B"
        )
        
        n <- length(sets)
        
        tagList(
          lapply(seq_len(n), function(i) {
            nm_safe <- make.names(names(sets)[i])
            colourpicker::colourInput(
              inputId = ns(paste0("venn_set_fill__", nm_safe)),
              label = paste(names(sets)[i], "colour"),
              value = default_cols[((i - 1) %% length(default_cols)) + 1]
            )
          })
        )
      })

      # MUST come after the output above is defined -- outputOptions() looks the
      # name up in the session's output registry and stops with "... is not in
      # list of output objects" if it is called first.
      #
      # Why it is needed at all: this output now lives inside a
      # collapsible_panel() that starts CLOSED, and a closed <details> is
      # display:none, which Shiny treats as hidden and therefore suspends. Without
      # this the per-set colour pickers would not be created until the panel was
      # first expanded. The draw code already falls back to the same default_cols
      # when an input is NULL, so the drawn figure is identical either way; this
      # keeps the folded and unfolded states carrying the same set of inputs.
      outputOptions(output, "venn_set_colour_ui", suspendWhenHidden = FALSE)

      output$venn_plot_ui <- renderUI({
        dims <- plot_dims()
        
        plotOutput(
          outputId = ns("venn_plot"),
          width = "100%",
          height = paste0(dims$height_px, "px")
        )
      })
      
    # -----------------------------------------------------------------------
    # Venn plot — triggered by button
    # -----------------------------------------------------------------------
      venn_base <- eventReactive(input$draw_venn, {
        sets <- valid_sets()
        
        shiny::validate(
          shiny::need(length(sets) >= 2, "⚠️ Select at least 2 lists."),
          shiny::need(length(sets) <= 5, "⚠️ Maximum 5 sets for a Venn diagram.")
        )
        
        sets
      }, ignoreNULL = FALSE)
      
      build_venn_plot <- reactive({
        sets <- venn_base()
        req(length(sets) >= 2)
        
        st <- style_inputs()
        set_cols <- st$set_cols[names(sets)]
        set_cols[is.na(set_cols)] <- "#999999"   # fallback: no picker rendered yet for this name
        
        p <- ggvenn::ggvenn(
          data = sets,
          fill_color = set_cols,
          fill_alpha = st$alpha,
          stroke_color = set_cols,
          stroke_size = st$edge_lwd,
          set_name_color = set_cols,
          set_name_size = st$set_textsize,
          text_color = st$label_col,
          text_size = st$region_textsize,
          show_percentage = identical(st$label_mode, "percent") || identical(st$label_mode, "both"),
          show_count = identical(st$label_mode, "count") || identical(st$label_mode, "both"),
          digits = st$pct_digits
        ) +
          ggplot2::scale_x_continuous(
            expand = ggplot2::expansion(mult = 0.02)
          ) +
          ggplot2::scale_y_continuous(
            expand = ggplot2::expansion(mult = 0.02)
          ) +
          ggplot2::theme_void() +
          ggplot2::theme(
            legend.position = "none",
            plot.margin = ggplot2::margin(t = 2, r = 2, b = 2, l = 2),
            plot.background = ggplot2::element_rect(fill = "white", colour = "white"),
            panel.background = ggplot2::element_rect(fill = "white", colour = "white")
          )
        
        p
      })
      
      # output START
      output$venn_plot <- renderPlot(
        {
          req(build_venn_plot())
          build_venn_plot()
        },
        width = function() plot_dims()$width_px,
        height = function() plot_dims()$height_px,
        res = 200,
        execOnResize = TRUE
      )

    # -----------------------------------------------------------------------
    # Region genes via ggVennDiagram::process_region_data (dropdown-driven).
    # -----------------------------------------------------------------------
    region_data <- reactive({
      sets <- valid_sets()
      vd <- ggVennDiagram::Venn(sets)
      rdata <- ggVennDiagram::process_region_data(vd)
      
      rdata <- as.data.frame(rdata)
      rdata <- rdata[rdata$count > 0, , drop = FALSE]
      
      if (nrow(rdata) == 0) {
        return(NULL)
      }
      
      rdata
    })
    
    output$region_selector_ui <- renderUI({
      rd <- region_data()
      if (is.null(rd) || nrow(rd) == 0) {
        return(helpText("No non-empty intersections available."))
      }
      
      labels <- paste0(rd$name, " (", rd$count, " genes)")
      
      selectInput(
        inputId = ns("selected_region"),
        label = "Intersection to inspect:",
        choices = stats::setNames(as.character(rd$id), labels),
        selected = as.character(rd$id[1]),
        width = "100%"
      )
    })
    
    clicked_region <- reactive({
      rd <- region_data()
      req(!is.null(rd), nrow(rd) > 0)
      req(input$selected_region)
      
      hit <- which(as.character(rd$id) == input$selected_region)
      req(length(hit) == 1)
      
      list(
        label = as.character(rd$name[[hit]]),
        genes = as.character(unlist(rd$item[[hit]]))
      )
    })

    # -----------------------------------------------------------------------
    # Region header
    # -----------------------------------------------------------------------
    output$region_header_ui <- renderUI({
      cr <- clicked_region()
      if (is.null(cr)) {
        return(helpText("👆 Select an intersection to inspect its genes."))
      }
      tagList(
        hr(),
        tags$p(
          style = "font-weight:bold; color:#1565a8; font-size:1.05em;",
          paste0("📍 Region: ", cr$label, "  —  ", length(cr$genes), " gene(s)")
        )
      )
    })

    # -----------------------------------------------------------------------
    # Region gene table
    # -----------------------------------------------------------------------
    output$region_table <- DT::renderDT({
      cr <- clicked_region()
      req(!is.null(cr), length(cr$genes) > 0)
      DT::datatable(
        data.frame(GeneID = cr$genes, stringsAsFactors = FALSE),
        rownames   = FALSE,
        extensions = "Buttons",
        options    = list(
          dom        = "Btp",
          pageLength = 25,
          buttons    = list(
            list(
              extend        = "copy",
              title         = NULL,
              header        = FALSE,
              exportOptions = list(
                columns  = 0,
                modifier = list(page = "all", selected = NULL)
              )
            ),
            list(
              extend        = "csv",
              exportOptions = list(modifier = list(page = "all", selected = NULL))
            )
          )
        )
      )
    }, server = FALSE)

    # -----------------------------------------------------------------------
    # Region download (.txt, one gene per line)
    # -----------------------------------------------------------------------
    output$region_dl_ui <- renderUI({
      cr <- clicked_region()
      if (is.null(cr) || length(cr$genes) == 0) return(NULL)
      tagList(
        br(),
        downloadButton(
          ns("dl_region_txt"),
          paste0("⬇ Download region genes (", length(cr$genes), ")"),
          class = "btn-sm btn-success"
        )
      )
    })

    output$dl_region_txt <- downloadHandler(
      filename = function() {
        cr   <- tryCatch(clicked_region(), error = function(e) NULL)
        safe <- if (is.null(cr)) "none" else gsub("[^A-Za-z0-9_]", "_", cr$label %||% "region")
        paste0("venn_region_", safe, "_",
               format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt")
      },
      content = function(file) {
        cr <- tryCatch(clicked_region(), error = function(e) NULL)
        if (is.null(cr)) {
          writeLines("No intersection selected \u2014 draw the Venn diagram first.", file)
        } else {
          writeLines(cr$genes, file)
        }
      }
    )

    # -----------------------------------------------------------------------
    # PNG download
    # -----------------------------------------------------------------------
    output$dl_venn_png <- downloadHandler(
      filename = function() {
        paste0("venn_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png")
      },
      content = function(file) {
        p    <- tryCatch(build_venn_plot(), error = function(e) .venn_placeholder_plot())
        dims <- plot_dims()
        
        ggplot2::ggsave(
          filename = file,
          plot = p,
          width = dims$width_px,
          height = dims$height_px,
          units = "px",
          dpi = 200,
          bg = "white",
          limitsize = FALSE
        )
      }
    )
    
    # Exposed so app.R can call it
    return(list(
      show = show_modal
    ))
    
  })
}
