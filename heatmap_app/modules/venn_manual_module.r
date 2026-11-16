# =============================================================================
# venn_manual_module.R
# Independent manual Venn diagram module for the application.
# Two UI entry points:
#   - vennManualSidebarUI(id)  → placed in sidebarPanel
#   - vennManualMainUI(id)     → placed in mainPanel
# Server:
#   - vennManualServer(id)
# Does NOT touch venn_module_v2.r or any other existing module.
# =============================================================================

# ── Sidebar UI: all style / personalisation controls ─────────────────────────

vennManualSidebarUI <- function(id) {
  ns <- NS(id)

  tagList(
    h4("\U0001f7e2 Venn diagram", style = "margin-bottom: 4px;"),
    tags$p(
      style = "font-size: 0.85em; color: #666; margin-bottom: 8px;",
      "Paste gene lists in the main panel."
    ),

    # ── Set management ──────────────────────────────────────────────────────
    fluidRow(
      column(6, actionButton(ns("addset"),    "\u2795 Add set",
                             class = "btn-success btn-sm", width = "100%")),
      column(6, actionButton(ns("removeset"), "\u2796 Remove set",
                             class = "btn-warning btn-sm", width = "100%"))
    ),
    tags$small(uiOutput(ns("setcount_label"))),

    hr(),

    # ── Action buttons ──────────────────────────────────────────────────────
    actionButton(ns("drawvenn"), "\U0001f504 Draw / Refresh Venn",
                 class = "btn-primary btn-sm", width = "100%"),
    br(), br(),
    fluidRow(
      column(6, downloadButton(ns("dlvennpng"), "\u2b07 PNG",
                               class = "btn-default btn-sm", width = "100%")),
      column(6, downloadButton(ns("dlvennpdf"), "\u2b07 PDF",
                               class = "btn-default btn-sm", width = "100%"))
    ),
    br(),
    actionButton(ns("resetall"), "\U0001f9f9 Clear all lists",
                 class = "btn-danger btn-sm", width = "100%"),

    hr(),

    # ── Fill colours ────────────────────────────────────────────────────────
    h5("🎨 Set colours", style = "margin-bottom: 6px;"),
    uiOutput(ns("set_colour_ui")),
    sliderInput(ns("alpha"), "Fill opacity",
                min = 0.05, max = 1, value = 0.75, step = 0.05),

    hr(),

    # ── Edge / set labels ────────────────────────────────────────────────────
    h5("\u2712\ufe0f Labels & edges", style = "margin-bottom: 6px;"),
    fluidRow(
      column(6, colourpicker::colourInput(ns("label_col"), "Region text col.", value = "#111111"))
    ),
    fluidRow(
      column(6, sliderInput(ns("edge_lwd"), "Edge width", min = 0.2, max = 3, value = 0.7, step = 0.1)),
      column(6, sliderInput(ns("set_textsize"), "Set label size", min = 2, max = 12, value = 5, step = 0.5))
    ),
    sliderInput(ns("region_textsize"), "Region text size",
                min = 2, max = 10, value = 4, step = 0.5),

    hr(),

    # ── Region label format ──────────────────────────────────────────────────
    h5("\U0001f522 Region label content", style = "margin-bottom: 6px;"),
    selectInput(ns("label_mode"), NULL,
                choices  = c("Count"   = "count",
                             "Percent" = "percent",
                             "Both"    = "both"),
                selected = "both",
                width    = "100%"),
    sliderInput(ns("pct_digits"), "Percent digits",
                min = 0, max = 2, value = 1, step = 1),

    hr(),

    # ── Background & layout ──────────────────────────────────────────────────
    h5("\U0001f4d0 Layout", style = "margin-bottom: 6px;"),
    colourpicker::colourInput(ns("bg_col"), "Background", value = "#FFFFFF"),
    fluidRow(
      column(6, sliderInput(ns("plot_h"), "Plot height (px)",
                            min = 300, max = 1200, value = 560, step = 25)),
      column(6, sliderInput(ns("plot_w"), "Plot width (px)",
                            min = 400, max = 1400, value = 700, step = 25))
    ),
    sliderInput(ns("expand_xy"), "Axis expansion (labels margin)",
                min = 0.1, max = 0.6, value = 0.30, step = 0.05)
  )
}

# ── Main panel UI: gene-list input boxes + plot output ───────────────────────

vennManualMainUI <- function(id) {
  ns <- NS(id)

  tagList(
    # Header row
    fluidRow(
      column(
        12,
        div(
          style = "display: flex; align-items: center; gap: 12px; margin-bottom: 12px;",
          h3("\U0001f7e2 Venn diagram", style = "margin: 0;"),
          tags$span(
            style = "font-size: 0.9em; color: #555;",
            # 5, not 6: ggvenn::ggvenn() -- the renderer used below -- refuses
            # more than 5 sets ("n_sets >= min_set_num && n_sets <= max_set_num
            # is not TRUE"), so a 6th box made "Draw / Refresh Venn" throw.
            # venn_module_v2.r already capped at 5; the two now agree.
            "Paste one gene list per box \u2014 min 2, max 5 sets"
          )
        )
      )
    ),

    # Dynamic gene-list boxes
    uiOutput(ns("genebox_ui")),

    hr(),

    # Plot area
    uiOutput(ns("venn_plot_container")),

    hr(),

    # Intersection inspector
    uiOutput(ns("inspector_ui"))
  )
}


# ── Pure region computation (extracted so it's testable in isolation) ────────
# sets: named list of character vectors (one per gene list). Returns a
# data.frame with one row per non-empty EXCLUSIVE region (id, name, count,
# item = list-column of gene IDs in that region only). No Shiny dependency.
compute_venn_regions <- function(sets) {
  if (length(sets) < 2) return(NULL)
  
  set_names <- names(sets)
  n <- length(sets)
  
  all_combinations <- unlist(
    lapply(seq_len(n), function(k) combn(set_names, k, simplify = FALSE)),
    recursive = FALSE
  )
  
  res_list <- lapply(seq_along(all_combinations), function(i) {
    comb <- all_combinations[[i]]
    
    in_genes <- Reduce(intersect, sets[comb])
    
    other_sets <- setdiff(set_names, comb)
    if (length(other_sets) > 0) {
      other_genes <- unique(unlist(sets[other_sets], use.names = FALSE))
      in_genes <- setdiff(in_genes, other_genes)
    }
    
    in_genes <- sort(unique(in_genes))
    
    data.frame(
      id = paste(comb, collapse = "/"),
      name = paste(comb, collapse = " \u2229 "),
      count = length(in_genes),
      stringsAsFactors = FALSE
    ) |>
      transform(item = I(list(in_genes)))
  })
  
  rdata <- do.call(rbind, res_list)
  rdata <- rdata[rdata$count > 0, , drop = FALSE]
  
  if (nrow(rdata) == 0) return(NULL)
  rownames(rdata) <- NULL
  rdata
}

# ── Server ─────────────────────────────────────────────────────────────────────

vennManualServer <- function(id) {
  moduleServer(id, function(input, output, session) {

    ns <- session$ns
    
    # add one placeholder-plot helper, and make each download tolerant of the "nothing drawn yet" state
    .venn_placeholder_plot <- function() {
      ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0.5, y = 0.5,
                          label = "No Venn diagram drawn yet.\nFill in \u2265 2 gene lists and click \"Draw / Refresh Venn\".",
                          size = 5, color = "grey40") +
        ggplot2::theme_void()
    }

    # -- Reactive state --------------------------------------------------------
    nsets <- reactiveVal(2L)

    # Hard cap at 5: ggvenn::ggvenn() cannot draw 6 (see vennManualMainUI).
    observeEvent(input$addset,    { if (nsets() < 5L) nsets(nsets() + 1L) })
    observeEvent(input$removeset, { if (nsets() > 2L) nsets(nsets() - 1L) })

    output$setcount_label <- renderUI({
      tags$span(
        style = "color: #555;",
        paste0("Sets: ", nsets(), " / 5")
      )
    })
    
    output$set_colour_ui <- renderUI({
      n <- nsets()
      
      default_cols <- c(
        "#1F77B4", "#D62728", "#2CA02C",
        "#FF7F0E", "#9467BD", "#8C564B"
      )
      
      tagList(
        lapply(seq_len(n), function(i) {
          colourpicker::colourInput(
            inputId = ns(paste0("set_fill_", i)),
            label = paste("Set", i, "colour"),
            value = default_cols[i]
          )
        })
      )
    })

    # Reset: clear all text areas by incrementing a counter that forces re-render
    reset_counter <- reactiveVal(0L)
    
    observeEvent(input$resetall, {
      n <- nsets()
      
      for (i in seq_len(n)) {
        updateTextInput(
          session = session,
          inputId = paste0("setname_", i),
          value = paste0("Set_", i)
        )
        updateTextAreaInput(
          session = session,
          inputId = paste0("setgenes_", i),
          value = ""
        )
      }
    })
    
    # Clear individual set inputs
    lapply(seq_len(6), function(i) {
      observeEvent(input[[paste0("clear_set_", i)]], {
        req(i <= nsets())
        
        updateTextInput(
          session = session,
          inputId = paste0("setname_", i),
          value = paste0("Set_", i)
        )
        updateTextAreaInput(
          session = session,
          inputId = paste0("setgenes_", i),
          value = ""
        )
      }, ignoreInit = TRUE)
    })

    # -- Dynamic gene-list boxes -----------------------------------------------
    output$genebox_ui <- renderUI({

      force(reset_counter())   # re-renders when reset is clicked

      n <- nsets()
      # Colour palette for box headers (mild, publication-friendly)
      header_cols <- c("#2171B5", "#CB181D", "#238B45", "#D94801",
                       "#6A3D9A", "#7B4173")

      fluidRow(
        lapply(seq_len(n), function(i) {
          
          box_col <- header_cols[i]
          
          column(
            width = 4, # change how many box in one line
            
            div(
              style = paste0(
                "border: 2px solid ", box_col, "; border-radius: 6px; ",
                "margin-bottom: 14px; overflow: hidden;"
              ),
              
              # Box header
              div(
                style = paste0(
                  "background: ", box_col, "; color: #fff; ",
                  "padding: 6px 12px; display: flex; ",
                  "align-items: center; gap: 10px;"
                ),
                tags$span(
                  style = "font-weight: bold; font-size: 0.95em;",
                  paste("Set", i)
                ),
                div(
                  style = "flex: 1;",
                  textInput(
                    inputId = ns(paste0("setname_", i)),
                    label = NULL,
                    # isolate(): re-read the CURRENT value so a legitimate
                    # re-render (Add set / Remove set / reset) preserves what the
                    # user typed -- WITHOUT making this renderUI depend on the
                    # input it is creating. Without the isolate, every keystroke
                    # invalidated genebox_ui, rebuilt all the boxes as fresh HTML,
                    # and the browser replaced the <input> being typed into, so
                    # focus was lost after ~1s and the user had to click again.
                    # Same idiom and same reason as app.r's annotation colour
                    # pickers (app.r:1210).
                    value = isolate(
                      if (!is.null(input[[paste0("setname_", i)]])) {
                        input[[paste0("setname_", i)]]
                      } else {
                        paste0("Set_", i)
                      }
                    ),
                    placeholder = "Set name",
                    width = "100%"
                  )
                )
              ),
              
              # Text area
              div(
                style = "padding: 8px;",
                textAreaInput(
                  inputId = ns(paste0("setgenes_", i)),
                  label = NULL,
                  # isolate() for the same reason as the set-name box above: this
                  # textarea had the identical focus-loss problem while typing.
                  value = isolate(
                    if (!is.null(input[[paste0("setgenes_", i)]])) {
                      input[[paste0("setgenes_", i)]]
                    } else {
                      ""
                    }
                  ),
                  placeholder = paste0(
                    "Paste gene IDs here (one per line, or comma-separated)\n",
                    "Example:\nGENE1\nGENE2\nGENE3"
                  ),
                  rows = 4,
                  width = "100%"
                ),
                actionButton(
                  inputId = ns(paste0("clear_set_", i)),
                  label = "Clear",
                  class = "btn btn-default btn-xs",
                  width = "100%"
                )
              )
            )
          )
        })
      )
    })

    # -- Parse gene lists ------------------------------------------------------
    parsed_sets <- reactive({
      n   <- nsets()
      out <- list()

      for (i in seq_len(n)) {
        nm      <- input[[paste0("setname_", i)]]
        raw_txt <- input[[paste0("setgenes_", i)]]

        nm <- trimws(
          if (is.null(nm) || !nzchar(trimws(nm))) paste0("Set_", i) else nm
        )

        genes <- unlist(
          strsplit(
            x     = if (is.null(raw_txt)) "" else raw_txt,
            split = "[,\n\r\t ;]+"
          ),
          use.names = FALSE
        )
        genes <- trimws(genes)
        genes <- canon_ids(genes)
        genes <- unique(genes)
        genes <- genes[nzchar(genes)]

        if (length(genes) > 0) {
          # Deduplicate set names
          key    <- nm
          suffix <- 1L
          while (key %in% names(out)) {
            key    <- paste0(nm, "_", suffix)
            suffix <- suffix + 1L
          }
          attr(genes, "box_index") <- i   # remember which colour box this came from
          out[[key]] <- genes
        }
      }
      out
    })

    valid_sets <- reactive({
      s <- parsed_sets()
      s[vapply(s, function(x) length(x) >= 1L, logical(1))]
    })

    # -- Build ggplot Venn (reactive, redraws on style changes too) -----------
    # eventReactive on drawvenn for the gene lists; style inputs are observed
    # via a separate reactive so that changing colours/sizes refreshes live.

    venn_base <- eventReactive(input$drawvenn, {
      sets <- valid_sets()
      shiny::validate(
        shiny::need(length(sets) >= 2,
                    "Please fill in at least 2 gene lists, then click Draw / Refresh Venn."),
        # Belt and braces behind the nsets() cap above: if a 6th set ever
        # reaches here, say so plainly instead of letting ggvenn() throw a
        # raw "n_sets >= min_set_num" error into the plot panel.
        shiny::need(length(sets) <= 5,
                    "A Venn diagram can show at most 5 sets. Remove a list, or use the intersection table below.")
      )
      sets
    }, ignoreNULL = FALSE)

    # Style inputs collected as a reactive (to allow live updates)
    style_inputs <- reactive({
      n <- nsets()
      
      set_cols <- vapply(seq_len(n), function(i) {
        col_i <- input[[paste0("set_fill_", i)]]
        if (is.null(col_i) || !nzchar(col_i)) {
          c("#1F77B4", "#D62728", "#2CA02C", "#FF7F0E", "#9467BD", "#8C564B")[i]
        } else {
          col_i
        }
      }, character(1))
      
      list(
        set_cols           = set_cols,
        alpha              = input$alpha,
        label_col          = input$label_col,
        edge_lwd           = input$edge_lwd,
        set_textsize       = input$set_textsize,
        region_textsize    = input$region_textsize,
        label_mode         = input$label_mode,
        pct_digits         = input$pct_digits,
        bg_col             = input$bg_col,
        expand_xy          = input$expand_xy
      )
    })

    build_venn_plot <- reactive({
      sets <- venn_base()
      shiny::req(length(sets) >= 2)
      
      st <- style_inputs()
      
      set_names <- names(sets)
      box_idx   <- vapply(sets, function(x) {
        bi <- attr(x, "box_index")
        if (is.null(bi)) NA_integer_ else bi
      }, integer(1))
      set_cols  <- st$set_cols[box_idx]
      
      total <- length(unique(unlist(sets)))
      
      # Build base ggvenn plot
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
        ggplot2::coord_equal() +
        ggplot2::scale_x_continuous(
          expand = ggplot2::expansion(mult = st$expand_xy)
        ) +
        ggplot2::scale_y_continuous(
          expand = ggplot2::expansion(mult = st$expand_xy)
        ) +
        ggplot2::theme_void() +
        ggplot2::theme(
          legend.position = "none",
          plot.margin = ggplot2::margin(20, 25, 20, 25),
          plot.background = ggplot2::element_rect(
            fill = st$bg_col,
            colour = st$bg_col
          ),
          panel.background = ggplot2::element_rect(
            fill = st$bg_col,
            colour = st$bg_col
          )
        )
      p
    })

    # -- Plot container (dynamic height / width) --------------------------------
    output$venn_plot_container <- renderUI({
      shiny::req(venn_base())
      h <- if (!is.null(input$plot_h)) input$plot_h else 560
      plotOutput(ns("vennplot"),
                 width  = "100%",
                 height = paste0(h, "px"))
    })

    output$vennplot <- renderPlot({
      build_venn_plot()
    }, res = 120, execOnResize = TRUE)

    # -- Intersection inspector ────────────────────────────────────────────────
    region_table <- reactive({
      sets <- venn_base()
      shiny::req(length(sets) >= 2)
      compute_venn_regions(sets)
    })

    output$inspector_ui <- renderUI({
      rd <- region_table()
      if (is.null(rd)) return(NULL)

      labels <- paste0(rd$name, " \u2014 ", rd$count, " genes")
      tagList(
        h4("\U0001f50d Inspect intersection"),
        fluidRow(
          column(
            6,
            selectInput(
              inputId  = ns("sel_region"),
              label    = "Select region",
              choices  = stats::setNames(as.character(rd$id), labels),
              selected = as.character(rd$id[1]),
              width    = "100%"
            )
          ),
          column(
            3,
            br(),
            downloadButton(ns("dl_region_txt"),
                           "\u2b07 Gene list (.txt)",
                           class = "btn-success btn-sm",
                           width = "100%")
          )
        ),
        DT::DTOutput(ns("region_dt"))
      )
    })

    clicked_region <- reactive({
      rd <- region_table()
      shiny::req(!is.null(rd), nrow(rd) > 0, input$sel_region)
      hit <- which(as.character(rd$id) == input$sel_region)
      shiny::req(length(hit) == 1)
      list(
        label = as.character(rd$name[hit]),
        genes = as.character(unlist(rd$item[hit]))
      )
    })

    output$region_dt <- DT::renderDT({
      cr <- clicked_region()
      shiny::req(cr, length(cr$genes) > 0)
      
      DT::datatable(
        data.frame(GeneID = sort(cr$genes), stringsAsFactors = FALSE),
        rownames = FALSE,
        extensions = "Buttons",
        options = list(
          dom = "Btp",
          pageLength = 30,
          buttons = list(
            list(
              extend = "copy",
              title = NULL,
              header = FALSE,
              exportOptions = list(
                columns = 0,
                modifier = list(page = "all", selected = NULL)
              )
            ),
            list(
              extend = "csv",
              exportOptions = list(
                modifier = list(page = "all", selected = NULL)
              )
            )
          )
        )
      )
    }, server = FALSE)

    # -- Downloads ─────────────────────────────────────────────────────────────
    output$dl_region_txt <- downloadHandler(
      filename = function() {
        cr   <- tryCatch(clicked_region(), error = function(e) NULL)
        safe <- if (is.null(cr)) "none" else gsub("[^A-Za-z0-9]+", "_", cr$label)
        paste0("venn_region_", safe, "_",
               format(Sys.time(), "%Y%m%d%H%M%S"), ".txt")
      },
      content = function(file) {
        cr <- tryCatch(clicked_region(), error = function(e) NULL)
        if (is.null(cr)) {
          writeLines("No intersection selected \u2014 draw the Venn diagram first.", file)
        } else {
          writeLines(sort(cr$genes), file)
        }
      }
    )

    output$dlvennpng <- downloadHandler(
      filename = function() {
        paste0("venn_", format(Sys.time(), "%Y%m%d%H%M%S"), ".png")
      },
      content = function(file) {
        p <- tryCatch(build_venn_plot(), error = function(e) .venn_placeholder_plot())
        h_in <- max(4, (if (!is.null(input$plot_h)) input$plot_h else 560) / 96)
        w_in <- max(5, (if (!is.null(input$plot_w)) input$plot_w else 700) / 96)
        ggplot2::ggsave(
          filename  = file,
          plot      = p,
          width     = w_in,
          height    = h_in,
          units     = "in",
          dpi       = 400,
          bg        = if (!is.null(input$bg_col)) input$bg_col else "#FFFFFF",
          limitsize = FALSE
        )
      }
    )

    output$dlvennpdf <- downloadHandler(
      filename = function() {
        paste0("venn_", format(Sys.time(), "%Y%m%d%H%M%S"), ".pdf")
      },
      content = function(file) {
        p <- tryCatch(build_venn_plot(), error = function(e) .venn_placeholder_plot())
        h_in <- max(4, (if (!is.null(input$plot_h)) input$plot_h else 560) / 96)
        w_in <- max(5, (if (!is.null(input$plot_w)) input$plot_w else 700) / 96)
        ggplot2::ggsave(
          filename  = file,
          plot      = p,
          width     = w_in,
          height    = h_in,
          units     = "in",
          bg        = if (!is.null(input$bg_col)) input$bg_col else "#FFFFFF",
          limitsize = FALSE
        )
      }
    )

  })
}
