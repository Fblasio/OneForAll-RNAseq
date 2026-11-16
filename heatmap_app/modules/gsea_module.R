# =============================================================================
# gsea_module.R  (goes in modules/, sourced from global.R after go_module.R)
# -----------------------------------------------------------------------------
# GSEA panel: rank-based enrichment on the FULL differential-expression result.
#
# Design mirrors go_module.r: a controls UI (left column) + an output UI
# (tabs: summary table, NES plot, and a classic running-score plot). All the
# heavy lifting lives in utils/gsea_functions.R; this module only wires inputs
# to run_gsea() and renders the outputs.
#
# IMPORTANT — what to feed it:
#   gseaServer("gsea_de", de_results_r = rnaseq_de$de_results_all)
#   Pass de_results_ALL (every gene), NOT de_sig_genes. GSEA ranks the whole
#   list; handing it only the significant genes defeats the point.
#
# Wiring recap (see app.r instructions):
#   UI (main panel):  gseaControlsUI("gsea_de")  +  gseaOutputUI("gsea_de")
#   Server:           gseaServer("gsea_de", de_results_r = rnaseq_de$de_results_all)
# =============================================================================


gseaControlsUI <- function(id) {
  ns <- NS(id)
  tagList(
    h5("GSEA settings"),
    selectInput(ns("database"), "Gene sets",
                choices = c("GO" = "GO", "KEGG pathways" = "KEGG"),
                selected = "GO"),
    conditionalPanel(
      condition = paste0("input['", ns("database"), "'] == 'GO'"),
      selectInput(ns("ontology"), "GO ontology",
                  choices = c("BP", "MF", "CC", "ALL"), selected = "BP")
    ),
    selectInput(ns("rank_metric"), "Rank genes by",
                choices = c(
                  "Wald statistic (recommended)" = "stat",
                  "log2 fold change (shrunken)"  = "log2fc",
                  "sign(log2FC) x -log10(p)"     = "signed_p"
                ),
                selected = "stat"),
    helpText(HTML(paste0(
      "<b>Wald statistic</b> = log2FC / SE. Signed, effect-size aware, single model. ",
      "Available whether or not shrinkage is on.<br>",
      "<b>log2 fold change</b> on a shrunken result ranks by effect size; ashr makes ",
      "log2FC comparable across genes of different precision.<br>",
      "<b>signed -log10(p)</b> ranks by significance only: a highly expressed gene ",
      "changing 5% outranks a lowly expressed gene changing 16-fold. Kept for ",
      "compatibility."
    ))),
    selectInput(ns("padj_method"), "p.adjust method",
                choices = c("BH", "BY", "bonferroni", "holm",
                            "hochberg", "hommel", "none"),
                selected = "BH"),
    numericInput(ns("pval"), "p.adjust cutoff", value = 0.05,
                 min = 0, max = 1, step = 0.01),
    fluidRow(
      column(6, numericInput(ns("min_size"),
                             lbl("Min set size", "Gene sets with fewer genes than this are not tested: a 3-gene set gives an unstable score and still costs you a hypothesis in the correction. 10 is standard."),
                             value = 10, min = 1, step = 1)),
      column(6, numericInput(ns("max_size"),
                             lbl("Max set size", "Sets larger than this are skipped - very broad sets are almost always mildly enriched and crowd the top of the table. 500 is standard. Together with Min set size this changes every adjusted p-value."),
                             value = 500, min = 10, step = 10))
    ),
    checkboxInput(ns("use_cache"), "Use cache", value = TRUE),
    actionButton(ns("run"), "\u25b6 Run GSEA",
                 class = "btn-info btn-sm", width = "100%"),
    br(), br(),
    uiOutput(ns("rank_info"))
  )
}


gseaOutputUI <- function(id) {
  ns <- NS(id)
  tabsetPanel(
    tabPanel(
      "GSEA table",
      br(),
      downloadButton(ns("dl_table"), "\u2b07 Download CSV",
                     class = "btn-sm btn-success"),
      br(), br(),
      helpText("NES > 0: set enriched at the top of the ranking (up); ",
               "NES < 0: enriched at the bottom (down)."),
      DTOutput(ns("table"))
    ),
    tabPanel(
      "NES plot",
      br(),
      fluidRow(
        column(4, sliderInput(ns("top_n"), "Top N sets (by |NES|)",
                              5, 60, 20, 5, ticks = FALSE)),
        column(4, selectInput(ns("plot_style"), "Style",
                              choices = c("Bars" = "bar", "Dots" = "dot"),
                              selected = "bar")),
        column(4, sliderInput(ns("plot_fs"), "Font size",
                              8, 24, 13, 1, ticks = FALSE))
      ),

      # ── Personalization panel ──────────────────────────────────────────
      # Native HTML <details>/<summary> instead of a JS accordion library:
      # closed by default (no visual clutter until the user opts in), one
      # click to open, zero new dependencies, and immune to whatever Shiny
      # theme/bootstrap version the app is running. Built via htmltools::tag()
      # rather than tags$details()/tags$summary() so it works regardless of
      # htmltools version (tag() is the generic constructor, guaranteed to
      # exist; the tags$xxx() shorthands for newer HTML5 elements aren't in
      # every version).
      htmltools::tag("details", list(
        style = "background:#f8f9fa; border:1px solid #dee2e6; border-radius:4px;
                 padding:6px 10px 10px 10px; margin:10px 0;",
        htmltools::tag("summary", list(
          style = "cursor:pointer; font-weight:bold; color:#1F3864; padding:4px 0;",
          "\U0001f3a8 Personalize plot style"
        )),
        div(
          style = "margin-top:8px;",
          fluidRow(
            column(4, selectInput(ns("style_font_family"), "Font",
                                  choices = c("Sans-serif" = "sans",
                                              "Serif" = "serif",
                                              "Monospace" = "mono"),
                                  selected = "sans")),
            column(4, selectInput(ns("style_theme"), "Theme",
                                  choices = c("Black & white" = "bw",
                                              "Minimal" = "minimal",
                                              "Classic" = "classic"),
                                  selected = "bw")),
            column(4, numericInput(ns("style_label_wrap"), "Wrap labels at (chars)",
                                   value = 45, min = 0, max = 120, step = 5))
          ),
          fluidRow(
            column(6, textInput(ns("style_title"), "Custom title",
                                placeholder = "(default: no title)")),
            column(6, textInput(ns("style_xlab"), "Custom x-axis label",
                                placeholder = "(default: Normalized Enrichment Score)"))
          ),
          fluidRow(
            column(4, selectInput(ns("style_color_by"), "Color by",
                                  choices = c("Significance (gradient)" = "significance",
                                              "Direction (up vs down)"  = "direction",
                                              "Single color"            = "flat"),
                                  selected = "significance"))
          ),
          conditionalPanel(
            condition = paste0("input['", ns("style_color_by"), "'] == 'significance'"),
            fluidRow(
              column(4, selectInput(ns("style_color_metric"), "Color scale by",
                                    choices = c("p.adjust" = "p.adjust",
                                                "-log10(p.adjust)" = "neglog10_padj"),
                                    selected = "p.adjust")),
              column(4, colourpicker::colourInput(ns("style_color_low"),  "Low",  value = "#3B7DD8")),
              column(4, colourpicker::colourInput(ns("style_color_high"), "High", value = "#C0392B"))
            )
          ),
          conditionalPanel(
            condition = paste0("input['", ns("style_color_by"), "'] == 'direction'"),
            fluidRow(
              column(6, colourpicker::colourInput(ns("style_up_color"),   "Up (NES > 0)",   value = "#C0392B")),
              column(6, colourpicker::colourInput(ns("style_down_color"), "Down (NES < 0)", value = "#3B7DD8"))
            )
          ),
          conditionalPanel(
            condition = paste0("input['", ns("style_color_by"), "'] == 'flat'"),
            fluidRow(
              column(4, colourpicker::colourInput(ns("style_flat_color"), "Bar / point color", value = "#3B7DD8"))
            )
          )
        )
      )),

      plotOutput(ns("nes_plot"), height = "600px"),
      br(),
      downloadButton(ns("dl_nes_plot"), "\u2b07 Download plot (PNG)", class = "btn-sm btn-success")
    ),
    tabPanel(
      "Running score",
      br(),
      selectizeInput(ns("gs_pick"), "Gene set to display",
                     choices = NULL,
                     options = list(placeholder = "Run GSEA, then pick a set...")),

      # ── Personalization panel — same closed-by-default <details> pattern
      # as the NES plot tab above, kept internally consistent within GSEA.
      htmltools::tag("details", list(
        style = "background:#f8f9fa; border:1px solid #dee2e6; border-radius:4px;
                 padding:6px 10px 10px 10px; margin:10px 0;",
        htmltools::tag("summary", list(
          style = "cursor:pointer; font-weight:bold; color:#1F3864; padding:4px 0;",
          "\U0001f3a8 Personalize plot style"
        )),
        div(
          style = "margin-top:8px;",
          fluidRow(
            column(4, colourpicker::colourInput(ns("run_color"), "Running score line color", value = "#2CA02C")),
            column(4, selectInput(ns("run_es_geom"), "Running score geometry",
                                  choices = c("Line" = "line", "Dots" = "dot"),
                                  selected = "line")),
            column(4, numericInput(ns("run_base_size"), "Font size", value = 11, min = 6, max = 24, step = 1))
          ),
          fluidRow(
            column(6, checkboxGroupInput(ns("run_subplots"), "Panels to show",
                                         choices = c("Running score" = "1",
                                                     "Hits (gene positions)" = "2",
                                                     "Ranked metric" = "3"),
                                         selected = c("1", "2", "3"), inline = TRUE)),
            column(6, checkboxInput(ns("run_pvalue_table"), "Show NES / p-value table", value = FALSE))
          )
        )
      )),

      plotOutput(ns("running_plot"), height = "500px"),
      br(),
      downloadButton(ns("dl_running_plot"), "\u2b07 Download plot (PNG)", class = "btn-sm btn-success")
    )
  )
}


gseaServer <- function(id, de_results_r) {
  moduleServer(id, function(input, output, session) {

    # ---- ranked vector, rebuilt whenever DE results or the metric change ----
    rank_vec <- reactive({
      de <- tryCatch(de_results_r(), error = function(e) NULL)
      build_gsea_rank(de, metric = input$rank_metric %||% "stat")
    })

    output$rank_info <- renderUI({
      rv <- rank_vec()
      if (is.null(rv))
        return(helpText("No ranked list yet \u2014 run a DE contrast first."))
      helpText(sprintf("Ranked list: %d genes. Top: %s (%.2f)  Bottom: %s (%.2f)",
                       length(rv), names(rv)[1], rv[1],
                       names(rv)[length(rv)], rv[length(rv)]))
    })

    # ---- run GSEA on button press ------------------------------------------
    result <- eventReactive(input$run, {
      rv <- rank_vec()
      # NOTE: must be shiny::validate()/shiny::need() here, fully qualified.
      # This app's global.R does library(jsonlite) *after* library(shiny),
      # and jsonlite also exports a function called validate() (for
      # validating JSON strings; its argument is literally named `txt`,
      # which is exactly what shows up in the "is.character(txt) is not
      # TRUE" crash if this is left as a bare, unqualified validate() call).
      # jsonlite::validate() then masks shiny::validate() on the search
      # path, so an unqualified call silently resolves to the wrong
      # function and crashes with that exact error the moment this runs.
      shiny::validate(shiny::need(!is.null(rv),
                    "No ranked gene list. Run a differential-expression contrast first."))
      # with_progress_messages() rather than withProgress(): fgsea's multilevel
      # step with eps = 0 is the slowest thing in this panel, and it narrates
      # its own progress. Relaying that is the difference between "working" and
      # "frozen". See log_utils.R.
      with_progress_messages({
        run_gsea(
          rank_vec    = rv,
          database    = input$database %||% "GO",
          ontology    = input$ontology %||% "BP",
          padj_method = input$padj_method %||% "BH",
          pval_cutoff = input$pval %||% 0.05,
          min_size    = input$min_size %||% 10,
          max_size    = input$max_size %||% 500,
          rank_metric = input$rank_metric %||% "signed_p",
          use_cache   = isTRUE(input$use_cache)
        )
      }, "Running GSEA (fgsea)")
    })

    # ---- surface run warnings -----------------------------------------------
    # Mirrors go_module.r:234-239 and kegg_module.R:90-99. Without this the GSEA
    # panel had no warning channel at all, so anything run_gsea() reported was
    # discarded before reaching the user or the session log.
    observeEvent(result(), {
      for (w in result()$warnings %||% character(0)) {
        showNotification(paste0("⚠️ ", w), type = "warning", duration = 14)
        log_warning(w)
      }
    })

    # ---- summary table ------------------------------------------------------
    output$table <- renderDT({
      res <- result()
      req(!is.null(res))
      tbl <- res$table

      # empty_result() payload: single Message column -> show the reason.
      if (identical(colnames(tbl), "Message")) {
        return(datatable(tbl, rownames = FALSE, options = list(dom = "t")))
      }

      show_cols <- intersect(
        c("ID", "Description", "setSize", "NES", "enrichmentScore",
          "pvalue", "p.adjust", "qvalue"),
        colnames(tbl)
      )
      disp <- tbl[, show_cols, drop = FALSE]
      num  <- intersect(c("NES", "enrichmentScore", "pvalue", "p.adjust", "qvalue"),
                        colnames(disp))
      for (cc in num) disp[[cc]] <- signif(as.numeric(disp[[cc]]), 3)

      datatable(disp, rownames = FALSE, selection = "single",
                filter = "top",
                options = list(scrollX = TRUE, pageLength = 15))
    }, server = FALSE)

    output$dl_table <- downloadHandler(
      filename = function()
        paste0("GSEA_", input$database %||% "GO", "_",
               format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
      content = function(file) {
        res <- result()
        # Same '#'-comment header as the GO/KEGG exports, so a GSEA table that
        # dropped part of its ranking to KEGG ID mapping says so in the file.
        if (is.null(res)) utils::write.csv(data.frame(), file, row.names = FALSE)
        else              write_result_csv(res, file)
      }
    )

    # ---- NES summary plot ---------------------------------------------------
    build_nes_plot <- function() {
      res <- result()
      req(!is.null(res))
      make_gsea_summary_plot(
        res$table,
        top_n        = input$top_n %||% 20,
        font_size    = input$plot_fs %||% 13,
        style        = input$plot_style %||% "bar",
        # ── personalization panel (collapsed by default; these all have
        # sane defaults so the plot looks identical until the user opens
        # the panel and changes something) ──
        font_family  = input$style_font_family %||% "sans",
        theme_style  = input$style_theme %||% "bw",
        label_wrap   = input$style_label_wrap %||% 45,
        title        = input$style_title %||% "",
        x_lab        = input$style_xlab %||% "",
        color_by     = input$style_color_by %||% "significance",
        color_metric = input$style_color_metric %||% "p.adjust",
        color_low    = input$style_color_low %||% "#3B7DD8",
        color_high   = input$style_color_high %||% "#C0392B",
        up_color     = input$style_up_color %||% "#C0392B",
        down_color   = input$style_down_color %||% "#3B7DD8",
        flat_color   = input$style_flat_color %||% "#3B7DD8"
      )
    }

    output$nes_plot <- renderPlot({ build_nes_plot() })

    output$dl_nes_plot <- downloadHandler(
      filename = function()
        paste0("GSEA_NES_plot_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
      content = function(file) {
        ggplot2::ggsave(file, plot = build_nes_plot(), width = 9, height = 7, dpi = 300)
      }
    )

    # ---- populate the running-score picker from the result ------------------
    observeEvent(result(), {
      res <- result()
      tbl <- res$table
      if (is.null(tbl) || identical(colnames(tbl), "Message") ||
          !all(c("ID", "Description") %in% colnames(tbl))) {
        updateSelectizeInput(session, "gs_pick", choices = character(0),
                             server = TRUE)
        return()
      }
      choices <- setNames(tbl$ID, paste0(tbl$Description, "  (", tbl$ID, ")"))
      updateSelectizeInput(session, "gs_pick", choices = choices,
                           selected = tbl$ID[1], server = TRUE)
    })

    # ---- classic running-score plot (enrichplot::gseaplot2) -----------------
    build_running_plot <- function() {
      res <- result()
      req(!is.null(res), !is.null(res$object))
      gs <- input$gs_pick
      req(gs, gs %in% res$object@result$ID)

      subplots_sel <- as.integer(input$run_subplots %||% c("1", "2", "3"))
      # Panel 1 (the running-score line) is the base the other panels attach
      # to in gseaplot2's internal layout — always include it even if the
      # user unchecked it, rather than letting an unusual combination throw.
      subplots_sel <- sort(unique(c(1L, subplots_sel)))

      p <- tryCatch(
        enrichplot::gseaplot2(
          res$object, geneSetID = gs,
          title = res$object@result$Description[match(gs, res$object@result$ID)],
          color        = input$run_color %||% "#2CA02C",
          base_size    = input$run_base_size %||% 11,
          subplots     = subplots_sel,
          pvalue_table = isTRUE(input$run_pvalue_table),
          ES_geom      = input$run_es_geom %||% "line"
        ),
        error = function(e) e
      )
      shiny::validate(shiny::need(!inherits(p, "error"),
                    paste("Could not draw the running-score plot with this panel combination:",
                          if (inherits(p, "error")) conditionMessage(p) else "")))
      p
    }

    output$running_plot <- renderPlot({ build_running_plot() })

    output$dl_running_plot <- downloadHandler(
      filename = function()
        paste0("GSEA_running_score_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
      content = function(file) {
        ggplot2::ggsave(file, plot = build_running_plot(), width = 8, height = 7, dpi = 300)
      }
    )

    # ── Accumulate EVERY run this session (for full-history reproducibility) ──
    # Mirrors de_comparison_log.R / the GO module: each successful GSEA run is
    # appended with BOTH its params and the exact ranked vector used, so the
    # reproducibility bundle can regenerate all runs, not just the current one.
    gsea_history <- reactiveVal(list())
    # NO ignoreInit, deliberately -- same reasoning as go_module.r:620-626.
    # result() is an eventReactive that ABORTS at init (input$run == 0 is a null
    # event), so ignoreInit's skip was consumed by the user's FIRST Run GSEA
    # click rather than by init: run 1 was never recorded and the bundle was
    # permanently one GSEA run short. The is.null(res$params) guard below
    # already makes an init-time fire a no-op.
    observeEvent(result(), {
      res <- result()
      if (is.null(res) || is.null(res$params)) return()
      rv <- tryCatch(rank_vec(), error = function(e) NULL)
      gsea_history(c(gsea_history(), list(list(
        params = res$params,
        ranked_ids    = if (is.null(rv)) character(0) else names(rv),
        ranked_scores = if (is.null(rv)) numeric(0)   else unname(rv)
      ))))
    })

    # Expose params for the reproducibility bundle (optional wiring in app.r).
    # `ranked` lets the repro module save the exact ranked vector to disk so
    # rerun_GSEA.R can rebuild it (see gene_lists/gsea_1_ranked_* in the bundle).
    # `history` is the accumulating list of every run (params + ranked vector).
    list(
      result     = reactive(result()),
      ranked     = reactive(rank_vec()),
      history    = reactive(gsea_history()),
      gsea_params = reactive({
        res <- tryCatch(result(), error = function(e) NULL)
        if (is.null(res)) NULL else res$params
      })
    )
  })
}
