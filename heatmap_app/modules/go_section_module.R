# =============================================================================
# go_section_module.R
# -----------------------------------------------------------------------------
# Standalone "GO Analysis" feature panel — replaces the old placeholder
# "Comparison" section in the sidebar/main-panel feature switcher.
#
# DESIGN GOAL: this module owns nothing about GO enrichment itself. It only
# supplies a gene list (+ optional universe) and then hands off to the
# existing, unmodified go_module.r building blocks:
#   goControlsUI() / goTableUI() / goPlotUI() / goHeatmapUI() / goServer()
# via a nested sub-module ("enrich"), the same way app.r already wires up
# the "go" (heatmap-brush) and "go_de" (DEG) panels. That means:
#   - Any future change to GO enrichment logic (enrichment_functions.R) or
#     to the plot/table/heatmap UI (go_module.r) is automatically picked up
#     here too, with zero duplication.
#   - This file adds ONLY the parts that are genuinely specific to a
#     "paste your own gene list" entry point: the gene-list box, an optional
#     custom background/universe box, and two convenience buttons to pull
#     genes from other parts of the app (current heatmap brush selection,
#     current DE-significant gene list) when those are available.
#
# Wiring (see app.r):
#   UI:     goSectionSidebarUI("go_section")   -> sidebar panel
#           goSectionMainUI("go_section")      -> main panel (tabs)
#   Server: goSectionServer("go_section",
#                            scaled_matrix_rv = scaled_matrix_rv,   # optional
#                            gene_labels_rv   = gene_labels_rv,     # optional
#                            col_fun_rv       = col_fun_rv,         # optional
#                            col_order_rv     = col_order_rv,       # optional
#                            deg_genes_rv     = deg_genes_rv,       # optional
#                            heatmap_selected_genes = selected_genes,   # optional bridge
#                            de_sig_genes_rv        = rnaseq_de$de_sig_genes)  # optional bridge
#
# All server arguments are optional (default NULL) and degrade gracefully:
#   - Without scaled_matrix_rv/gene_labels_rv/col_fun_rv/col_order_rv, the
#     "GO Heatmap" sub-tab simply shows its existing empty-state message
#     (go_module.r already handles a NULL matrix, no changes needed there).
#   - Without heatmap_selected_genes / de_sig_genes_rv, the two bridge
#     buttons stay visible but show a notification explaining nothing is
#     available yet, instead of erroring or being silently removed.
# =============================================================================

# -----------------------------------------------------------------------------
# UI — sidebar: gene list input, optional universe, then the reused GO controls
# -----------------------------------------------------------------------------
goSectionSidebarUI <- function(id) {
  ns <- NS(id)
  tagList(
    h4("\U0001f9ec Gene Ontology Analysis"),
    helpText("Paste a gene list (or pull one from the Heatmap / RNAseq DE ",
             "panels below) and run GO enrichment on it directly \u2014 no ",
             "heatmap required."),

    textAreaInput(ns("gene_list_text"),
                  label       = "Gene IDs (one per line)",
                  placeholder = GENE_ID_PLACEHOLDER,
                  rows        = 6),
    fluidRow(
      column(6, actionButton(ns("load_genes"), "\U0001f4e5 Load gene list",
                             class = "btn-primary btn-sm", width = "100%")),
      column(6, actionButton(ns("clear_genes"), "\U0001f9f9 Clear",
                             class = "btn-default btn-sm", width = "100%"))
    ),
    uiOutput(ns("gene_count_info")),

    br(),
    tags$b("Pull genes from elsewhere in the app:"),
    fluidRow(
      column(6, actionButton(ns("use_heatmap_sel"), "\U0001f525 Heatmap selection",
                             class = "btn-info btn-sm", width = "100%")),
      column(6, actionButton(ns("use_de_sig"), "\U0001f9ec DE-significant genes",
                             class = "btn-info btn-sm", width = "100%"))
    ),
    hr(),

    h5("\U0001f310 Gene Universe (background)"),
    radioButtons(ns("universe_mode"), "Background gene set for GO enrichment",
                 choices  = c("Default (whole-genome background)" = "default",
                              "Custom pasted gene list"            = "custom"),
                 selected = "default"),
    conditionalPanel(
      condition = paste0("input['", ns("universe_mode"), "'] == 'custom'"),
      textAreaInput(ns("universe_text"),
                    label       = "Universe gene IDs (one per line)",
                    placeholder = GENE_ID_PLACEHOLDER,
                    rows        = 5),
      helpText("Accepts newlines, spaces or commas."),
      uiOutput(ns("universe_info"))
    ),
    hr(),

    goControlsUI(ns("enrich"))
  )
}

# -----------------------------------------------------------------------------
# UI — main panel: current-selection header + the reused GO table/plot/heatmap
# -----------------------------------------------------------------------------
goSectionMainUI <- function(id) {
  ns <- NS(id)
  tagList(
    uiOutput(ns("header_info")),
    tabsetPanel(
      tabPanel("GO Table",   goTableUI(ns("enrich"))),
      tabPanel("GO Plot",    goPlotUI(ns("enrich")))
      # "GO Heatmap" is deliberately NOT included here. It shows the expression
      # of a clicked term's genes using scaled_matrix_rv, which in this
      # gene-list-first section is either NULL (no heatmap generated -> tab is
      # just empty) or, worse, belongs to an unrelated earlier heatmap session
      # -> genuinely misleading. goServer()/go_module.r are untouched, so the
      # heatmap-brush ("go") and DEG ("go_de") panels still show it normally,
      # where scaled_matrix_rv is always the data the gene list actually came
      # from. See the chat writeup for the "Term Genes" tab alternative if a
      # matrix-free per-term gene view is wanted here later.
    )
  )
}

# -----------------------------------------------------------------------------
# SERVER
# -----------------------------------------------------------------------------
goSectionServer <- function(id,
                            scaled_matrix_rv       = NULL,
                            # Pass-through only; consumed by goServer() below.
                            cluster_matrix_rv      = NULL,
                            linkage_rv             = NULL,
                            gene_labels_rv         = NULL,
                            col_fun_rv             = NULL,
                            col_order_rv           = NULL,
                            deg_genes_rv           = NULL,
                            heatmap_selected_genes = NULL,
                            de_sig_genes_rv        = NULL,
                            # Optional, NULL-defaulted like every other argument
                            # here: without it the module behaves exactly as
                            # before, it just cannot refresh its placeholders.
                            species_rv             = NULL) {
  moduleServer(id, function(input, output, session) {

    # ── Keep the ID-format hints in step with the active species ────────────
    # goSectionSidebarUI() sits inside a conditionalPanel in app.r's static
    # sidebar, so GENE_ID_PLACEHOLDER is read once at startup. Without this
    # observer the panel would keep showing the previous organism's ID format
    # after a species switch -- in the one panel whose whole purpose is pasting
    # gene IDs by hand.
    # `value` is left NULL so a list the user has already pasted is not cleared
    # out from under them.
    if (!is.null(species_rv)) {
      observeEvent(species_rv(), {
        updateTextAreaInput(session, "gene_list_text", placeholder = GENE_ID_PLACEHOLDER)
        updateTextAreaInput(session, "universe_text",  placeholder = GENE_ID_PLACEHOLDER)
      }, ignoreNULL = FALSE)
    }

    # ── Gene list: parsed LIVE from the textarea ────────────────────────────
    # Deliberately NOT a reactiveVal gated behind a "Load" button: Run GO
    # (inside the nested go_module.r controls) reads this reactive at the
    # moment it's clicked, so whatever is currently in the box is what gets
    # used — no separate "did you remember to click Load first?" step that
    # can silently leave Run GO with nothing to run on.
    selected_genes <- reactive({
      ids <- parse_id_textarea(clean_pasted_ids(input$gene_list_text))
      if (length(ids) == 0) return(NULL)
      ids
    })

    # "Load gene list" is now just an explicit confirmation/validation step
    # for the user (shows how many IDs were recognized) — Run GO does not
    # depend on this having been clicked.
    observeEvent(input$load_genes, {
      genes <- selected_genes()
      if (is.null(genes)) {
        showNotification("\u274c No valid gene IDs found in the text box.",
                          type = "error", duration = 6)
      } else {
        showNotification(paste0("\u2705 ", length(genes), " gene(s) recognized. ",
                                 "You can now click Run GO."),
                          type = "message", duration = 4)
      }
    })

    observeEvent(input$clear_genes, {
      updateTextAreaInput(session, "gene_list_text", value = "")
    })

    # ── Bridge buttons: pull genes from other parts of the app ─────────────
    # Both just populate the textarea; selected_genes() picks the change up
    # automatically via the reactive above.
    observeEvent(input$use_heatmap_sel, {
      genes <- if (is.function(heatmap_selected_genes)) heatmap_selected_genes() else NULL
      if (is.null(genes) || length(genes) == 0) {
        showNotification(
          "\u2139\ufe0f No heatmap gene selection available yet \u2014 generate a heatmap and brush a cluster first.",
          type = "warning", duration = 6)
        return()
      }
      genes <- canon_ids(genes)
      updateTextAreaInput(session, "gene_list_text", value = paste(genes, collapse = "\n"))
      showNotification(paste0("\u2705 ", length(genes), " gene(s) pulled from the heatmap selection. ",
                               "Click Run GO when ready."),
                        type = "message", duration = 5)
    })

    observeEvent(input$use_de_sig, {
      genes <- if (is.function(de_sig_genes_rv)) tryCatch(de_sig_genes_rv(), error = function(e) NULL) else NULL
      if (is.null(genes) || length(genes) == 0) {
        showNotification(
          "\u2139\ufe0f No DE-significant gene list available yet \u2014 run a differential expression analysis first.",
          type = "warning", duration = 6)
        return()
      }
      genes <- canon_ids(genes)
      updateTextAreaInput(session, "gene_list_text", value = paste(genes, collapse = "\n"))
      showNotification(paste0("\u2705 ", length(genes), " DE-significant gene(s) pulled in. ",
                               "Click Run GO when ready."),
                        type = "message", duration = 5)
    })

    output$gene_count_info <- renderUI({
      genes <- selected_genes()
      if (is.null(genes) || length(genes) == 0) {
        return(tags$p(style = "color:gray; font-size:0.85em; font-style:italic;",
                      "No gene list loaded yet."))
      }
      tags$p(style = "color:#2c7bb6; font-weight:bold; font-size:0.9em;",
             icon("circle-info"), paste0(" ", length(genes), " gene(s) loaded."))
    })

    # ── Universe (background) ───────────────────────────────────────────────
    universe_genes <- reactive({
      if (!identical(input$universe_mode, "custom")) return(NULL)
      ids <- parse_id_textarea(clean_pasted_ids(input$universe_text))
      if (length(ids) == 0) return(NULL)
      ids
    })

    output$universe_info <- renderUI({
      u <- universe_genes()
      if (is.null(u)) {
        return(tags$p(style = "color:#999; font-size:0.85em;",
                      "Paste at least one gene ID above."))
      }
      tags$p(style = "color:#2c7bb6; font-weight:bold; font-size:0.85em;",
             icon("circle-info"), paste0(" ", length(u), " genes in custom universe"))
    })

    # ── Header shown above the Table/Plot/Heatmap tabs ─────────────────────
    output$header_info <- renderUI({
      genes <- selected_genes()
      n     <- if (is.null(genes)) 0 else length(genes)
      tagList(
        h4(paste0("\U0001f52c GO Enrichment \u2014 ", n, " gene(s) in list")),
        if (!identical(input$universe_mode %||% "default", "default") && !is.null(universe_genes())) {
          tags$p(style = "color:#2c7bb6; font-size:0.9em;",
                 icon("circle-info"),
                 paste0(" Custom universe active: ", length(universe_genes()), " genes"))
        }
      )
    })

    # ── Hand off to the existing, unmodified GO module ─────────────────────
    go_handles <- goServer(
      "enrich",
      selected_genes   = selected_genes,
      scaled_matrix_rv = scaled_matrix_rv,
      cluster_matrix_rv = cluster_matrix_rv,
      linkage_rv       = linkage_rv,
      gene_labels_rv   = gene_labels_rv,
      universe         = universe_genes,
      col_fun_rv       = col_fun_rv,
      col_order_rv     = col_order_rv,
      deg_genes_rv     = deg_genes_rv,
      context          = "go_section"
    )

    list(
      result         = go_handles$result,
      selected_genes = selected_genes
    )
  })
}
