options(shiny.maxRequestSize = 500 * 1024^2)
source("global.R")

# clean_pasted_ids() moved to consistency_utils.R (loaded via global.R) so
# it's visible to module files as well as this script's own server function
# — see consistency_utils.R for the function and the reason it had to move.

ui <- fluidPage(
  shinyjs::useShinyjs(),
  # Without this the browser tab and any bookmark show the bare URL.
  tags$head(
    tags$title("One for All"),
    # Re-measure plots when a collapsible_panel() is expanded.
    #
    # A closed <details> is display:none, so anything inside reports a width of
    # 0. Shiny sizes renderPlot() from that client-reported width, so a plot
    # whose FIRST render happens while the panel is folded comes back blank or
    # at the wrong size, and stays that way after expanding -- nothing tells it
    # to look again. Reachable in one sequence: fold "Exploratory Data
    # Analysis", then build the DESeq2 object, then unfold.
    #
    # Firing a resize on expand makes every output re-report its size. Cheap,
    # and it covers every current and future collapsible_panel() that holds a
    # plot, not just the EDA block.
    #
    # Capture phase (the `true`): the `toggle` event does NOT bubble, so a
    # listener on document only sees it on the way down.
    tags$script(HTML(
      "document.addEventListener('toggle', function(e){",
      "  if (e.target && e.target.tagName === 'DETAILS' && e.target.open) {",
      "    window.dispatchEvent(new Event('resize'));",
      "  }",
      "}, true);"
    )),
    # Styling for help_tip() / lbl() (global.R). Pure CSS so no package is
    # added; the text lives in data-tip and is revealed on hover or keyboard
    # focus. The bubble is placed BELOW the icon rather than beside it because
    # the sidebar is only 3 columns wide, and it is given a high z-index so it
    # can overflow the well panel instead of being clipped by it.
    tags$style(HTML("
.cheski-help {
  display:inline-flex; align-items:center; justify-content:center;
  width:15px; height:15px; margin-left:5px; position:relative;
  border-radius:50%; border:1px solid #7a8ba8;
  background:#eef2f8; color:#33507a;
  font-size:10px; font-weight:700; line-height:1;
  cursor:help; vertical-align:middle; user-select:none;
}
.cheski-help::after {
  content: attr(data-tip);
  position:absolute; top:20px; left:-6px;
  width:300px; max-width:70vw; padding:8px 10px;
  background:#1F3864; color:#ffffff; border-radius:4px;
  font-size:12px; font-weight:400; font-style:normal; line-height:1.45;
  white-space:normal; text-align:left;
  opacity:0; visibility:hidden; transition:opacity .12s ease;
  z-index:9999; pointer-events:none;
  box-shadow:0 2px 10px rgba(0,0,0,.28);
}
.cheski-help:hover::after, .cheski-help:focus::after { opacity:1; visibility:visible; }

/* Keeps 'Generate Heatmap' reachable without scrolling the whole sidebar.
   The negative side margins bleed the strip out to the edges of the
   surrounding .well so it reads as a footer bar rather than a floating
   button; the padding puts the well's own inset back. position:sticky is
   inert when the sidebar fits on screen, and needs no ancestor to scroll -
   but it IS cancelled by an ancestor with overflow:hidden/auto, so keep the
   sidebar free of those. */
.cheski-sticky-run {
  position: sticky; bottom: 0; z-index: 20;
  margin: 0 -19px -19px; padding: 12px 19px;
  background: #f5f5f5; border-top: 1px solid #ddd;
  box-shadow: 0 -5px 12px -8px rgba(0,0,0,.35);
}
"))
  ),
  div(
    h2(HTML("\U0001f525 One for All \U0001f525")),
    h3(HTML("Clustering, Heatmaps, Enrichment and Systems Knowledge Interface")),
    style = "line-height: 1.2; text-align: center;"
  ),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      # Feature selector is displayed above all option panels.
      # Each button updates the active feature used by conditional panels.
      div(
        style = "margin-bottom: 12px;",
        h4("🔬 Select Feature", style = "margin-bottom: 8px;"),
        div(
          style = "display: flex; gap: 6px; margin-bottom: 6px;",
          actionButton("feat_heatmap", "🔥 Heatmap", class = "btn btn-default btn-sm feature-btn active-feature", width = "100%"),
          actionButton("feat_rnaseq", "🧬 RNAseq Analysis", class = "btn btn-default btn-sm feature-btn", width = "100%"),
          actionButton("feat_go", "🧬 GO Analysis", class = "btn btn-default btn-sm feature-btn", width = "100%"),
        ),
        
        # Row 2: Venn diagram (second line, full width)
        div(
          style = "display: flex; gap: 6px;",
          actionButton("feat_venn", "🟢 Venn diagram", class = "btn btn-default btn-sm feature-btn", width = "100%")
        ),
        
        tags$style(HTML("
.feature-btn { border: 2px solid #ccc; background: #f8f9fa; }
.active-feature { border: 2px solid #337ab7 !important; background: #d9edf7 !important; font-weight: bold; }
"))
      ),
      # Active feature is tracked client-side and exposed as input$active_feature.
      tags$script(HTML("
var activeFeature = 'heatmap';
function setFeature(f) {
  activeFeature = f;
  Shiny.setInputValue('active_feature', f, {priority: 'event'});
  document.querySelectorAll('.feature-btn').forEach(function(b) {
    b.classList.remove('active-feature');
  });
  document.getElementById('feat_' + f).classList.add('active-feature');
}
$(document).on('click', '#feat_heatmap', function() { setFeature('heatmap'); });
$(document).on('click', '#feat_rnaseq', function() { setFeature('rnaseq'); });
$(document).on('click', '#feat_go', function() { setFeature('go'); });
$(document).on('click', '#feat_venn', function() { setFeature('venn'); });
")),
      conditionalPanel(
        condition = "input.active_feature !== 'venn'",
        wellPanel(
          style = "background:#f0f4ff; border:1px solid #c5d0e6; padding:10px;",
          # This well panel is shared by the Heatmap, RNAseq and GO features, but
          # only the first two are a numbered sequence -- so the badge appears
          # there and the GO sidebar keeps the plain heading it always had,
          # rather than opening with a lone "1" that leads nowhere.
          conditionalPanel(
            condition = "input.active_feature !== 'go'",
            step_header(1, "Choose species")
          ),
          conditionalPanel(
            condition = "input.active_feature === 'go'",
            h5("\U0001f331 Active Species")
          ),
          uiOutput("species_label_ui"),
          selectInput(
            "species_switch",
            label = NULL,
            choices = c("Arabidopsis thaliana" = "arabidopsis",
                        "Solanum lycopersicum" = "tomato"),
            selected = ACTIVE_SPECIES
          ),
          actionButton(
            "apply_species", "🔄 Switch species & reload",
            class = "btn-primary btn-sm", width = "100%"
          )
        )
      ),
      hr(),
      conditionalPanel(
        condition = "input.active_feature === 'heatmap' || input.active_feature === undefined || input.active_feature === null",
      step_header(2, "Upload data"),
      fileInput("file", "Upload CSV or Excel",
                accept = c(".csv", ".xlsx", ".xls")),
      helpText("Col 1 = gene IDs. All remaining columns = expression values."),
      actionButton("demo_btn", "\U0001f9ea Load Arabidopsis demo (synthetic)", class = "btn-default btn-sm"),
      checkboxInput("is_raw_counts",
                    label = lbl("These are raw counts", "Tick for raw, un-normalised counts - applies DESeq2 VST first. Left off on real counts, highly expressed genes dominate both the colours and the clustering. Untick for TPM / CPM / FPKM / VST."),
                    value = FALSE),
      conditionalPanel(
        condition = "input.is_raw_counts == true",
        helpText(
          "\U0001f9ec Raw counts detected. The matrix will be automatically",
          "normalised with DESeq2 VST (blind = TRUE) before heatmap rendering."
        )
      ),
      hr(),
      step_header(3, "Filter genes"),
      textAreaInput("gene_filter",
                    label = "Gene IDs to display (one per line)",
                    placeholder = "AT1G01010\nAT2G05510\nAT3G18780",
                    rows = 6),
      actionButton("clear_gene_filter", "\U0001f9f9 Clear",
                   class = "btn-default btn-sm", width = "100%"),
      helpText("Leave empty to show all. Accepts newlines, spaces or commas."),
      numericInput("top_var_n",
                   label = "…or keep only the top N most variable genes (0 = all)",
                   value = 0, min = 0, step = 100),
      helpText("Ranked by variance across conditions on the VST/expression scale, ",
               "never on the display transform. Applied after the gene ID filter."),
      hr(),
      # "Gene Name Mapping" lives in the Labels tab of the appearance panel,
      # because it decides what is written down the side of the figure. Its div
      # keeps the id "map_gene_names_block": the species observer calls
      # shinyjs::show/hide on that id, so renaming it would silently disable the
      # hide-when-no-symbols behaviour.
      # \u2500\u2500 Experiment setup: groups + annotation factors + assignments \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500
      # These three sections are the bulk of the sidebar -- for a 4-group,
      # 2-factor design they expand to roughly 24 controls (4 name+selectize
      # pairs, 2 name+levels pairs plus one colour picker per level, and
      # n_groups x n_annots assignment dropdowns). They are also pure SETUP:
      # once the experiment is described you stop touching them, while the
      # controls you do keep adjusting (filters, normalization, clustering,
      # Generate) sit below and are pushed off-screen by them.
      #
      # Wrapped in one addressable div so the whole block can be folded away
      # after a successful Generate, leaving a one-line summary. Nothing is
      # removed or disabled: hiding is purely visual, every input keeps its
      # value and stays readable by the server exactly as before.
      # The step number is carried in the LABEL STRING, not by step_header(),
      # because this is an actionLink whose label the server rebuilds on every
      # fold (see the setup_collapsed observer). Both places must spell it the
      # same way or the number disappears the first time the panel is toggled.
      actionLink("toggle_setup", "\u25be 4 \u00b7 \U0001f4e6 Experiment setup",
                 style = paste0("font-size:15px; font-weight:bold; color:#1F3864; ",
                                "text-decoration:none; display:inline-block; ",
                                "margin-bottom:4px;")),
      uiOutput("setup_summary"),
      div(
        id = "setup_block",
        step_header("4a", "Define groups"),
        helpText("Each group = one heatmap column (mean of replicates)."),
        uiOutput("group_ui"),
        fluidRow(
          column(6, actionButton("add_group", "\u2795 Add group",
                                 class = "btn-success btn-sm", width = "100%")),
          column(6, actionButton("remove_group", "\u2796 Remove last",
                                 class = "btn-warning btn-sm", width = "100%"))
        ),
        hr(),
        step_header("4b", "Annotation categories"),
        helpText("Define factors and colors for top bars."),
        uiOutput("annot_cat_ui"),
        fluidRow(
          column(6, actionButton("add_annot", "\u2795 Add factor",
                                 class = "btn-success btn-sm", width = "100%")),
          column(6, actionButton("remove_annot", "\u2796 Remove last",
                                 class = "btn-warning btn-sm", width = "100%"))
        ),
        hr(),
        step_header("4c", "Assign groups to levels"),
        helpText("Say which level of each factor every group belongs to. ",
                 "This is what the coloured annotation bars above the heatmap show."),
        uiOutput("annot_assign_ui")
      ),
      hr(),
      # Heatmap colours, z-limits, label and annotation-bar options and the
      # figure export moved to the "Appearance & export" panel under the
      # heatmap itself (main panel). They are judged by looking at the figure,
      # so they belong next to it -- and they were roughly 400px of sidebar
      # sitting between the analysis controls and the Generate button.
      step_header(5, "Display & clustering"),
      selectInput("normalization", "Normalization method",
                  choices  = normalization_choices,
                  selected = "zscore_row"),
      # Reference condition selector — visible only when Log2FC method is active.
      # The conditionalPanel is evaluated client-side (no server round-trip).
      conditionalPanel(
        condition = "input.normalization == 'log2fc'",
        uiOutput("ref_condition_ui"),
        helpText("Reference column(s) are hidden from the heatmap (all-zero by definition).")
      ),
      checkboxInput("cluster_rows",
                    lbl("Cluster rows (Pearson)", "Groups genes by the SHAPE of their profile (1 - Pearson), not by absolute level. Always computed on the VST matrix, never on the display transform, so the colour scale can never change the dendrogram."),
                    value = TRUE),
      conditionalPanel(
        condition = "input.cluster_rows == true",
        # Choices come from CLUSTER_LINKAGE_CHOICES (normalization_utils.R) so the
        # UI can never offer a method cluster_rows_or_reason() would reject.
        selectInput("linkage", "Linkage method",
                    choices  = CLUSTER_LINKAGE_CHOICES,
                    selected = "complete"),
        checkboxInput("cut_tree",
                      lbl("Cut dendrogram into k clusters", "Cuts the row dendrogram into k blocks, numbered 1-k from the top. Each gets a button that sends its genes to enrichment, plus a CSV export. Use the elbow and silhouette panels below to choose k."),
                      value = FALSE),
        conditionalPanel(
          condition = "input.cut_tree == true",
          sliderInput("k_clusters", "Number of clusters (k)",
                      min = 2, max = 10, value = 3, step = 1),
          optimalKUI("optimal_k")
        )
      ),
      checkboxInput("cluster_cols",
                    lbl("Cluster columns (Spearman)", "Reorders the conditions by profile similarity (1 - Spearman, complete linkage - independent of the row Linkage setting). Useful as a replicate check; switch off to keep a fixed order for a figure."),
                    value = FALSE),
      hr(),
      step_header(6, "Highlight DEGs", optional = TRUE),
      # Label reworded only to avoid repeating the step header verbatim one line
      # above it; the input id and behaviour are unchanged.
      checkboxInput("show_degs", "\U0001f9ec Highlight a list of DEGs", value = FALSE),
      conditionalPanel(
        condition = "input.show_degs == true",
        textAreaInput("deg_list",
                      label = "DEG gene IDs (one per line)",
                      placeholder = GENE_ID_PLACEHOLDER,
                      rows = 5),
        actionButton("clear_deg_list", "\U0001f9f9 Clear",
                     class = "btn-default btn-sm", width = "100%"),
        helpText("These genes will appear bold and magenta in all heatmaps.")
      ),
      hr(),
      hr(),
      step_header(7, "GO gene universe", optional = TRUE),
      radioButtons("universe_mode",
                   lbl("Background gene set for GO enrichment", "The gene list every GO and KEGG p-value is tested against. Default is the whole genome, including genes your tissue cannot express"),
                   choices = c("Default (whole-genome background)"        = "default",
                               "Expression threshold (derived from data)" = "threshold",
                               "Custom pasted gene list"                  = "custom"),
                   selected = "default"),
      conditionalPanel(
        condition = "input.universe_mode == 'threshold'",
        helpText("A gene enters the background when it reaches the threshold below in at ",
                 "least as many samples as there are replicates in your smallest group."),
        numericInput("universe_threshold",
                     label = lbl("Min. count threshold per sample", "A gene joins the background if it reaches this value in at least as many samples as your smallest group. Counted on the values as uploaded, before any transform. 10 is typical for raw counts."),
                     value = 10, min = 0, step = 1),
        uiOutput("universe_info")
      ),
      conditionalPanel(
        condition = "input.universe_mode == 'custom'",
        textAreaInput("universe_custom_genes",
                      label = "Universe gene IDs (one per line)",
                      placeholder = GENE_ID_PLACEHOLDER,
                      rows = 6),
        actionButton("clear_universe_custom", "\U0001f9f9 Clear",
                     class = "btn-default btn-sm", width = "100%"),
        helpText("Accepts newlines, spaces or commas. Used as the GO background instead of the whole genome."),
        uiOutput("universe_custom_info")
      ),
      hr(),
      # \u2500\u2500 Session tools (log, reproducibility bundle, cache) \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500
      # Folded by default: all three are per-session utilities, touched once if
      # at all, and unfolded they would occupy the stretch of sidebar directly
      # above the Generate button. "Export Figure" is deliberately not here --
      # it belongs with the appearance controls in the main panel, next to the
      # figure it exports.
      actionLink("toggle_session", "\u25b8 \U0001f9f0 Session tools",
                 style = paste0("font-size:15px; font-weight:bold; color:#1F3864; ",
                                "text-decoration:none; display:inline-block; ",
                                "margin-bottom:4px;")),
      shinyjs::hidden(div(
        id = "session_block",
        h4("\U0001f4dc Session Log"),
        helpText("Everything this session did, with its warnings, as a text file."),
        appLogUI("applog"),
        hr(),
        h4("\U0001f4c1 Reproducibility"),
        helpText("One .zip holding every heatmap, GO, KEGG, GSEA and DE run of ",
                 "this session, plus scripts that reproduce each of them."),
        reproDownloadUI("repro"),
        hr(),
        # The "?" goes on the heading, not on the button: putting it inside an
        # actionButton's label would make it part of the clickable area.
        h4(lbl("\U0001f5d1\ufe0f Cache", "GO, KEGG and GSEA results are cached per gene list + background + cutoffs, so a repeat run returns instantly. Clear only after switching species or replacing the annotation. Nothing is lost - the next run just recomputes.")),
        actionButton("clear_cache", "Clear enrichment cache",
                     class = "btn-sm btn-danger", width = "100%")
      )),
      hr(),
      # Pinned to the bottom of the viewport while the sidebar is scrolled, so
      # re-running never requires scrolling past the ~60 controls above it.
      # position:sticky achieves this with no JS and no layout change: when the
      # sidebar is shorter than the window the button renders in place as normal.
      div(class = "cheski-sticky-run",
          uiOutput("stale_badge"),
          # The number rides in the label rather than a step_header() above it:
          # this div is position:sticky, so anything added here permanently eats
          # vertical space at the bottom of the viewport.
          actionButton("go_btn", "8 · Generate Heatmap",
                       icon  = icon("fire"),
                       class = "btn-primary btn-lg",
                       width = "100%"))
      ),
      
      conditionalPanel(
        condition = "input.active_feature === 'rnaseq'",
        # start_step continues the numbering begun by the shared "1 Choose
        # species" panel above: the DDS module owns 2-4, the DE module 5-6.
        rnaseqDdsUI("rnaseq_dds", start_step = 2),
        hr(),
        rnaseqDeUI("rnaseq_de", start_step = 5)
      ),
      
      conditionalPanel(
        condition = "input.active_feature === 'go'",
        goSectionSidebarUI("go_section")
      ),
      
      conditionalPanel(
        condition = "input.active_feature === 'venn'",
        vennManualSidebarUI("venn_manual")
      )
    ),
    mainPanel(
      width = 9,
      conditionalPanel(
        condition = "input.active_feature === 'heatmap' || input.active_feature === undefined || input.active_feature === null",
        tags$script(HTML("
$(document).one('shiny:idle', function() {
var cb = document.getElementById('ht_shiny_show_row_names_checkbox');
if (cb && cb.checked) {
cb.checked = false;
$(cb).trigger('change');
}
});
")),
        InteractiveComplexHeatmapOutput("ht_shiny",
                                        height1 = 700,
                                        width1 = 800,
                                        action = "click",
                                        response = c("click", "brush")),
        hr(),
        # ── Appearance & export ─────────────────────────────────────────────
        # Everything here is judged by LOOKING at the figure, so it lives next
        # to the figure rather than in the setup column. Still Generate-driven:
        # these inputs are read inside the go_btn observer exactly as before, so
        # a change takes effect on the next "Generate Heatmap". Shown by default
        # rather than hidden, because a tabsetPanel only ever renders one tab and
        # costs far less height here than the ~400px these controls took in the
        # sidebar -- and because controls nobody can find are worse than clutter.
        actionLink("toggle_style", "▾ \U0001f3a8 Appearance & export",
                   style = paste0("font-size:15px; font-weight:bold; color:#1F3864; ",
                                  "text-decoration:none; display:inline-block; ",
                                  "margin-bottom:6px;")),
        div(
          id = "style_block",
          tabsetPanel(
            tabPanel(
              "\U0001f3a8 Colours",
              br(),
              helpText("The colour ramp of the heatmap body, and the range it spans."),
              fluidRow(
                column(3, colourpicker::colourInput("col_low",  "Low",  value = "#00008b")),
                column(3, colourpicker::colourInput("col_mid",  "Mid",  value = "white")),
                column(3, colourpicker::colourInput("col_high", "High", value = "#FFFF00")),
                column(3, sliderInput("zlim", "Z-score scale limits",
                                      min = -10, max = 10, value = c(-1, 1), step = 0.5))
              )
            ),
            tabPanel(
              "\U0001f524 Labels",
              br(),
              helpText("What is written down the side of the heatmap, and how big."),
              fluidRow(
                column(4, checkboxInput("show_labels", "Show gene labels", value = TRUE)),
                column(4, sliderInput("fontsize_labels", "Label font size",
                                      min = 8, max = 24, value = 14, step = 1)),
                # The species observer calls
                # shinyjs::show/hide("map_gene_names_block") when the active
                # annotation has no symbol column, so this id is load-bearing.
                column(4,
                       div(id = "map_gene_names_block",
                           checkboxInput("map_gene_names",
                                         "Map gene IDs to gene symbols",
                                         value = FALSE)),
                       uiOutput("map_gene_names_help"))
              )
            ),
            tabPanel(
              "\U0001f3f7️ Annotation bars",
              br(),
              helpText("The coloured strips above the heatmap: how tall they are, ",
                       "and which colour each level of each factor gets."),
              fluidRow(
                column(4, sliderInput("anno_bar_size", "Annotation bar height (mm)",
                                      min = 2, max = 15, value = 5, step = 1))
              ),
              uiOutput("annot_colors_ui")
            ),
            tabPanel(
              "⬇ Export",
              br(),
              helpText("Vector PDF at the exact size you set — for figures going into a paper."),
              fluidRow(
                column(3, numericInput("main_ht_w", "Width (in)",
                                       value = 10, min = 2, max = 40, step = 0.5)),
                column(3, numericInput("main_ht_h", "Height (in)",
                                       value = 12, min = 2, max = 60, step = 0.5)),
                column(4, br(),
                       downloadButton("dl_main_ht_pdf", "⬇ Main heatmap (PDF)",
                                      class = "btn-sm btn-success"))
              )
            )
          )
        ),
        hr(),
        uiOutput("cluster_tools_panel"),
        uiOutput("enrichment_panel"),
        # "cart_heatmap_panel" must appear exactly once on the page. The scroll
        # helper in the go_cart_signal observer calls
        # document.getElementById('cart_heatmap_panel'), which returns only the
        # FIRST match -- so a second element sharing this id would silently
        # capture the scroll and the panel would never be revealed.
        uiOutput("cart_heatmap_panel")
      ),
      conditionalPanel(
        condition = "input.active_feature == 'rnaseq'",
        # The whole EDA block folds, so the DE table below can be reached without
        # scrolling past five tabs of diagnostics. section = TRUE keeps the
        # summary at h3 size, matching "Differential expression table" below --
        # otherwise the folded section's title would read as a sub-heading of it.
        # Expanded by default: the fold is an option, not a new default.
        collapsible_panel(
          "Exploratory Data Analysis",
          tabsetPanel(
            tabPanel("PCA",                  rnaseqVizPlotUI("rnaseq_viz")),
            tabPanel("PCA (batch-corrected)", rnaseqEdaPcaCorrectedUI("rnaseq_eda")),
            tabPanel("MDS",                  rnaseqEdaMdsUI("rnaseq_eda")),
            tabPanel("Sample correlation",   rnaseqEdaCorUI("rnaseq_eda")),
            tabPanel("Batch diagnostics",    rnaseqEdaBatchUI("rnaseq_eda"))
          ),
          open = TRUE, heading = "h3"
        ),
        br(),
        h3("Differential expression table"),
        fluidRow(
          column(
            12,
            div(
              style = "display: flex; gap: 8px; align-items: center; margin-bottom: 10px;",
              actionButton(
                "open_de_table",
                "📋 Table",
                class = "btn-default btn-sm"
              ),
              actionButton(
                "open_go_de",
                "🧬 GO Enrichment (DEGs)",
                class = "btn-primary btn-sm"
              ),
              actionButton(
                "open_kegg_de",
                "🗺️ KEGG Enrichment (DEGs)",
                class = "btn-warning btn-sm"
              ),
              actionButton("open_gsea_de", "\U0001f4c8 GSEA (ranked)", class = "btn-info btn-sm"),
              actionButton(
                "open_de_hist",
                "💾 Save / View gene lists",
                class = "btn-success btn-sm"
              ),
              actionButton(
                "open_venn",
                "🔵 Venn from selections",
                class = "btn-info btn-sm"
              ),
              rnaseqVolcanoUI("rnaseq_volcano")
            ),
            vennModalUI("venn")
          )
        ),
        uiOutput("rnaseq_main_panel")
      ),
      
      # GO Analysis main panel — standalone Gene Ontology enrichment section
      conditionalPanel(
        condition = "input.active_feature == 'go'",
        goSectionMainUI("go_section")
      ),

      # Venn manual main panel — full module UI rendered here
      conditionalPanel(
        condition = "input.active_feature == 'venn'",
        vennManualMainUI("venn_manual")
      )
    )
  )
)
# ════════════════════════════════════════════════════════════════════════════════
server <- function(input, output, session) {
  
  register_autoclose_on_browser_close(session)

  n_groups <- reactiveVal(4)
  n_annots <- reactiveVal(2)

  # ── Experiment-setup fold state ────────────────────────────────────────────
  # FALSE = expanded. Flipped to TRUE automatically after a heatmap is
  # successfully generated (see the go_btn observer), on the assumption that
  # once the experiment is described you are done editing it. Purely visual:
  # every input below keeps its value and is read by the server unchanged.
  setup_collapsed <- reactiveVal(FALSE)

  observeEvent(input$toggle_setup, {
    setup_collapsed(!isTRUE(setup_collapsed()))
  }, ignoreInit = TRUE)

  # Single place that applies the state, so the caret, the panel and the
  # summary can never disagree. The link is created ONCE in the UI and only
  # relabelled here -- rebuilding it inside a renderUI would reset its click
  # counter on every fold and re-fire the observer above.
  observeEvent(setup_collapsed(), {
    collapsed <- isTRUE(setup_collapsed())
    shinyjs::toggle("setup_block", condition = !collapsed)
    updateActionButton(
      session, "toggle_setup",
      label = paste0(if (collapsed) "▸" else "▾",
                     " 4 · \U0001f4e6 Experiment setup"))
  }, ignoreNULL = FALSE)

  output$setup_summary <- renderUI({
    if (!isTRUE(setup_collapsed())) return(NULL)
    sel      <- lapply(seq_len(n_groups()), function(i) input[[paste0("grp_cols_", i)]])
    n_valid  <- sum(vapply(sel, function(x) !is.null(x) && length(x) > 0, logical(1)))
    n_used   <- length(unique(unlist(sel)))
    tagList(
      tags$p(style = "color:#1F3864; font-size:0.88em; margin:0 0 8px 16px;",
             sprintf("%d group(s) · %d sample(s) · %d annotation factor(s)",
                     n_valid, n_used, n_annots()),
             tags$span(style = "color:#999;", " — click the title to edit"))
    )
  })

  # ══ Analysis / appearance split ════════════════════════════════════════════
  # The heatmap is built in two halves so that a cosmetic change never triggers
  # a recompute. Building it entirely inside observeEvent(input$go_btn) would
  # mean re-running VST, filtering, group means, normalization and clustering
  # for a colour change -- 2-4s on 1000 genes, far more on a full matrix -- and
  # would append a heatmap_history entry each time, filling the reproducibility
  # bundle with near-duplicate matrix CSVs and rerun scripts.
  #
  # The two halves:
  #   analysis_rv  <- the expensive result (matrices, clustering, cut). Written
  #                   ONLY by the go_btn observer.
  #   draw observer -> rebuilds the Heatmap object from analysis_rv() plus the
  #                   live appearance inputs. Cheap, no recompute.
  #
  # The dividing line is exactly the one the UI already draws after phases 2-3:
  #   SIDEBAR  = analysis -> needs Generate, appends a history entry
  #   APPEARANCE PANEL = style -> applies live, updates the last entry in place
  # so "which controls need Generate?" has a visible answer rather than a
  # remembered one.
  analysis_rv <- reactiveVal(NULL)

  # Everything that changes the MATRIX. Hashed to decide (a) whether the figure
  # on screen is stale and (b) whether a Generate is a new analysis or a redraw
  # of the same one. Style inputs are deliberately absent: they cannot make a
  # figure stale because they are applied immediately.
  data_spec <- reactive({
    list(
      file          = if (is.null(input$file)) NULL else input$file$datapath,
      demo          = !is.null(demo_df()),
      is_raw_counts = isTRUE(input$is_raw_counts),
      gene_filter   = input$gene_filter,
      top_var_n     = input$top_var_n,
      normalization = input$normalization,
      ref_condition = input$ref_condition,
      cluster_rows  = isTRUE(input$cluster_rows),
      cluster_cols  = isTRUE(input$cluster_cols),
      linkage       = input$linkage,
      cut_tree      = isTRUE(input$cut_tree),
      k_clusters    = input$k_clusters,
      show_degs     = isTRUE(input$show_degs),
      deg_list      = if (isTRUE(input$show_degs)) input$deg_list else NULL,
      n_groups      = n_groups(),
      groups        = lapply(seq_len(n_groups()), function(i)
                        list(input[[paste0("grp_name_", i)]],
                             input[[paste0("grp_cols_", i)]])),
      n_annots      = n_annots(),
      annots        = lapply(seq_len(n_annots()), function(k)
                        list(input[[paste0("cat_name_", k)]],
                             input[[paste0("cat_levels_", k)]],
                             vapply(seq_len(n_groups()),
                                    function(i) input[[paste0("assign_", i, "_", k)]] %||% "", "")))
    )
  })
  data_digest <- reactive(digest::digest(data_spec()))

  # Says out loud which of the two states you are in. Without it the appearance
  # controls apply instantly while the sidebar ones silently do not, and a
  # figure can disagree with its own sidebar with nothing on screen saying so.
  output$stale_badge <- renderUI({
    a <- analysis_rv()
    if (is.null(a)) {
      return(tags$div(style = "font-size:0.82em; color:#777; margin-bottom:6px;",
                      "No heatmap yet — set your groups and click Generate."))
    }
    if (identical(a$data_digest, data_digest())) {
      tags$div(style = "font-size:0.82em; color:#3c763d; margin-bottom:6px;",
               "✓ Figure matches the current settings.")
    } else {
      tags$div(style = paste0("font-size:0.82em; color:#a94442; font-weight:bold; ",
                              "margin-bottom:6px;"),
               "⚠ Analysis settings changed — click Generate to update.")
    }
  })

  # ── Appearance panel (main panel) and Session tools (sidebar) folds ────────
  # Same one-reactiveVal-per-fold pattern as setup_collapsed above: the link is
  # built once in the UI and only relabelled here, so its click counter is never
  # reset. Appearance starts OPEN (it is new, and a tabsetPanel is cheap);
  # Session tools starts CLOSED (its UI is wrapped in shinyjs::hidden()).
  style_open   <- reactiveVal(TRUE)
  session_open <- reactiveVal(FALSE)

  observeEvent(input$toggle_style,   { style_open(!isTRUE(style_open())) },
               ignoreInit = TRUE)
  observeEvent(input$toggle_session, { session_open(!isTRUE(session_open())) },
               ignoreInit = TRUE)

  observeEvent(style_open(), {
    open <- isTRUE(style_open())
    shinyjs::toggle("style_block", condition = open)
    updateActionButton(session, "toggle_style",
                       label = paste0(if (open) "▾" else "▸",
                                      " \U0001f3a8 Appearance & export"))
  }, ignoreNULL = FALSE)

  observeEvent(session_open(), {
    open <- isTRUE(session_open())
    shinyjs::toggle("session_block", condition = open)
    updateActionButton(session, "toggle_session",
                       label = paste0(if (open) "▾" else "▸",
                                      " \U0001f9f0 Session tools"))
  }, ignoreNULL = FALSE)

  species_rv <- reactiveVal(ACTIVE_SPECIES)
  
  observeEvent(input$apply_species, {
    req(input$species_switch)
    new_sp <- input$species_switch
    if (new_sp != species_rv()) {
      withProgress(message = paste0("Loading ", new_sp, " config..."), {
        tryCatch({
          cfg <- get_species_config(new_sp)
          list2env(cfg, envir = .GlobalEnv)
          assign("ACTIVE_SPECIES", new_sp, envir = .GlobalEnv)
          clear_enrichment_cache()   # avoid serving stale cross-species results
          scaled_matrix_rv(NULL)     # matrix was canonicalized under old rules
          selected_genes(NULL)
          cart_heatmap_genes(NULL)
          # ── Everything else keyed to the OLD species ──────────────────────
          # Only the four above were cleared, so the rest survived a switch and
          # kept describing the previous organism. The concrete failures:
          #   * cluster_assign_rv kept the Clusters panel rendered, and its
          #     buttons re-injected the old organism's gene IDs into
          #     selected_genes() -- which then reached GO/KEGG under the NEW
          #     annotation and came back empty, reported as "No terms found
          #     below the chosen cutoffs" rather than as an ID mismatch.
          #   * main_ht_rv made "Download main heatmap (PDF)" export the
          #     previous organism's figure.
          #   * heatmap_history would have mixed two organisms in one bundle.
          # Listed explicitly rather than looped, so adding a new reactiveVal
          # forces a deliberate decision about whether it survives a switch.
          cluster_matrix_rv(NULL)
          cluster_assign_rv(NULL)
          gene_labels_rv(NULL)
          current_gene_ids(NULL)
          col_fun_rv(NULL)
          col_order_rv(NULL)
          main_ht_rv(NULL)
          cluster_stability_rv(NULL)
          heatmap_history(list())
          if (!is.null(rnaseq_handles$reset)) rnaseq_handles$reset()
          species_rv(new_sp)
          showNotification(
            paste0("\u2705 Species switched to ", cfg$SPECIES_COMMON),
            type = "message", duration = 4
          )
        }, error = function(e) {
          showNotification(
            paste0("\u274c Failed to load species config: ", e$message),
            type = "error", duration = 8
          )
        })
      })
    }
  })
  observeEvent(input$add_group,    { if (n_groups() < MAX_GROUPS) n_groups(n_groups() + 1) })
  output$species_label_ui <- renderUI({
    species_rv()  # reactive dependency
    tags$p(
      style = "font-weight:bold; color:#1F3864; margin-bottom:6px;",
      SPECIES_COMMON
    )
  })
  
  output$map_gene_names_help <- renderUI({
    species_rv()  # reactive dependency
    if (!has_gene_symbols()) {
      # State why the control is absent, rather than leaving a silent gap.
      return(tags$small(
        style = "color:#8a6d3b;",
        paste0(ORG_DB_PKG, " carries no gene-symbol column, so ", SPECIES_COMMON,
               " gene IDs cannot be mapped to symbols. IDs are shown as-is.")))
    }
    tagList(
      tags$small(paste0("Uses ", ORG_DB_PKG, ". Gene ID kept as fallback."))
    )
  })

  # ── Species-dependent sidebar refresh ──────────────────────────────────────
  # The sidebar is built once at app startup, so anything read from a species
  # global at UI-construction time keeps describing the previous organism after
  # "Switch species & reload". Panels rendered through renderUI (the help text
  # above, the enrichment panel hosting tairGenesUI) re-render on their own;
  # the two controls below do not, and are refreshed here.
  #
  # Symbol mapping: "Map gene IDs to gene symbols" is a no-op whenever ORG_DB
  # carries no symbol column, which is the case for org.Slycopersicum.eg.db.
  # The control is hidden AND forced off, so a heatmap recorded after a switch
  # cannot claim a mapping that never happened.
  #
  # ID placeholders: the two gene-ID textareas otherwise keep the previous
  # organism's example format (AT1G01010... vs Solyc01g005000.1...). Only the
  # hint is refreshed; updateTextAreaInput() with `value` left NULL preserves
  # whatever the user has already typed, which is why this is an update rather
  # than a re-render of the panel.
  observeEvent(species_rv(), {
    if (has_gene_symbols()) {
      shinyjs::show("map_gene_names_block")
    } else {
      shinyjs::hide("map_gene_names_block")
      updateCheckboxInput(session, "map_gene_names", value = FALSE)
    }
    for (id in c("deg_list", "universe_custom_genes")) {
      updateTextAreaInput(session, id, placeholder = GENE_ID_PLACEHOLDER)
    }
  }, ignoreNULL = FALSE)
  observeEvent(input$remove_group, { if (n_groups() > 1)          n_groups(n_groups() - 1) })
  observeEvent(input$add_annot,    { if (n_annots() < MAX_ANNOTS) n_annots(n_annots() + 1) })
  observeEvent(input$remove_annot, { if (n_annots() > 1)          n_annots(n_annots() - 1) })

  selected_genes     <- reactiveVal(NULL)
  current_gene_ids   <- reactiveVal(NULL)
  scaled_matrix_rv   <- reactiveVal(NULL)
  # cluster_matrix_rv: VST/expression-scale (pre-display-normalization)
  # condition-mean matrix, used ONLY as input to hclust()/cor(). Kept
  # separate from scaled_matrix_rv (display/color only) so a cosmetic choice
  # in "Display transform" can never change a dendrogram. Row/column-aligned
  # with scaled_matrix_rv at all times 
  cluster_matrix_rv  <- reactiveVal(NULL)
  gene_labels_rv     <- reactiveVal(NULL)
  cart_heatmap_genes <- reactiveVal(NULL)
  col_order_rv       <- reactiveVal(NULL)
  col_fun_rv         <- reactiveVal(NULL)
  # cluster_assign_rv: named character vector (gene_id -> cluster number) from
  # the last cutree(). Held at this level, not inside the go_btn observer, so
  # the partition is usable as well as visible: per-cluster enrichment, CSV
  # export and the reproducibility bundle all read it. NULL whenever
  # "Cut dendrogram" is off or clustering was refused.
  cluster_assign_rv  <- reactiveVal(NULL)
  # main_ht_rv: the last main Heatmap object, kept so it can be re-drawn into a
  # PDF device at an arbitrary size without recomputing the analysis.
  main_ht_rv         <- reactiveVal(NULL)
  # heatmap_history: one entry per successful "Generate Heatmap", each holding
  # the config AND the plotted matrix as they were at that moment. Every other
  # reactiveVal here is overwritten by the next run; this one accumulates, which
  # is what lets the bundle describe every figure of the session rather than the
  # last, and describe each with the settings that produced it rather than with
  # whatever the sidebar says at download time.
  # A list, like go_history/kegg_history/gsea_history, so N runs give N records.
  # Memory: one display matrix per run (genes x conditions, conditions <= 12).
  # A 5000-gene heatmap is ~0.5 MB, so a long session costs a few MB.
  heatmap_history    <- reactiveVal(list())
  scale_label_rv     <- reactiveVal("Row Z-score")
  demo_df            <- reactiveVal(NULL)

  

  deg_genes_rv <- reactive({
    if (!isTRUE(input$show_degs)) return(character(0))
    
    raw <- input$deg_list
    
    # message("DEBUG raw deg_list (nchar): ", nchar(raw))
    # message("DEBUG raw deg_list (dput): ")
    # message(paste(dput(raw), collapse = ""))
    
    if (is.null(raw) || !nzchar(trimws(raw))) return(character(0))
    
    parse_id_textarea(clean_pasted_ids(raw))
  })
  
  observeEvent(input$clear_gene_filter, {
    updateTextAreaInput(session, "gene_filter", value = "")
  })
  
  observeEvent(input$clear_deg_list, {
    updateTextAreaInput(session, "deg_list", value = "")
  })
  
  mat_data <- reactive({
    # Gene IDs are read into a separate vector rather than straight into row
    # names. A data.frame cannot hold duplicate row names, so reading them as
    # row names makes a repeated gene ID -- ordinary in a per-transcript export
    # -- fail at the READ step, before validation, with a bare
    # "duplicate 'row.names' are not allowed". The aggregation below is designed
    # to handle exactly that case, so the read must let it through.
    if (is.null(input$file) && !is.null(demo_df())) {
      df  <- demo_df()
      ids <- rownames(df)
      import_notes <- character(0)
    } else {
      req(input$file)
      imported     <- read_counts_file(input$file$datapath, input$file$name)
      ids          <- imported$ids
      df           <- imported$values
      import_notes <- imported$notes
      # Placeholder names only, so validate_count_matrix() can run its own ID
      # checks. The aggregation below uses `ids`, never these.
      rownames(df) <- make.unique(ids, sep = "..dup")
    }

    for (n in import_notes) {
      showNotification(paste0("ℹ️ ", n), type = "message", duration = 10)
      log_warning(n)
    }

    v <- validate_count_matrix(df, require_counts = isTRUE(input$is_raw_counts),
                               source_label = "uploaded file")
    notify_and_log_validation(v)  # shows notifications AND writes them to the
    # downloadable session log (utils/log_utils.R)
    req(v$ok)                     # stop the reactive cleanly if unusable
    # Only COLUMNS are ever dropped by validation, so `ids` is still aligned
    # with v$matrix row for row.
    df <- as.data.frame(v$matrix, check.names = FALSE)

    collapse_counts_by_id(df, ids)
  })
  
  rep_cols <- reactive({
    req(mat_data())
    colnames(mat_data())
  })
  
  universe_genes <- reactive({
    mode <- input$universe_mode %||% "default"
    
    if (identical(mode, "custom")) {
      raw <- input$universe_custom_genes
      if (is.null(raw) || !nzchar(trimws(raw))) return(NULL)
      ids <- parse_id_textarea(clean_pasted_ids(raw))
      if (length(ids) == 0) return(NULL)
      return(ids)
    }
    
    if (identical(mode, "threshold")) {
      req(mat_data())
      
      df <- mat_data()
      expr <- as.matrix(df)
      storage.mode(expr) <- "numeric"
      
      thr <- as.numeric(input$universe_threshold)
      if (is.na(thr)) thr <- 0
      
      grp_ids <- grep("^grp_cols_", names(input), value = TRUE)
      group_sizes <- vapply(grp_ids, function(id) {
        sel <- input[[id]]
        if (is.null(sel) || length(sel) == 0) 0L else length(sel)
      }, integer(1))
      group_sizes <- group_sizes[group_sizes > 0]
      min_group_size_heatmap <- if (length(group_sizes) > 0) min(group_sizes) else 1L
      
      keep <- rowSums(expr >= thr, na.rm = TRUE) >= min_group_size_heatmap
      
      return(rownames(expr)[keep])
    }
    
    NULL   # "default": no custom universe passed; GO uses the whole-genome background
  })
  

  # ── Dynamic reference-condition selector ──────────────────────────────────────
  # Re-derives choices from current group names each time they are edited.
  # Visible only when normalization == "log2fc" (conditionalPanel in UI).
  output$ref_condition_ui <- renderUI({
    n         <- n_groups()
    grp_names <- vapply(seq_len(n), function(i) {
      input[[paste0("grp_name_", i)]] %||% paste("Group", i)
    }, "")
    # Keep the reference the user chose across a re-render. This block lists the
    # GROUP NAMES as its choices, so it re-renders whenever a group is renamed or
    # added -- and it used to re-seed `selected` from grp_names[1] every time,
    # silently moving the log2FC baseline back to group 1 after an unrelated
    # edit. Every value in the figure is relative to that baseline, so the next
    # Generate produced a different matrix from the one that was asked for. The
    # result was never inconsistent (the notification and the recorded config
    # both named the reference actually used) -- it simply was not the user's
    # choice any more.
    #
    # isolate(): required, not cosmetic. An un-isolated read of `ref_condition`
    # -- an input this block CREATES -- would rebuild the selector every time it
    # changed, closing the dropdown mid-selection.
    #
    # intersect(): required too. `multiple = TRUE`, and the choices ARE the group
    # names, so a renamed or removed group leaves a stale entry behind. Passing a
    # value that is not among `choices` selects NOTHING, and
    # log2fc_vs_reference() stop()s on an empty reference -- hence the explicit
    # fallback to group 1.
    #
    # Residual case, deliberately not handled here: renaming THE REFERENCE GROUP
    # ITSELF still falls back to group 1, because the reference is tracked by
    # name. Tracking it by index would mean changing normalize_matrix() and
    # drop_reference_columns(), which both match references by name.
    prev <- intersect(isolate(input$ref_condition), grp_names)
    selectInput(
      inputId  = "ref_condition",
      label    = "Reference condition(s)",
      choices  = grp_names,
      selected = if (length(prev)) prev else grp_names[1],
      multiple = TRUE,
      width    = "100%"
    )
  })

  # ── Preview matrix for optimal-k (reactive, no Generate button needed) ────────
  preview_matrix_for_optk <- reactive({
    req(mat_data())
    df       <- mat_data()
    gene_ids <- rownames(df)
    expr_mat <- as.matrix(df)
    storage.mode(expr_mat) <- "numeric"

    # VST for raw counts — silent in preview reactive (avoids notification spam).
    if (isTRUE(input$is_raw_counts)) {
      expr_mat <- tryCatch({
        counts_int <- round(expr_mat)
        counts_int[counts_int < 0] <- 0L
        storage.mode(counts_int) <- "integer"
        dds <- DESeq2::DESeqDataSetFromMatrix(
          countData = counts_int,
          colData   = data.frame(row.names = colnames(counts_int),
                                 grp       = factor(rep("A", ncol(counts_int)))),
          design    = ~ 1)
        # vst() has TWO guards that both stop(): nrow(dds) < nsub, and
        # sum(baseMean > 5) < nsub. Default nsub = 1000, and nothing pre-filters
        # this matrix, so a sparse or shallow matrix trips the second guard at
        # any row count and a small gene list trips the first. Falling back to
        # varianceStabilizingTransformation() (no nsub, fits the full trend)
        # keeps the data variance-stabilised instead of dropping to raw counts.
        vsd <- tryCatch(
          DESeq2::vst(dds, blind = TRUE, nsub = min(1000, nrow(dds))),
          error = function(e) DESeq2::varianceStabilizingTransformation(dds, blind = TRUE)
        )
        vst_out <- as.matrix(SummarizedExperiment::assay(vsd))
        vst_out
      }, error = function(e) expr_mat)  # fallback: use unchanged matrix
    }

    gene_filter_raw <- input$gene_filter
    if (!is.null(gene_filter_raw) && nzchar(trimws(gene_filter_raw))) {
      filter_ids <- parse_id_textarea(clean_pasted_ids(gene_filter_raw))
      keep <- ids_in(gene_ids, filter_ids)
      
      if (sum(keep) == 0) {
        log_warning(sprintf(
          "k-suggestion preview: 0 of %d filtered gene ID(s) matched the uploaded matrix.",
          length(filter_ids)))
        return(NULL)
      }
      
      df <- df[keep, , drop = FALSE]
      gene_ids <- gene_ids[keep]
      expr_mat <- expr_mat[keep, , drop = FALSE]
    }

    n <- n_groups()
    
    # Same source as the go_btn path and as output$ref_condition_ui, so the
    # colnames set below are the exact strings input$ref_condition holds.
    group_names <- vapply(seq_len(n),function(i) input[[paste0("grp_name_", i)]] %||% paste("Group", i), "")
    
    group_means <- lapply(seq_len(n), function(i) {
      sel <- input[[paste0("grp_cols_", i)]]
      if (is.null(sel) || length(sel) == 0) return(NULL)
      rowMeans(expr_mat[, sel, drop = FALSE], na.rm = TRUE)
    })
    valid <- !sapply(group_means, is.null)
    if (sum(valid) == 0) {
      log_warning("k-suggestion preview: no groups have columns assigned yet.")
      return(NULL)
    }
    valid_idx <- which(valid)

    # Same restriction as the go_btn path: only columns assigned to a group take
    # part in the analysis.
    used_cols <- unique(unlist(lapply(valid_idx, function(i) {
      input[[paste0("grp_cols_", i)]]
    })))

    mean_mat <- do.call(cbind, group_means[valid])
    rownames(mean_mat) <- gene_ids
    colnames(mean_mat) <- group_names[valid]

    # Remove constant rows (crash silhouette / kmeans)
    keep <- apply(mean_mat, 1, function(x) {
      v <- var(x, na.rm = TRUE)
      !is.na(v) && v >= 1e-10
    })
    if (sum(keep) < 4) {
      log_warning(sprintf(
        "k-suggestion preview: only %d gene(s) have non-zero variance across groups (need >= 4) — elbow/silhouette not shown.",
        sum(keep)))
      return(NULL)
    }
    mean_mat <- mean_mat[keep, , drop = FALSE]

    # Same top-N filter as the go_btn path, so the suggested k is scored on the
    # same gene set the heatmap will actually cluster. Applied after the
    # constant-row drop above, matching the order in the go_btn observer.
    top_n_optk <- suppressWarnings(as.integer(input$top_var_n %||% 0))
    if (!is.na(top_n_optk) && top_n_optk > 0 && top_n_optk < nrow(mean_mat)) {
      mean_mat <- mean_mat[top_variable_rows(mean_mat, top_n_optk), , drop = FALSE]
    }

    # ── k-suggestion clustering substrate ─────────────────────────────────────
    # mean_mat above is used only to SELECT rows (constant-row drop, top-N), the
    # same way the go_btn path does. The substrate returned is the sample-level
    # matrix, ROW Z-SCORED, because that is exactly what the heatmap clusters:
    # scoring k on a different matrix would score a different tree than the one
    # being cut. Kept in step with the `clust_mat` block in the go_btn observer
    # -- if one side changes, so must the other. The log2fc reference-column drop
    # does not apply here: the columns are samples, not group means, and the
    # go_btn path clusters on all assigned samples.
    tryCatch({
      sub <- expr_mat[match(rownames(mean_mat), rownames(expr_mat)), used_cols, drop = FALSE]
      zz  <- t(scale(t(sub)))
      if (anyNA(zz) || !all(is.finite(zz))) sub else zz
    }, error = function(e) NULL)
  })

  output$universe_info <- renderUI({
    req(input$file)
    req(identical(input$universe_mode, "threshold"))
    u <- universe_genes()
    if (is.null(u)) return(NULL)
    tags$p(
      style = "color: #2c7bb6; font-weight: bold; margin-top: 4px;",
      icon("circle-info"),
      paste0(" ", length(u), " genes in universe")
    )
  })
  
  output$universe_custom_info <- renderUI({
    req(identical(input$universe_mode, "custom"))
    u <- universe_genes()
    if (is.null(u)) {
      return(tags$p(style = "color: #999; font-size: 0.85em; margin-top: 4px;",
                    "Paste at least one gene ID above."))
    }
    tags$p(
      style = "color: #2c7bb6; font-weight: bold; margin-top: 4px;",
      icon("circle-info"),
      paste0(" ", length(u), " genes in custom universe")
    )
  })
  
  observeEvent(input$clear_universe_custom, {
    updateTextAreaInput(session, "universe_custom_genes", value = "")
  })

  output$group_ui <- renderUI({
    req(rep_cols())
    lapply(seq_len(n_groups()), function(i) {
      wellPanel(
        style = "padding: 8px; margin-bottom: 6px;",
        fluidRow(
          # isolate(): re-read the CURRENT values so a re-render (Add group /
          # Remove last / a new upload) keeps what the user set. Without it this
          # block re-seeded every box from its hard-coded default, so clicking
          # "Add group" renamed every group back to Group1..GroupN and cleared
          # the replicate selections. isolate() is what stops the read from
          # making this renderUI depend on the inputs it creates -- that
          # dependency would rebuild the panel on every keystroke and take the
          # cursor with it. Same idiom as annot_colors_ui below.
          #
          # `selected` on the selectize also removes a race: the mutual-exclusion
          # observer further down calls updateSelectizeInput(selected = ...) and
          # used to be the only thing that could restore a selection after a
          # re-render. It still enforces exclusivity; it is no longer what keeps
          # the selection alive.
          #
          # A NEW upload changes rep_cols(), so the stored names are absent from
          # `choices` and Shiny selects nothing -- which correctly clears the
          # groups when the matrix changes.
          column(5, textInput(paste0("grp_name_", i),
                              label = paste("Group", i, "name"),
                              value = isolate(input[[paste0("grp_name_", i)]]) %||%
                                      paste0("Group", i))),
          column(7, selectizeInput(paste0("grp_cols_", i),
                                   label    = "Select replicates",
                                   choices  = rep_cols(),
                                   selected = isolate(input[[paste0("grp_cols_", i)]]),
                                   multiple = TRUE,
                                   options  = list(
                                     plugins     = list("remove_button"),
                                     placeholder = "Click to select..."
                                   )))
        )
      )
    })
  })
  
  observe({
    all_cols <- rep_cols()
    req(all_cols)
    n <- n_groups()
    all_selected <- lapply(seq_len(n), function(i) {
      val <- input[[paste0("grp_cols_", i)]]
      if (is.null(val)) character(0) else val
    })
    for (i in seq_len(n)) {
      others    <- unlist(all_selected[-i])
      avail     <- setdiff(all_cols, others)
      sel_clean <- intersect(all_selected[[i]], avail)
      updateSelectizeInput(session,
                           inputId  = paste0("grp_cols_", i),
                           choices  = avail,
                           selected = sel_clean)
    }
  })

  output$annot_cat_ui <- renderUI({
    lapply(seq_len(n_annots()), function(k) {
      wellPanel(
        style = "padding: 8px; margin-bottom: 4px;",
        fluidRow(
          # isolate(): same reason as group_ui above. The values below are only
          # EXAMPLES -- everyone renames them -- so re-seeding them from the
          # hard-coded vector on every re-render destroyed real input: clicking
          # "Add factor" put "Temperature" / "T22,T32" back over whatever the
          # user had typed. The defaults now apply only when the input does not
          # exist yet.
          #
          # %||% (global.R:42) also treats "" as absent, so deliberately
          # clearing a box and then re-rendering brings the example back. That is
          # the accepted trade for matching the idiom used by annot_colors_ui;
          # swap in an explicit is.null() test if an empty box should stay empty.
          column(5, textInput(paste0("cat_name_", k),
                              label = paste("Factor", k, "name"),
                              value = isolate(input[[paste0("cat_name_", k)]]) %||%
                                      c("Temperature","Treatment","Genotype","Condition")[k])),
          column(7, textInput(paste0("cat_levels_", k),
                              label = "Levels (comma-separated)",
                              value = isolate(input[[paste0("cat_levels_", k)]]) %||%
                                      c("T22,T32","mock,B82T","WT,mutant","control,treated")[k]))
        )
        # The per-level colour pickers are NOT here: they live in the
        # "Annotation bars" tab of the appearance panel. They are the most
        # numerous purely-cosmetic control in the app -- one per level, so a
        # third factor with six levels adds six more -- and belong next to the
        # figure they affect rather than in the setup column. Their inputIds
        # remain cat_col_<k>_<j>, which is what get_annot_colors() and the
        # reproducibility bundle read.
      )
    })
  })

  # All factors' level colours in one output, rendered next to the figure.
  # Consolidated from the previous observe() + local() loop that assigned one
  # output per factor: with a single container there is nothing to keep in sync
  # when n_annots() changes. The isolate() is kept for the reason it was there
  # originally -- it re-reads each picker's CURRENT value so editing one
  # factor's levels does not reset colours the user already chose.
  output$annot_colors_ui <- renderUI({
    blocks <- lapply(seq_len(n_annots()), function(kk) {
      lvls <- trimws(strsplit(input[[paste0("cat_levels_", kk)]] %||% "", ",")[[1]])
      lvls <- lvls[nzchar(lvls)]
      if (length(lvls) == 0) return(NULL)
      def_pal  <- annotation_palettes[[((kk - 1) %% length(annotation_palettes)) + 1]]
      cat_name <- input[[paste0("cat_name_", kk)]] %||% paste("Factor", kk)
      tagList(
        tags$small(tags$b(paste0("\U0001f3a8 ", cat_name, " — level colours:"))),
        fluidRow(lapply(seq_along(lvls), function(j) {
          def_col <- if (j <= length(def_pal)) def_pal[j] else "#AAAAAA"
          column(
            max(2, floor(12 / min(length(lvls), 6))),
            colourpicker::colourInput(
              inputId = paste0("cat_col_", kk, "_", j),
              label   = lvls[j],
              value   = isolate(
                input[[paste0("cat_col_", kk, "_", j)]] %||% def_col
              )
            )
          )
        }))
      )
    })
    blocks <- Filter(Negate(is.null), blocks)
    if (length(blocks) == 0)
      return(helpText("Define at least one annotation factor with levels in the sidebar."))
    do.call(tagList, blocks)
  })

  output$annot_assign_ui <- renderUI({
    lapply(seq_len(n_groups()), function(i) {
      grp_name <- input[[paste0("grp_name_", i)]] %||% paste("Group", i)
      selects  <- lapply(seq_len(n_annots()), function(k) {
        lvls <- trimws(strsplit(input[[paste0("cat_levels_", k)]] %||% "", ",")[[1]])
        lvls <- lvls[nzchar(lvls)]
        if (length(lvls) == 0) lvls <- "?"
        column(
          floor(12 / max(n_annots(), 1)),
          # selected + isolate(): this block re-renders on EVERY keystroke in a
          # group name, a factor name or a levels string (it reads all three),
          # and selectInput() with no `selected` marks the FIRST choice. So
          # editing any of those silently reset every group's level assignment
          # to lvls[1] -- and that vector is what draws the coloured annotation
          # bars and what lands in the reproducibility bundle, so the figure
          # could end up contradicting its own column names with nothing said.
          #
          # isolate() is required, not cosmetic: an un-isolated read of an input
          # this block CREATES would make the panel re-render on every dropdown
          # change and drop focus mid-interaction.
          #
          # If the levels change so the stored value is no longer among
          # `choices`, Shiny marks nothing selected and the browser falls back to
          # the first option -- the previous behaviour, and the sensible one.
          selectInput(paste0("assign_", i, "_", k),
                      label    = input[[paste0("cat_name_", k)]] %||% paste("Cat", k),
                      choices  = lvls,
                      selected = isolate(input[[paste0("assign_", i, "_", k)]]),
                      width    = "100%")
        )
      })
      wellPanel(
        style = "padding: 6px; margin-bottom: 4px;",
        tags$b(grp_name),
        do.call(fluidRow, selects)
      )
    })
  })

  get_annot_colors <- function(lvls, k) {
    cols <- vapply(seq_along(lvls), function(j) {
      def_pal <- annotation_palettes[[((k - 1) %% length(annotation_palettes)) + 1]]
      def_col <- if (j <= length(def_pal)) def_pal[j] else "#AAAAAA"
      input[[paste0("cat_col_", k, "_", j)]] %||% def_col
    }, "")
    setNames(cols, lvls)
  }

  observeEvent(input$clear_cache, {
    clear_enrichment_cache()
    showNotification("\U0001f5d1\ufe0f Enrichment cache cleared.", type = "message", duration = 3)
  })
  
  observeEvent(input$go_btn, {
    req(mat_data())
    df       <- mat_data()
    gene_ids <- rownames(df)
    expr_mat <- as.matrix(df)
    storage.mode(expr_mat) <- "numeric"
    
    vst_applied <- FALSE

    # VST for raw counts — with user notifications.
    if (isTRUE(input$is_raw_counts)) {
      showNotification(
        "\U0001f9ec Applying DESeq2 VST normalisation\u2026",
        type = "message", id = "vst_msg", duration = NULL)
      expr_mat <- tryCatch({
        counts_int <- round(expr_mat)
        counts_int[counts_int < 0] <- 0L
        storage.mode(counts_int) <- "integer"
        dds <- DESeq2::DESeqDataSetFromMatrix(
          countData = counts_int,
          colData   = data.frame(row.names = colnames(counts_int),
                                 grp       = factor(rep("A", ncol(counts_int)))),
          design    = ~ 1)
        # See the note on the other vst() call site: nsub has two stop()
        # conditions and neither is about the user's data being unusable, so
        # fall back to the full-trend estimator rather than to raw counts.
        # The static "Applying VST..." notification above says it started; this
        # says what it is doing, including the fallback to the full-trend
        # estimator, which is the slow branch. See with_progress_messages() in
        # log_utils.R. The tryCatch/fallback logic is untouched.
        vsd <- with_progress_messages(
          tryCatch(
            DESeq2::vst(dds, blind = TRUE, nsub = min(1000, nrow(dds))),
            error = function(e) DESeq2::varianceStabilizingTransformation(dds, blind = TRUE)
          ),
          "Applying DESeq2 VST")
        vst_out <- as.matrix(SummarizedExperiment::assay(vsd))
        removeNotification("vst_msg")
        showNotification("\u2705 DESeq2 VST normalisation applied.",
                         type = "message", duration = 5)
        
        vst_applied <- TRUE
        vst_out
      }, error = function(e) {
        removeNotification("vst_msg")
        showNotification(
          paste0("\u274c VST failed: ", conditionMessage(e),
                 " \u2014 raw values used."),
          type = "error", duration = 10)
        expr_mat
      })
    }

    gene_filter_raw <- input$gene_filter
    if (!is.null(gene_filter_raw) && nzchar(trimws(gene_filter_raw))) {
      filter_ids <- parse_id_textarea(clean_pasted_ids(gene_filter_raw))
      keep <- ids_in(gene_ids, filter_ids)
      
      if (sum(keep) == 0) {
        showNotification("❌ No genes matched. Check your IDs.", type = "error")
        return()
      }
      
      df <- df[keep, , drop = FALSE]
      gene_ids <- gene_ids[keep]
      expr_mat <- expr_mat[keep, , drop = FALSE]
    }

    n           <- n_groups()
    group_names <- vapply(seq_len(n),
                          function(i) input[[paste0("grp_name_", i)]] %||% paste("Group", i), "")
    group_means <- lapply(seq_len(n), function(i) {
      sel <- input[[paste0("grp_cols_", i)]]
      if (is.null(sel) || length(sel) == 0) return(NULL)
      rowMeans(expr_mat[, sel, drop = FALSE], na.rm = TRUE)
    })
    valid       <- !sapply(group_means, is.null)
    group_means <- group_means[valid]
    group_names <- group_names[valid]
    valid_idx   <- which(valid)

    if (length(group_means) == 0) {
      showNotification("\u274c Please assign at least one group with replicates.", type = "error")
      return()
    }

    # expr_mat still holds every column of the uploaded matrix, but only the
    # columns assigned to a group take part in the analysis. Restrict it here so
    # the sample-level substrate covers exactly the samples the heatmap is built
    # from, and no others.
    used_cols <- unique(unlist(lapply(valid_idx, function(i) {
      input[[paste0("grp_cols_", i)]]
    })))
    expr_mat <- expr_mat[, used_cols, drop = FALSE]
    
    # Group names become the heatmap's column names, and every downstream
    # consumer matches columns BY NAME: input$ref_condition, col_order_rv(),
    # and the sub-heatmap reorder guard (length(ord_valid) == ncol(sub_mat)).
    # Duplicates make all three ambiguous and fail SILENTLY — the ref drop
    # removes both, the reorder guard quietly skips. Refuse instead of drawing
    # a figure whose columns don't mean what the legend says.
    dup_names <- unique(group_names[duplicated(group_names)])
    if (length(dup_names) > 0) {
        showNotification(paste0("\u274c Duplicate group name(s): ", paste(dup_names, collapse = ", "),
                               ". Group names are used as column IDs \u2014 please make them unique."),
        type = "error", duration = 10)
        return()
      }

    mean_mat           <- do.call(cbind, group_means)
    colnames(mean_mat) <- group_names
    rownames(mean_mat) <- gene_ids

    # ── Top-N most variable genes ────────────────────────────────────────────
    # Applied HERE, before mean_mat_raw is taken and before normalize_matrix(),
    # so the raw and scaled matrices are born already aligned — no extra
    # synchronized drop is needed.
    # The ranking is computed on mean_mat, which is still VST/expression scale.
    # Ranking on the display matrix would be a category error of the same kind
    # drop_zero_variance_rows() guards against: with "Z-score per gene" every
    # row has variance 1, so the "most variable genes" would be an arbitrary
    # tie, and a cosmetic choice would silently decide which genes are shown.
    top_n <- suppressWarnings(as.integer(input$top_var_n %||% 0))
    if (!is.na(top_n) && top_n > 0 && top_n < nrow(mean_mat)) {
      keep_var <- top_variable_rows(mean_mat, top_n)
      showNotification(
        paste0("\u2139\ufe0f Kept the ", length(keep_var), " most variable gene(s) of ",
               nrow(mean_mat), " (ranked by variance across conditions on the ",
               "VST/expression scale)."),
        type = "message", duration = 6)
      mean_mat <- mean_mat[keep_var, , drop = FALSE]
      gene_ids <- gene_ids[keep_var]
      expr_mat <- expr_mat[keep_var, , drop = FALSE]
    }

    # ── Clustering substrate (kept separate from display normalization) ───────
    mean_mat_raw <- mean_mat

    # ── Normalise ──────────────────────────────────────────────────────────────
    if (input$normalization %in% c("zscore_row", "zscore_col", "log2zscore")) {
      # ── Every z-score is computed on the SAMPLE matrix, then averaged ────────
      # The order matters, and not only for zscore_col. Averaging is linear, so
      # the CENTRE is the same whichever comes first, but the SD is not:
      # averaging replicates shrinks variance, so a denominator taken from group
      # means is smaller than one taken from samples. Standardising the group
      # means instead of the samples therefore gives different numbers.
      #
      # Standardising the actual samples is also the more defensible of the two:
      # a per-sample z-score exists to remove per-sample scale differences, and
      # doing it after averaging corrects the averaged profile rather than the
      # samples that differed. Both axes now use the same substrate, so a user
      # switching between them compares transforms of the same data.
      #
      # base::scale() omits NAs column-wise, so a missing value costs only its
      # own cell rather than blanking a whole sample.
      base_mat <- if (identical(input$normalization, "log2zscore") && !vst_applied) {
        log2(expr_mat + 1)
      } else {
        expr_mat
      }
      # zscore_col standardises each SAMPLE across genes; the other two
      # standardise each GENE across samples.
      z_samples <- if (identical(input$normalization, "zscore_col")) {
        scale(base_mat)
      } else {
        t(scale(t(base_mat)))
      }
      mean_mat_scaled <- do.call(cbind, lapply(valid_idx, function(i) {
        sel <- input[[paste0("grp_cols_", i)]]
        rowMeans(z_samples[, sel, drop = FALSE], na.rm = TRUE)
      }))
      colnames(mean_mat_scaled) <- group_names
      rownames(mean_mat_scaled) <- gene_ids
    } else {
    mean_mat_scaled <- tryCatch(
      normalize_matrix(mean_mat,
                       method        = input$normalization,
                       ref_condition = input$ref_condition,
                       already_log_scale = vst_applied),
      
      error = function(e) {
        showNotification(paste("\u274c Normalization error:", e$message), type = "error")
        NULL
      }
    )
    }
    if (is.null(mean_mat_scaled)) return()
    
    # ── Inform user when a requested log transform was a no-op ────────────────
    # When raw counts are uploaded, VST (log-like) is applied upstream and
    # already_log_scale = TRUE, so normalize_matrix() SKIPS the extra log2 to
    # avoid double log-transformation. That decision is silent inside the pure
    # utility, so mirror it here (same condition) so the user understands why a
    # "Log2" choice shows VST values rather than re-logged ones. log2fc is NOT
    # included: there the flag switches ratio->difference, a correct transform,
    # not a skipped step.
    if (vst_applied &&
        input$normalization %in% c("log2", "log2zscore")) {
      msg <- if (identical(input$normalization, "log2")) {
        paste0("These are raw counts, so DESeq2 VST (log-like scale) was already ",
               "applied. The extra Log2 step was skipped to avoid double ",
               "log-transformation \u2014 the heatmap shows VST values.")
      } else {
        paste0("These are raw counts, so DESeq2 VST (log-like scale) was already ",
               "applied. The extra Log2 step was skipped to avoid double ",
               "log-transformation \u2014 the heatmap shows VST values, then ",
               "z-scored per gene.")
      }
      showNotification(paste0("\u2139\ufe0f ", msg), type = "message", duration = 10)
      log_warning(msg)
    }

    # ── Remove NaN / Inf rows (drop_invalid_rows() is unit-tested) ─────────────
    step1           <- drop_invalid_rows(mean_mat_raw, mean_mat_scaled, gene_ids)
    mean_mat_raw    <- step1$raw
    mean_mat_scaled <- step1$scaled
    gene_ids        <- step1$gene_ids
    if (length(step1$dropped) > 0) {
      showNotification(
        paste0("\u26a0\ufe0f ", length(step1$dropped), " gene(s) removed (constant expression): ",
               paste(step1$dropped, collapse = ", ")),
        type = "warning", duration = 10)
    }
    if (nrow(mean_mat_scaled) == 0) {
      showNotification("\u274c No genes remaining after filtering.", type = "error")
      return()
    }

    # ── Drop reference columns for Log2FC (drop_reference_columns() is unit-tested) ──
    # display_col_mask: logical vector over the valid-group columns; TRUE = keep.
    # Reused inside withProgress to align per-column annotation vectors.
    step2 <- drop_reference_columns(mean_mat_raw, mean_mat_scaled,
                                    method        = input$normalization,
                                    ref_condition = input$ref_condition)
    if (!is.null(step2$error)) {
      showNotification(paste("\u274c", step2$error), type = "error")
      return()
    }
    mean_mat_raw      <- step2$raw
    mean_mat_scaled   <- step2$scaled
    display_col_mask  <- step2$display_col_mask
    if (length(step2$dropped) > 0) {
      showNotification(
        paste0("\u2139\ufe0f Reference column(s) hidden: ",
               paste(step2$dropped, collapse = ", ")),
        type = "message", duration = 8)
    }

    # ══ ANALYSIS HALF ═════════════════════════════════════════════════════════
    # Everything from here to analysis_rv() depends only on SIDEBAR inputs and
    # produces the matrices, the dendrograms and the cut. How the figure LOOKS
    # is not decided here -- see the appearance observer below.
    analysis <- tryCatch(
      withProgress(message = "Running analysis...", value = 0, {

        setProgress(0.15, detail = "Checking data size...")
        mem_msgs <- check_memory(length(gene_ids), ncol(expr_mat))
        for (msg in mem_msgs)
          showNotification(msg, type = "warning", duration = 8)

        setProgress(0.30, detail = "Dropping constant genes...")
        # gene_labels moved to the draw step (they depend on map_gene_names,
        # which is live now), so gene_ids is passed as the carrier vector here.
        # drop_zero_variance_rows() only transports that vector alongside the
        # matrices -- it never uses it to decide what to drop.
        step3           <- drop_zero_variance_rows(mean_mat_raw, mean_mat_scaled,
                                                   gene_ids, gene_ids)
        mean_mat_raw    <- step3$raw
        mean_mat_scaled <- step3$scaled
        gene_ids        <- step3$gene_ids
        if (length(step3$dropped) > 0) {
          showNotification(
            paste0("⚠️ ", length(step3$dropped), " gene(s) removed (zero variance)."),
            type = "warning", duration = 8)
        }
        if (nrow(mean_mat_scaled) == 0) {
          showNotification("❌ No genes remaining after filtering.", type = "error")
          return(NULL)
        }

        # The drop_* steps above have pruned gene_ids; re-align the sample-level
        # matrix to match before it is used as the clustering substrate.
        expr_mat <- expr_mat[match(gene_ids, rownames(expr_mat)), , drop = FALSE]

        # ── Clustering substrate: VST -> row z-score -> cluster ────────────────
        # The pipeline is VST, then z-score per gene ACROSS SAMPLES, then
        # clustering, then averaging for display. Standardising here puts every
        # gene on the same footing before any distance is computed, so the tree
        # answers "which genes share a PATTERN" rather than "which genes sit at
        # the same expression LEVEL".
        #
        # Honest note, so nobody flips this back and forth: with the Pearson
        # correlation distance this function uses, the row tree is IDENTICAL
        # either way. Pearson is invariant under a per-row affine map
        # x -> (x - a)/b, and the row z-score is exactly that -- measured on the
        # real data, the two distance matrices differ by 4.4e-16 and the merge
        # tree, the leaf order and cutree() at k = 2/4/6/10 all match exactly.
        # What this line buys is not a different answer but an explicit
        # substrate: if the distance is ever changed to Euclidean, standardising
        # first stops the tree from stratifying by abundance (measured on the
        # same data: without it, Euclidean puts the single most-expressed gene
        # in a cluster of its own).
        #
        # Fallback: a gene that is constant across every sample has sd = 0 and
        # z-scores to NaN, which cluster_rows_or_reason() would refuse outright
        # -- losing the whole dendrogram over one gene. drop_zero_variance_rows()
        # above already removes those, so this cannot normally fire; when it
        # does, the unscaled matrix gives the same tree anyway.
        # Written as three statements rather than one tryCatch with an assignment
        # in its handler: a `<-` inside error = function(e) would land in that
        # handler's own frame, so the flag would never come back.
        .zz <- tryCatch(t(scale(t(expr_mat))), error = function(e) NULL)
        clust_z_applied <- !is.null(.zz) && !anyNA(.zz) && all(is.finite(.zz))
        clust_mat <- if (clust_z_applied) .zz else expr_mat

        setProgress(0.45, detail = "Storing matrices...")
        mat_named           <- mean_mat_scaled
        rownames(mat_named) <- gene_ids
        scaled_matrix_rv(mat_named)
        mat_raw_named           <- clust_mat
        rownames(mat_raw_named) <- gene_ids
        cluster_matrix_rv(mat_raw_named)

        # ── Colour-scale ANCHORS (not the colours themselves) ─────────────────
        # WHERE the ramp is pinned depends on the matrix, so it is decided here.
        # WHICH colours sit on those anchors is style, read on every draw.
        # NULL means "use the zlim slider", i.e. the fully live path.
        col_breaks_fixed <- NULL
        if (identical(input$normalization, "none") ||
            identical(input$normalization, "log2")) {
          mat_min <- as.numeric(quantile(mean_mat_scaled, 0.02, na.rm = TRUE))
          mat_max <- as.numeric(quantile(mean_mat_scaled, 0.98, na.rm = TRUE))
          if (!is.finite(mat_min) || !is.finite(mat_max) || mat_min == mat_max) {
            mat_min <- min(mean_mat_scaled, na.rm = TRUE)
            mat_max <- max(mean_mat_scaled, na.rm = TRUE)
          }
          col_breaks_fixed <- c(mat_min, (mat_min + mat_max) / 2, mat_max)

        } else if (identical(input$normalization, "zscore_col")) {
          # zscore_col standardises each SAMPLE across genes: centred on 0 by
          # construction, but strongly right-skewed (measured on the root data:
          # skewness +1.45 Arabidopsis / +1.00 tomato, median -0.27 / -0.17).
          # On the symmetric +-zlim scale the bulk of the map was pushed into
          # the low colour while a long right tail saturated, so both ends of
          # the scale carried no information. The midpoint stays pinned at 0 --
          # unlike none/log2, zero here is a real reference (the sample mean) --
          # and only the two ARMS become data-driven, independently.
          z_lo <- as.numeric(quantile(mean_mat_scaled, 0.02, na.rm = TRUE))
          z_hi <- as.numeric(quantile(mean_mat_scaled, 0.98, na.rm = TRUE))
          if (!is.finite(z_lo) || !is.finite(z_hi) || z_lo >= 0 || z_hi <= 0) {
            z_span <- max(abs(range(mean_mat_scaled, na.rm = TRUE)))
            if (!is.finite(z_span) || z_span == 0) z_span <- 1
            z_lo <- -z_span; z_hi <- z_span
          }
          col_breaks_fixed <- c(z_lo, 0, z_hi)
          msg <- sprintf(paste0(
            "Column Z-score is centred on 0 but right-skewed, so the colour ",
            "scale was set from the data (%.2f … 0 … %.2f, robust 2%%/98%% ",
            "quantiles) instead of the 'Z-score scale limits' slider, which ",
            "would have saturated most of the map."), z_lo, z_hi)
          showNotification(paste0("ℹ️ ", msg), type = "message", duration = 10)
          log_warning(msg)
        }

        scale_label <- switch(input$normalization,
                              zscore_row = "Row Z-score",
                              zscore_col = "Column Z-score",
                              log2       = "log2(x+1)",
                              log2zscore = "log2 Z-score",
                              log2fc     = "log2FC vs reference",
                              none       = "Expression",
                              "Value")
        scale_label_rv(scale_label)

        setProgress(0.60, detail = "Clustering columns...")
        # Substrate is mean_mat_raw (VST group means), never mean_mat_scaled: a
        # display transform must not change which conditions are judged similar.
        #
        # NOT z-scored, deliberately -- unlike the row substrate above. Centring
        # every gene makes each row sum to ~0 across the conditions, which forces
        # the average correlation BETWEEN columns to -1/(k-1) by construction.
        # Measured on the real data with 4 conditions: the between-condition
        # Spearman correlations collapse from 0.92 .. 0.97 to -0.66 .. +0.25.
        # That is a compositional artefact, not similarity, so the row z-score
        # stops at the row tree and does not travel to this axis.
        col_cluster_arg <- FALSE
        if (isTRUE(input$cluster_cols) && ncol(mean_mat_raw) > 1) {
          hc_col <- hclust(as.dist(1 - cor(mean_mat_raw, method = "spearman")),
                           method = "complete")
          col_cluster_arg <- hc_col
          col_order_rv(colnames(mean_mat_raw)[hc_col$order])
        } else {
          col_order_rv(colnames(mean_mat_scaled))
        }

        setProgress(0.75, detail = "Clustering rows...")
        n_genes          <- nrow(mean_mat_scaled)
        left_anno        <- NULL
        hc               <- NULL
        cluster_rows_arg <- FALSE
        row_split_arg    <- NULL
        cluster_assign_rv(NULL)   # cleared unless a cut actually happens below
        if (isTRUE(input$cluster_rows)) {
          cl <- cluster_rows_or_reason(clust_mat, linkage = input$linkage %||% "complete")
          if (isFALSE(cl$hc)) {
            showNotification(paste0("ℹ️ ", cl$reason,
                                    " Rows are shown unclustered."),
                             type = "warning", duration = 8)
            log_warning(cl$reason)
          } else {
            hc <- cl$hc
          }
        }

        if (!is.null(hc)) {
          if (isTRUE(input$cut_tree)) {
            k         <- max(2, min(input$k_clusters, n_genes))
            cluster_v <- cutree(hc, k = k)
            # Renumber 1..k in DENDROGRAM order so "Cluster 1" is the top block
            # of the figure. cutree() numbers by original row order, which comes
            # out scrambled once row_split + hc impose the dendrogram order.
            unique_in_display <- unique(cluster_v[hc$order])
            cluster_remap     <- setNames(seq_along(unique_in_display),
                                          as.character(unique_in_display))
            cluster_orig      <- as.character(cluster_remap[as.character(cluster_v)])
            pal_clust <- c("#E07B54","#5B9BD5","#6DBF8A","#F4C842",
                           "#A87CC7","#E05C7A","#7B9E87","#C97D6E",
                           "#7D9CC0","#D4A85A")
            clust_col <- setNames(pal_clust[seq_len(k)], as.character(seq_len(k)))
            left_anno <- rowAnnotation(
              Cluster = cluster_orig,
              col     = list(Cluster = clust_col),
              annotation_name_gp      = gpar(fontsize = 11, fontface = "bold"),
              # "bottom", not "top". The column annotation draws its own names
              # ("Temperature", "Genotype") to the LEFT of the annotation bars,
              # which puts them in the same top-left corner this label occupied,
              # and the two overprinted each other. Moving this one under the
              # cluster strip is the only placement where nothing else is drawn.
              annotation_name_side    = "bottom",
              simple_anno_size        = unit(5, "mm"),
              annotation_legend_param = list(Cluster = list(
                title_gp    = gpar(fontsize = 12, fontface = "bold"),
                labels_gp   = gpar(fontsize = 11),
                grid_width  = unit(5, "mm"),
                grid_height = unit(5, "mm"))))
            # row_split must be a single number when cluster_rows is an hclust;
            # passing k lets ComplexHeatmap cut the same tree, so the slices
            # match cluster_v above and the dendrogram is kept.
            row_split_arg    <- k
            cluster_rows_arg <- hc
            cluster_assign_rv(setNames(cluster_orig, gene_ids))
          } else {
            cluster_rows_arg <- hc
          }
        }

        display_order <- if (!is.null(hc)) hc$order else seq_len(n_genes)

        # A brush selection refers to the matrix it was made on. Cleared only
        # when the gene SET changed, so a pure restyle keeps the selection.
        if (!identical(gene_ids, current_gene_ids())) selected_genes(NULL)
        current_gene_ids(gene_ids)

        # DEG highlighting is a SIDEBAR control, so it is frozen here with the
        # rest of the analysis rather than applied live -- keeping the rule
        # "sidebar needs Generate, appearance panel is live" exact.
        is_deg <- normalize_gene_ids(gene_ids) %in% deg_genes_rv()

        prov_caption <- build_provenance_caption(
          scale_label   = scale_label,
          row_clustered = !is.null(hc),
          linkage       = input$linkage %||% "complete",
          col_clustered = !isFALSE(col_cluster_arg),
          n_genes       = n_genes,
          vst_applied   = vst_applied,
          colour_substrate = normalization_substrate(input$normalization),
          # Named per axis, because they are no longer the same matrix.
          row_substrate = paste0(if (vst_applied) "VST " else "",
                                 if (clust_z_applied) "z-scored samples" else "samples"),
          col_substrate = paste0(if (vst_applied) "VST " else "", "group means"))

        # Annotation STRUCTURE only (name, levels, per-displayed-column level).
        # The per-level colours are style and are read fresh on every draw, so
        # they are deliberately not frozen here.
        annot_spec <- lapply(seq_len(n_annots()), function(k) {
          lvls <- trimws(strsplit(input[[paste0("cat_levels_", k)]] %||% "", ",")[[1]])
          lvls <- lvls[nzchar(lvls)]
          vec  <- vapply(valid_idx,
                         function(i) input[[paste0("assign_", i, "_", k)]] %||% "?", "")
          vec  <- vec[display_col_mask]
          list(name = input[[paste0("cat_name_", k)]] %||% paste("Cat", k),
               levels = lvls, vec = vec, k = k)
        })

        list(
          data_digest      = isolate(data_digest()),
          generated_at     = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
          mat_scaled       = mat_named,
          gene_ids         = gene_ids,
          display_order    = display_order,
          is_deg           = is_deg,
          cluster_rows_arg = cluster_rows_arg,
          col_cluster_arg  = col_cluster_arg,
          row_split_arg    = row_split_arg,
          left_anno        = left_anno,
          prov_caption     = prov_caption,
          scale_label      = scale_label,
          col_breaks_fixed = col_breaks_fixed,
          n_genes          = n_genes,
          annot_spec       = annot_spec,
          # The DATA half of the manifest record. The STYLE half is written by
          # the appearance observer onto this same entry.
          config_data = list(
            normalization  = input$normalization,
            ref_condition  = input$ref_condition,
            cluster_rows   = input$cluster_rows, cluster_cols = input$cluster_cols,
            linkage        = input$linkage,
            row_distance   = "pearson",
            row_linkage    = input$linkage,
            col_distance   = "spearman",
            col_linkage    = "complete",
            cut_tree       = input$cut_tree, k = input$k_clusters,
            gene_filter    = input$gene_filter, is_raw_counts = input$is_raw_counts,
            top_var_n      = input$top_var_n,
            vst_applied    = vst_applied,
            scale_label    = scale_label,
            colour_substrate = normalization_substrate(input$normalization),
            n_genes        = nrow(mat_named),
            column_order   = colnames(mat_named),
            groups = lapply(seq_len(n_groups()),
                            function(i) list(name = input[[paste0("grp_name_", i)]],
                                             cols = input[[paste0("grp_cols_", i)]])),
            clusters = {
              ca <- cluster_assign_rv()
              if (is.null(ca)) NULL else as.list(table(factor(as.integer(ca))))
            }
          )
        )
      }),
      error = function(e) {
        showNotification(paste0("❌ Heatmap error: ", conditionMessage(e)),
                         type = "error", duration = 10)
        NULL
      }
    )

    if (is.null(analysis)) return()

    # Fold the setup block away now that a figure exists. After the guard on
    # purpose: a refused or failed run leaves the section open, because that is
    # exactly when the user still needs it.
    setup_collapsed(TRUE)

    # ── History: APPEND only for a new analysis ───────────────────────────────
    # Clicking Generate after changing nothing but colours yields the same
    # digest and therefore no second record. Without this test the bundle would
    # fill with cosmetic near-duplicates -- one matrix CSV and one rerun script
    # per click rather than per analysis. The STYLE half of `config` is written
    # by the appearance observer, onto whichever entry is last.
    local({
      h  <- heatmap_history()
      dg <- analysis$data_digest
      if (length(h) == 0 || !identical(h[[length(h)]]$data_digest, dg)) {
        heatmap_history(c(h, list(list(
          generated_at = analysis$generated_at,
          species      = ACTIVE_SPECIES,
          data_digest  = dg,
          config       = analysis$config_data,
          matrix       = analysis$mat_scaled))))
      }
    })

    analysis_rv(analysis)
  })

  # ══ APPEARANCE HALF ════════════════════════════════════════════════════════
  # Redraws whenever the analysis changes OR any appearance control moves. No
  # recompute: it reads the frozen matrices out of analysis_rv().
  #
  # Debounced, because the font-size and z-limit sliders fire continuously
  # while being dragged and every redraw re-renders the interactive widget.
  style_spec_raw <- reactive({
    a <- analysis_rv()
    list(
      col_low = input$col_low, col_mid = input$col_mid, col_high = input$col_high,
      zlim            = input$zlim,
      show_labels     = isTRUE(input$show_labels),
      fontsize_labels = input$fontsize_labels,
      anno_bar_size   = input$anno_bar_size,
      map_gene_names  = isTRUE(input$map_gene_names),
      # Reading the level colours here is what makes this reactive DEPEND on
      # them. Their values are re-read through get_annot_colors() inside the
      # draw, so the dependency and the value can never disagree.
      cat_cols = if (is.null(a)) NULL else lapply(a$annot_spec, function(s)
        vapply(seq_along(s$levels),
               function(j) input[[paste0("cat_col_", s$k, "_", j)]] %||% "", ""))
    )
  })
  style_spec <- debounce(style_spec_raw, 450)

  observe({
    a <- analysis_rv()
    style_spec()          # dependency ONLY -- it decides WHEN to redraw
    req(!is.null(a))

    # Everything below is isolated: this observer WRITES heatmap_history,
    # col_fun_rv, gene_labels_rv and main_ht_rv, and a reactive read of any of
    # them would make it re-trigger itself.
    isolate(tryCatch({
      # Values are re-read FRESH here rather than taken from the debounced
      # snapshot. debounce() emits its first value at app START -- before any
      # input has arrived and before analysis_rv() exists -- so drawing from
      # that snapshot meant the first redraw after Generate painted with an
      # EMPTY style: zlim was NULL, c(NULL[1], 0, NULL[2]) collapsed to a
      # length-1 vector, and colorRamp2() threw "Length of breaks should be
      # equal to length of colors". The debounce still does its job -- it
      # decides WHEN this observer runs -- it just no longer decides WITH WHAT.
      s <- style_spec_raw()

      # Colours placed on the anchors the analysis chose, or on the zlim slider
      # when it left them open. heatmap_color_breaks() (utils/color_utils.R) owns
      # every guard: a cleared (NA) slider, reversed or identical handles, and --
      # the one that used to produce a wrong figure silently -- a window that
      # does not contain 0, where hardcoding the midpoint to 0 gave colorRamp2()
      # non-increasing breaks and a scale that ran low -> BELOW low -> high.
      # The same rule is emitted into rerun_heatmap_NN.R, so the bundle keeps
      # reproducing what was on screen.
      brk     <- heatmap_color_breaks(s$zlim, fixed = a$col_breaks_fixed)
      col_fun <- colorRamp2(brk$breaks, c(s$col_low %||% "#00008b",
                                          s$col_mid %||% "white",
                                          s$col_high %||% "#FFFF00"))
      col_fun_rv(col_fun)

      # Shown with a fixed id so dragging the slider replaces the message rather
      # than stacking one per redraw, and removed as soon as the setting is
      # ordinary again. Logged too, so it travels in the session log.
      if (!is.null(brk$note)) {
        showNotification(paste0("ℹ️ ", brk$note), type = "warning",
                         duration = 12, id = "zlim_note")
        log_warning(brk$note)
      } else {
        removeNotification("zlim_note")
      }

      gene_ids    <- a$gene_ids
      gene_labels <- gene_ids
      if (isTRUE(s$map_gene_names)) {
        symbols <- map_gene_symbols(gene_ids)
        if (!is.null(symbols)) {
          norm_ids    <- normalize_gene_ids(gene_ids)
          gene_labels <- ifelse(is.na(symbols[norm_ids]), gene_ids, symbols[norm_ids])
        }
      }
      gene_labels_rv(setNames(gene_labels[a$display_order],
                              normalize_gene_ids(gene_ids[a$display_order])))

      # Top annotation: structure frozen by the analysis, colours read live.
      annot_list <- list(); color_list <- list(); legend_list <- list()
      for (sp in a$annot_spec) {
        if (length(sp$levels) == 0 || length(sp$vec) == 0) next
        annot_list[[sp$name]]  <- sp$vec
        color_list[[sp$name]]  <- get_annot_colors(sp$levels, sp$k)
        legend_list[[sp$name]] <- list(
          title_gp    = gpar(fontsize = 14, fontface = "bold"),
          labels_gp   = gpar(fontsize = 12),
          grid_width  = unit(5, "mm"),
          grid_height = unit(5, "mm"))
      }
      col_anno <- if (length(annot_list) == 0) NULL else do.call(
        HeatmapAnnotation,
        c(annot_list,
          list(col                     = color_list,
               annotation_name_gp      = gpar(fontsize = 12, fontface = "bold"),
               annotation_name_side    = "left",
               simple_anno_size        = unit(s$anno_bar_size %||% 5, "mm"),
               gap                     = unit(1, "mm"),
               annotation_legend_param = legend_list)))

      is_deg  <- a$is_deg
      fg_col  <- ifelse(is_deg, "#FF00FF", "black")
      fg_face <- ifelse(is_deg, "bold", "plain")
      fs      <- s$fontsize_labels %||% 14

      gene_name_anno <- if (isTRUE(s$show_labels)) {
        if (isTRUE(s$map_gene_names)) {
          gene_id_grey <- ifelse(gene_labels == gene_ids, "", gene_ids)
          rowAnnotation(
            space   = anno_empty(border = FALSE, width = unit(3, "mm")),
            symbol  = anno_text(gene_labels,
                                gp = gpar(fontsize = fs, fontface = fg_face, col = fg_col),
                                location = 0, just = "left"),
            space2  = anno_empty(border = FALSE, width = unit(3, "mm")),
            gene_id = anno_text(gene_id_grey,
                                gp = gpar(fontsize = fs * 0.85, fontface = fg_face,
                                          col = ifelse(is_deg, "#FF00FF", "grey50")),
                                location = 0, just = "left"),
            gap = unit(2, "mm"))
        } else {
          rowAnnotation(
            space   = anno_empty(border = FALSE, width = unit(3, "mm")),
            gene_id = anno_text(gene_ids,
                                gp = gpar(fontsize = fs, fontface = fg_face, col = fg_col),
                                location = 0, just = "left"),
            gap = unit(2, "mm"))
        }
      } else if (any(is_deg)) {
        # Labels hidden but DEGs flagged: keep a minimal colour strip so the
        # highlight is not silently lost.
        rowAnnotation(
          deg_mark = anno_simple(as.numeric(is_deg),
                                 col = c("0" = "white", "1" = "#FF00FF"),
                                 gp  = gpar(col = "grey70")),
          annotation_name_gp = gpar(fontsize = 9, fontface = "bold"),
          width = unit(3, "mm"))
      } else NULL

      ht <- tryCatch(
        Heatmap(a$mat_scaled,
                name                 = "zscore",
                col                  = col_fun,
                cluster_rows         = a$cluster_rows_arg,
                cluster_columns      = a$col_cluster_arg,
                show_row_names       = FALSE,
                show_column_names    = TRUE,
                column_names_gp      = gpar(fontsize = 14),
                column_title         = a$prov_caption,
                column_title_gp      = gpar(fontsize = 11, fontface = "plain"),
                top_annotation       = col_anno,
                left_annotation      = a$left_anno,
                right_annotation     = gene_name_anno,
                row_split            = a$row_split_arg,
                heatmap_legend_param = list(
                  title       = a$scale_label,
                  title_gp    = gpar(fontsize = 14, fontface = "bold"),
                  labels_gp   = gpar(fontsize = 12),
                  grid_width  = unit(5, "mm"),
                  grid_height = unit(5, "mm"))),
        error = function(e) {
          showNotification(paste0("❌ Heatmap error: ", conditionMessage(e)),
                           type = "error", duration = 10)
          NULL
        })

      if (!is.null(ht)) {
        main_ht_rv(ht)

        # ── History, strategy (b): update the LAST entry IN PLACE ────────────
        # Same analysis, different paint -- so the record is refreshed, never
        # duplicated. restyled_at sits alongside generated_at so the bundle
        # still states when the numbers were computed AND when the figure last
        # changed appearance, with no lost audit trail.
        h <- heatmap_history()
        if (length(h) > 0) {
          # Plain element-wise overwrite, NOT modifyList(). modifyList()
          # RECURSES into a value that is itself a list and merges it by name --
          # and `annotation_categories` is an UNNAMED list (one element per
          # factor), so the recursion matches nothing and silently leaves
          # whatever was written on the first draw. Every level colour then
          # stayed frozen at its palette default no matter what the user picked,
          # while `colors` (a named list) and the atomic fields updated
          # correctly, which made it look like the mechanism worked. The style
          # half is meant to replace the style half wholesale, so assign it.
          style_cfg <- list(
            # isTRUE(): the control is hidden for annotations with no symbol
            # column, so the input can be absent. NULL would serialise as
            # `null` and read as "unknown" rather than "not applied".
            map_gene_names = isTRUE(s$map_gene_names),
            zlim           = s$zlim,
            colors         = list(low = s$col_low, mid = s$col_mid, high = s$col_high),
            anno_bar_size  = s$anno_bar_size,
            annotation_categories = lapply(a$annot_spec, function(sp) list(
              name        = sp$name,
              levels      = sp$levels,
              colors      = if (length(sp$levels) > 0)
                              as.list(get_annot_colors(sp$levels, sp$k)) else NULL,
              assignments = if (length(sp$vec) == ncol(a$mat_scaled))
                              as.list(stats::setNames(sp$vec, colnames(a$mat_scaled)))
                            else as.list(sp$vec)))
          )
          cfg <- h[[length(h)]]$config
          for (nm in names(style_cfg)) cfg[[nm]] <- style_cfg[[nm]]
          h[[length(h)]]$config      <- cfg
          h[[length(h)]]$restyled_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
          heatmap_history(h)
        }

        makeInteractiveComplexHeatmap(
          input, output, session, ht, "ht_shiny",
          brush_action = function(df, input, output, session) {
            ids <- current_gene_ids()
            if (!is.null(df) && nrow(df) > 0 && !is.null(ids)) {
              rows <- unique(unlist(df$row_index))
              rows <- rows[rows >= 1 & rows <= length(ids)]
              if (length(rows) > 0) selected_genes(ids[rows])
            }
          })
      }
    },
    # A failure in here must not take the session down with it. An uncaught
    # error inside an observer aborts the whole reactive flush, so a single bad
    # redraw would silently break every unrelated panel queued behind it -- the
    # RNA-seq tab, the enrichment modules, the species switch. Caught, reported,
    # and the rest of the app carries on.
    error = function(e) {
      showNotification(paste0("❌ Heatmap render error: ", conditionMessage(e)),
                       type = "error", duration = 10)
      NULL
    }))
  })

  output$enrichment_panel <- renderUI({
    req(selected_genes())
    n <- length(selected_genes())
    tagList(
      h4(paste0("\U0001f52c Enrichment Analysis \u2014 ", n, " selected genes")),
      if (!identical(input$universe_mode %||% "default", "default") && !is.null(universe_genes())) {
        tags$p(
          style = "color: #2c7bb6; font-size: 0.9em;",
          icon("circle-info"),
          paste0(
            if (identical(input$universe_mode, "custom")) " Custom pasted universe active: "
            else " Expression-threshold universe active: ",
            length(universe_genes()), " genes"
          )
        )
      },
      fluidRow(
        column(3,
               goControlsUI("go"),
               hr(),
               keggControlsUI("kegg")
        ),
        column(9,
               tabsetPanel(
                 tabPanel("GO Table",          goTableUI("go")),
                 tabPanel("GO Plot",           goPlotUI("go")),
                 tabPanel("GO Heatmap",        goHeatmapUI("go")),
                 tabPanel("KEGG Table",        keggTableUI("kegg")),
                 tabPanel("KEGG Plot",         keggPlotUI("kegg")),
                 tabPanel("Selection History", selectionHistoryUI("hist")),
                 tabPanel("STRING Network",    stringUI("string")),
                 tabPanel("Selected Genes",    tairGenesUI("tair")),
                 tabPanel("Gene Detail",       geneDetailUI("gene_detail"))
               )
        )
      )
    )
  })

  # ── Module servers ─────────────────────────────────────────────────────────
  go_cart_signal <- reactiveVal(NULL)
  
  go_handles <- goServer("go",
                         selected_genes   = selected_genes,
                         scaled_matrix_rv = scaled_matrix_rv,
                         # Substrate + linkage, so the GO sub-heatmap clusters
                         # exactly like the main heatmap instead of falling back
                         # to ComplexHeatmap's euclidean-on-display default.
                         cluster_matrix_rv = cluster_matrix_rv,
                         linkage_rv       = reactive(input$linkage %||% "complete"),
                         scale_label_rv = scale_label_rv,
                         gene_labels_rv   = gene_labels_rv,
                         universe         = universe_genes,
                         on_new_heatmap   = go_cart_signal,
                         col_fun_rv       = col_fun_rv,
                         col_order_rv     = col_order_rv,
                         deg_genes_rv     = deg_genes_rv,
                         context = "heatmap"
  )

  observeEvent(go_cart_signal(), {
    sig <- go_cart_signal()
    if (!is.null(sig) && length(sig$genes) > 0) {
      cart_heatmap_genes(sig$genes)
      showNotification(
        paste0("\U0001f525 Custom heatmap generated with ", length(sig$genes),
               " genes. Scroll down to see it."),
        type = "message", duration = 5
      )
      shinyjs::runjs("
        setTimeout(function() {
          var el = document.getElementById('cart_heatmap_panel');
          if (el) { el.scrollIntoView({ behavior: 'smooth', block: 'start' }); }
        }, 300);
      ")
    }
  }, ignoreNULL = TRUE, ignoreInit = TRUE)

  kegg_handles <- keggServer("kegg",
             selected_genes   = selected_genes,
             scaled_matrix_rv = scaled_matrix_rv,
             deg_genes_rv     = deg_genes_rv,
             # Same reactive the GO panel already uses, so GO and KEGG can no
             # longer disagree about what the experiment tested.
             universe         = universe_genes)
  
  gsea_de_handles <- gseaServer("gsea_de", de_results_r = rnaseq_de$de_results_all)

  hist_out <- selectionHistoryServer("hist", current_genes = selected_genes)
  observeEvent(hist_out$loaded_genes(), {
    genes <- hist_out$loaded_genes()
    if (!is.null(genes) && length(genes) > 0)
      selected_genes(genes)
  })

  stringServer("string", selected_genes = selected_genes)

  tairGenesServer("tair",
                  selected_genes = selected_genes,
                  gene_labels_rv = gene_labels_rv)

  geneDetailServer("gene_detail",
                   selected_genes   = selected_genes,
                   scaled_matrix_rv = scaled_matrix_rv,
                   gene_labels_rv   = gene_labels_rv,
                   col_fun_rv       = col_fun_rv,
                   col_order_rv     = col_order_rv,
                   deg_genes_rv     = deg_genes_rv
  )

  # preview_matrix_for_optk is the RAW/VST preview matrix, not a display matrix,
  # which is why the argument is named cluster_matrix_rv. linkage_rv is passed
  # so the suggested k always scores the tree the heatmap will actually cut.
  optimalKServer("optimal_k",
                 cluster_matrix_rv = preview_matrix_for_optk,
                 linkage_rv        = reactive(input$linkage %||% "complete"))
  # connection to the three modules of the application (RNAseq)
  # Same start_step as rnaseqDdsUI() above: the design panel is rendered
  # server-side, so it needs the offset too.
  rnaseq_handles <- rnaseqDdsServer("rnaseq_dds", start_step = 2)
  
  appLogServer("applog")
  
  repro_collectors <- list(
    # EVERY heatmap generated this session, each frozen at generation time.
    # Was a single list built from live input$ values when the download button
    # was pressed, which both drifted from the plotted matrix and kept only the
    # last figure. See the freeze block in the go_btn observer.
    heatmap_config = function() {
      h <- tryCatch(heatmap_history(), error = function(e) NULL)
      if (is.null(h) || length(h) == 0) return(NULL)
      lapply(h, function(e) c(list(generated_at = e$generated_at,
                                   species      = e$species), e$config))
    },

    # One {matrix, config} pair per generated heatmap, so the bundle can write a
    # CSV and a rerun script for each. `config` is the SAME frozen object the
    # manifest reports, so the colour-scale branch in rerun_heatmap.R can no
    # longer disagree with the matrix sitting next to it: normalization decides
    # the anchors ("none"/"log2" get a 2%/98% quantile scale, everything else
    # stays centred on 0 with zlim), and both now come from the same snapshot.
    heatmap_matrix = function() {
      h <- tryCatch(heatmap_history(), error = function(e) NULL)
      if (is.null(h) || length(h) == 0) return(NULL)
      lapply(h, function(e) list(matrix = e$matrix, config = e$config))
    },
    gsea_params = function() {
      tryCatch(gsea_de_handles$gsea_params(), error = function(e) NULL)
    },
    
    gsea_history = function() {
      tryCatch(gsea_de_handles$history(), error = function(e) NULL)
    },
    
    de_params = function() list(
      design    = input[["rnaseq_dds-rnaseq_design_formula"]],
      mode      = input[["rnaseq_de-result_mode"]],
      factor    = input[["rnaseq_de-contrast_factor"]],
      test      = input[["rnaseq_de-contrast_num"]],
      reference = input[["rnaseq_de-contrast_den"]],
      alpha     = input[["rnaseq_de-alpha_cutoff"]]
    ),
    
    de_history = function() {
      tryCatch(de_comp_log$comparisons(), error = function(e) NULL)
    },
    
    gene_lists = function() {
      go_h <- tryCatch(go_handles$result(),    error = function(e) NULL)
      go_d <- tryCatch(go_de_handles$result(), error = function(e) NULL)
      # The partition itself, one gene_lists/cluster_N.txt per cluster.
      # heatmap_config records `k`, but not which genes k produced; without
      # these files the bundle would describe the cut without containing its
      # result. Emitted as plain character vectors, the shape repro_module.R
      # already writes for every other gene list.
      ca <- cluster_assign_rv()
      cluster_lists <- if (!is.null(ca)) {
        setNames(lapply(sort(unique(as.integer(ca))),
                        function(i) names(ca)[as.integer(ca) == i]),
                 paste0("cluster_", sort(unique(as.integer(ca)))))
      } else list()
      c(list(
        selected            = selected_genes(),
        degs                = deg_genes_rv(),
        cart                = cart_heatmap_genes(),
        de_sig              = tryCatch(rnaseq_de$de_sig_genes(), error = function(e) NULL),
        go_heatmap_gene_ids = go_h$params$gene_ids,
        go_heatmap_universe = go_h$params$universe,
        go_de_gene_ids      = go_d$params$gene_ids,
        go_de_universe      = go_d$params$universe
        # gsea_1_ranked_ids / gsea_1_ranked_scores removed. They came from
        # gsea_de_handles$ranked(), which is reactive(rank_vec()) -- a LIVE
        # vector that rebuilds whenever the DE result or the "Rank genes by"
        # dropdown changes. So the file named "gsea_1" held whatever the sidebar
        # said at DOWNLOAD time, not what any run used, and could describe a run
        # that never happened. write_rerun_GSEA() already writes
        # gene_lists/gsea_run_NN_ranked_ids.txt + _scores.txt from gsea_history,
        # frozen at run time, one pair per actual run -- which is what
        # rerun_GSEA.R reads. This pair was the only thing left in the bundle
        # that was not frozen when it happened.
      ), cluster_lists)
    },
    go_params = function() {
      extract <- function(h, label) {
        r <- tryCatch(h$result(), error = function(e) NULL)
        if (is.null(r) || is.null(r$params)) return(NULL)
        p <- r$params
        list(
          context         = label,
          ontology        = p$process,
          padj_method     = p$padj_method,
          pval_cutoff     = p$pval_cutoff,
          qval_cutoff     = p$qval_cutoff,
          simplify_cutoff = p$simplify_cutoff,
          species         = p$species,
          org_db_pkg      = p$org_db_pkg,
          key_type        = p$key_type,
          n_genes_input   = length(p$gene_ids),
          n_universe      = if (is.null(p$universe)) NA else length(p$universe)
        )
      }
      list(
        heatmap = extract(go_handles,    "go (heatmap panel)"),
        de      = extract(go_de_handles, "go_de (DE panel)")
      )
    },
    
    # Merged in the order the runs HAPPENED, not panel by panel. Concatenating
    # go_handles then go_de_handles listed every heatmap-panel run before every
    # DEG-panel run, so "run 2" in the manifest could be the third thing you did.
    # Each entry now carries a run_at stamp (see go_module.R), so the merge can
    # sort. seq_along() is the tie-break, which keeps insertion order for two
    # runs that landed in the same millisecond.
    #
    # Deliberately only these two panels: GO runs made from the standalone
    # "GO Analysis" feature are NOT part of the bundle, by design. The bundle's
    # README says so, so the omission is stated rather than silent.
    go_history = function() {
      h  <- tryCatch(go_handles$history(),    error = function(e) NULL)
      hd <- tryCatch(go_de_handles$history(), error = function(e) NULL)
      all <- c(h %||% list(), hd %||% list())
      if (length(all) == 0) return(all)
      all[order(vapply(all, function(e) e$run_at %||% "", ""), seq_along(all))]
    },
    
    # EVERY KEGG run this session, both panels, in the order they were run.
    # Was a snapshot of live widgets taken at download time, two of whose IDs
    # ("kegg-pval_cutoff"/"kegg-qval_cutoff") do not exist -- the widgets are
    # named pval/qval -- so those two fields serialised as null on every bundle
    # ever produced, and n_genes_input read the heatmap brush even for the DE
    # panel. Emitted under the existing `kegg_params` key rather than a new
    # `kegg_history` one because build_manifest() writes a fixed set of keys and
    # would silently drop an unrecognised collector name.
    # Sorted by run_at for the same reason as go_history above.
    kegg_params = function() {
      h  <- tryCatch(kegg_handles$history(),    error = function(e) NULL)
      hd <- tryCatch(kegg_de_handles$history(), error = function(e) NULL)
      all <- c(h %||% list(), hd %||% list())
      if (length(all) == 0) return(all)
      all[order(vapply(all, function(e) e$run_at %||% "", ""), seq_along(all))]
    },

    annotation_info = function() c(
      species        = ACTIVE_SPECIES,
      itag_version   = if (identical(ACTIVE_SPECIES, "tomato")) ITAG_VERSION %||% NA else NA,
      org_db         = ORG_DB_PKG,
      org_db_version = tryCatch(as.character(utils::packageVersion(ORG_DB_PKG)), error = function(e) NA),
      kegg_code      = KEGG_CODE,
      # ── ID-map provenance ───────────────────────────────────────────────────
      # For tomato this file is the sole authority for KEGG: every tomato KEGG
      # p-value is a function of it. Recording only the path would not be
      # enough -- two revisions of the same filename give different results and
      # otherwise identical bundles -- so the content is fingerprinted too.
      # digest is already a dependency (global.R:16).
      id_map_source  = if (is.na(ID_MAP_SOURCE)) NA_character_ else ID_MAP_SOURCE,
      id_map_md5     = tryCatch(
        if (is.na(ID_MAP_SOURCE)) NA_character_
        else digest::digest(file = ID_MAP_SOURCE, algo = "md5"),
        error = function(e) NA_character_),
      id_map_entrez_pairs = tryCatch(
        if (exists("ENTREZ_MAP") && !is.null(ENTREZ_MAP)) nrow(ENTREZ_MAP) else NA_integer_,
        error = function(e) NA_integer_),
      DESeq2         = tryCatch(as.character(utils::packageVersion("DESeq2")), error = function(e) NA)
    )
  )
  
  reproServer("repro",        collectors = repro_collectors)
  reproServer("repro_rnaseq", collectors = repro_collectors)

  appLogServer("applog_rnaseq")   # same underlying log as "applog" — just a second button
  
  demoServer(input, output, session, on_load = function(df) demo_df(df))
  
  # Universe for DEG enrichment: all gene IDs present in the DESeq2 object.
  # Returns NULL when no dds is available yet, which is safe for run_go_enrichment().
  universe_dds <- reactive({
    dds <- rnaseq_handles$dds()
    if (is.null(dds)) return(NULL)
    
    ids <- trimws(as.character(rownames(dds)))
    
    if (identical(ACTIVE_SPECIES, "tomato")) {
      ids <- canonicalize_tomato_ids(ids)
    }
    
    unique(ids)
  })
  
  # VST-normalised, row z-scored matrix derived from the DESeq2 object.
  # Used exclusively by the GO Heatmap tab in the DEG enrichment panel.
  # blind = TRUE: normalisation is for visualisation only,
  # not for differential expression testing.
  de_scaled_matrix <- reactive({
    dds <- rnaseq_handles$dds()
    req(!is.null(dds))
    
    tryCatch({
      # Same nsub guards as the other two call sites. This dds IS pre-filtered,
      # so it usually clears them, but a small or shallow experiment still can
      # not: fall back to the full-trend estimator rather than error out.
      vsd <- tryCatch(
        DESeq2::vst(dds, blind = TRUE, nsub = min(1000, nrow(dds))),
        error = function(e) DESeq2::varianceStabilizingTransformation(dds, blind = TRUE)
      )
      mat <- as.matrix(SummarizedExperiment::assay(vsd))
      
      rn <- trimws(as.character(rownames(mat)))
      if (identical(ACTIVE_SPECIES, "tomato")) {
        rn <- canonicalize_tomato_ids(rn)
      }
      rownames(mat) <- rn
      
      t(scale(t(mat)))
    }, error = function(e) NULL)
  })
  
  # Gene labels for DEG enrichment GO heatmap.
  # A simple identity mapping is sufficient in the RNA-seq context:
  # displayed label = gene ID.
  
  de_gene_labels <- reactive({
    dds <- rnaseq_handles$dds()
    if (is.null(dds)) return(setNames(character(0), character(0)))
    
    ids <- trimws(as.character(rownames(dds)))
    
    if (identical(ACTIVE_SPECIES, "tomato")) {
      ids <- canonicalize_tomato_ids(ids)
    }
    
    norm_ids <- normalize_gene_ids(ids)
    stats::setNames(ids, norm_ids)
  })
  
  rnaseqVizServer(
    "rnaseq_viz",
    dds_rv = rnaseq_handles$dds,
    metadata_rv = rnaseq_handles$metadata_saved
  )
  
  rnaseqEdaServer("rnaseq_eda",
                  dds_rv      = rnaseq_handles$dds,
                  metadata_rv = rnaseq_handles$metadata_saved)
  
  rnaseq_de <- rnaseqDeServer(
    "rnaseq_de",
    dds_rv = rnaseq_handles$dds,
    metadata_rv = rnaseq_handles$metadata_saved,
    # Passed so the staleness guard in de_sig_genes() depends on a REACTIVE
    # rather than on the ACTIVE_SPECIES global, which invalidates nothing.
    # Without this the DEG enrichment panel kept offering the previous
    # organism's gene IDs after "Switch species & reload".
    species_rv = species_rv
  )
  
  rnaseq_session_info <- reactive({
    c(
      species        = ACTIVE_SPECIES,
      species_common = SPECIES_COMMON,
      itag_version   = if (identical(ACTIVE_SPECIES, "tomato")) ITAG_VERSION %||% NA else NA,
      org_db         = ORG_DB_PKG,
      org_db_version = tryCatch(as.character(utils::packageVersion(ORG_DB_PKG)), error = function(e) NA),
      key_type       = KEY_TYPE,
      kegg_code      = KEGG_CODE,
      DESeq2         = tryCatch(as.character(utils::packageVersion("DESeq2")), error = function(e) NA),
      ashr           = tryCatch(as.character(utils::packageVersion("ashr")), error = function(e) NA)
    )
  })
  
  de_comp_log <- deComparisonLogServer(
    "de_comp_log",
    run_trigger     = rnaseq_de$run_trigger,
    de_results_all  = rnaseq_de$de_results_all,
    contrast_params = rnaseq_de$contrast_params,
    session_info    = rnaseq_session_info
  )
  
  # ── RNA-seq main panel state ───────────────────────────────────────────────
  # Allowed values: "table", "enrich", "kegg", "save"
  # Venn remains modal-based, so it does not need a persistent panel state.
  rnaseq_panel_rv <- reactiveVal("table")
  
  # DE gene list history: reuses the same selectionHistoryServer with a separate
  # instance ("de_hist") so heatmap and DE histories are independent.
  de_hist_out <- selectionHistoryServer(
    "de_hist",
    current_genes = rnaseq_de$de_sig_genes
  )
  
  venn_mod <- vennServer("venn", de_history_rv = de_hist_out$history)
  # Independent manual Venn module — no connection to rnaseq flow
  vennManualServer("venn_manual")
  
  # ── RNA-seq mutually exclusive panel navigation ────────────────────────────
  observeEvent(input$open_de_table, {
    rnaseq_panel_rv("table")
  })
  
  observeEvent(input$open_go_de, {
    req(rnaseq_de$de_sig_genes())
    rnaseq_panel_rv("enrich")
  })
  
  observeEvent(input$open_gsea_de, {
    rnaseq_panel_rv("gsea")
  })
  
  observeEvent(input$open_kegg_de, {
    req(rnaseq_de$de_sig_genes())
    rnaseq_panel_rv("kegg")
  })
  
  observeEvent(input$open_de_hist, {
    rnaseq_panel_rv("save")
  })
  
  observeEvent(input$open_venn, {
    venn_mod$show()
  })
  
  # ── RNA-seq central panel renderer ─────────────────────────────────────────
  output$rnaseq_main_panel <- renderUI({
    active_panel <- rnaseq_panel_rv()
    
    switch(
      active_panel,
      
      "table" = tagList(
        DT::DTOutput("rnaseq_de-de_table")
      ),
      
      "enrich" = {
        req(rnaseq_de$de_sig_genes())
        n <- length(rnaseq_de$de_sig_genes())
        
        tagList(
          hr(),
          h4(paste0("🧬 GO Enrichment Analysis on DEGs — ", n, " genes")),
          fluidRow(
            column(
              3,
              goControlsUI("go_de")
            ),
            column(
              9,
              tabsetPanel(
                tabPanel("GO Table", goTableUI("go_de")),
                tabPanel("GO Plot", goPlotUI("go_de")),
                tabPanel("GO Heatmap", goHeatmapUI("go_de"))
              )
            )
          )
        )
      },
      
      "gsea" = {
        req(rnaseq_de$de_results_all())
        tagList(
          hr(),
          h4("\U0001f4c8 GSEA \u2014 rank-based enrichment on the full DE result"),
          fluidRow(
            column(3, gseaControlsUI("gsea_de")),
            column(9, gseaOutputUI("gsea_de"))
          )
        )
      },
      
      "kegg" = {
        req(rnaseq_de$de_sig_genes())
        n <- length(rnaseq_de$de_sig_genes())
        
        tagList(
          hr(),
          h4(paste0("🗺️ KEGG Enrichment Analysis on DEGs — ", n, " genes")),
          fluidRow(
            column(
              3,
              keggControlsUI("kegg_de")
            ),
            column(
              9,
              tabsetPanel(
                tabPanel("KEGG Table", keggTableUI("kegg_de")),
                tabPanel("KEGG Plot", keggPlotUI("kegg_de"))
              )
            )
          )
        )
      },
      
      "save" = tagList(
        hr(),
        h4("🗂️ DE gene list history"),
        selectionHistoryUI("de_hist", enable_load = FALSE),
        hr(),
        # deComparisonLogServer("de_comp_log") has been running since app start
        # and its `comparisons` reactive already feeds the bundle's de_history,
        # but its UI half was never mounted anywhere -- so the per-comparison
        # table, the delete/clear buttons and the three session exports
        # (R script / JSON / CSV zip) were unreachable. Mounted here, next to
        # the other per-session RNA-seq records it belongs with.
        deComparisonLogUI("de_comp_log"),
        hr(),
        h4("\U0001f4c1 Reproducibility bundle (RNAseq)"),
        reproDownloadUI("repro_rnaseq")
      ),
      
      tagList(
        DT::DTOutput("rnaseq_de-de_table")
      )
    )
  })
  
  # ── Module servers for DEG enrichment ──────────────────────────────────────
  # scaled_matrix_rv is not passed here: DEG enrichment operates on gene IDs
  # only, not on expression values, so the heatmap matrix is not needed.
  # deg_genes_rv is also not passed: in this context every gene IS a DEG.
  go_de_handles <- goServer("go_de",
                            selected_genes = rnaseq_de$de_sig_genes,
                            universe       = universe_dds,      # all genes in the dds, not just the heatmap's
                            scaled_matrix_rv = de_scaled_matrix,
                            gene_labels_rv   = de_gene_labels,
                            context = "rnaseq",
                            scale_label_rv = reactive("Row Z-score")
  )
  
  kegg_de_handles <- keggServer("kegg_de",
             selected_genes = rnaseq_de$de_sig_genes,
             # rownames(dds) after the count filter — identical to the universe
             # goServer("go_de") is already given.
             universe       = universe_dds
  )
  
  rnaseqVolcanoServer(
    "rnaseq_volcano",
    de_results_r = rnaseq_de$de_results_all
  )

  # ── Standalone GO Analysis section (feature = "go") ─────────────────────────
  # Independent gene-list entry point into the same go_module.r machinery used
  # by the heatmap-brush ("go") and DEG ("go_de") panels above. scaled_matrix_rv
  # etc. are passed so its "GO Heatmap" sub-tab can show real expression data
  # whenever a main heatmap already exists; heatmap_selected_genes/de_sig_genes_rv
  # are optional one-click bridges so users don't have to retype gene lists that
  # already exist elsewhere in the app.
  goSectionServer(
    "go_section",
    scaled_matrix_rv       = scaled_matrix_rv,
    cluster_matrix_rv      = cluster_matrix_rv,
    linkage_rv             = reactive(input$linkage %||% "complete"),
    gene_labels_rv         = gene_labels_rv,
    col_fun_rv             = col_fun_rv,
    col_order_rv           = col_order_rv,
    deg_genes_rv           = deg_genes_rv,
    heatmap_selected_genes = selected_genes,
    de_sig_genes_rv        = rnaseq_de$de_sig_genes,
    # Lets the panel refresh its gene-ID placeholders after a species switch.
    species_rv             = species_rv
  )

  # ── Main heatmap PDF export ─────────────────────────────────────────────────
  # Vector export of the main figure at a user-controlled size, for figures
  # going into a paper. main_ht_rv holds the last built Heatmap object, so the
  # export re-draws exactly what is on screen rather than rebuilding it.
  output$dl_main_ht_pdf <- downloadHandler(
    filename = function() paste0("heatmap_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf"),
    content  = function(f) {
      ht <- main_ht_rv()
      if (is.null(ht)) {
        showNotification("\u274c Generate a heatmap first.", type = "error")
        return()
      }
      w <- input$main_ht_w %||% 10
      h <- input$main_ht_h %||% 12
      pdf(f, width = as.numeric(w), height = as.numeric(h))
      on.exit(dev.off(), add = TRUE)
      draw(ht)
    }
  )

  # ── Cluster tools ───────────────────────────────────────────────────────────
  # Per-cluster enrichment hand-off, membership export and stability estimation.
  # cluster_assign_rv is NULL unless the last run actually cut the tree, so this
  # whole panel simply does not render otherwise.
  cluster_stability_rv <- reactiveVal(NULL)

  # Recomputing membership invalidates any previous stability estimate.
  observeEvent(cluster_assign_rv(), { cluster_stability_rv(NULL) }, ignoreNULL = FALSE)

  output$cluster_tools_panel <- renderUI({
    ca <- cluster_assign_rv()
    req(ca)
    ks <- sort(unique(as.integer(ca)))
    tagList(
      hr(),
      h4("\U0001f9e9 Clusters"),
      tags$p(style = "color:#1F3864; font-size:0.9em;",
             icon("circle-info"),
             paste0(" ", length(ks), " clusters over ", length(ca),
                    " genes. Send one to the enrichment panel below, ",
                    "or export the membership.")),
      div(
        lapply(ks, function(i) {
          n_i <- sum(as.integer(ca) == i)
          actionButton(paste0("clust_sel_", i),
                       paste0("Cluster ", i, " (", n_i, ")"),
                       class = "btn-sm btn-default",
                       style = "margin: 2px;")
        })
      ),
      br(),
      fluidRow(
        column(4, downloadButton("dl_cluster_csv", "\u2b07 Membership (CSV)",
                                 class = "btn-sm btn-info", width = "100%")),
        column(4, actionButton("run_stability", "\U0001f9ea Cluster stability",
                               class = "btn-sm btn-warning", width = "100%")),
        column(4, helpText("Resamples genes and re-clusters; higher Jaccard = more reproducible."))
      ),
      uiOutput("cluster_stability_out")
    )
  })

  # One observer per possible cluster. k_clusters maxes out at 10 (UI slider),
  # so these are registered once at startup rather than re-created on each run —
  # re-registering observers inside renderUI leaks a new one per render.
  lapply(seq_len(10), function(i) {
    observeEvent(input[[paste0("clust_sel_", i)]], {
      ca <- cluster_assign_rv()
      req(ca)
      genes <- names(ca)[as.integer(ca) == i]
      req(length(genes) > 0)
      # Reuses the existing selection channel, so GO/KEGG/GSEA/STRING all pick
      # this up with no further wiring — same path the brush uses.
      selected_genes(genes)
      showNotification(paste0("\u2705 Cluster ", i, ": ", length(genes),
                              " genes sent to enrichment."),
                       type = "message", duration = 4)
    }, ignoreInit = TRUE)
  })

  output$dl_cluster_csv <- downloadHandler(
    filename = function() paste0("cluster_membership_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content  = function(f) {
      ca <- cluster_assign_rv()
      req(ca)
      gl  <- gene_labels_rv()
      ids <- names(ca)
      sym <- ids
      if (!is.null(gl) && length(gl) > 0) {
        idx <- ids_match(ids, names(gl))
        ok  <- !is.na(idx)
        if (any(ok)) sym[ok] <- as.character(gl[idx[ok]])
      }
      utils::write.csv(
        data.frame(gene_id = ids, gene_label = sym, cluster = as.integer(ca),
                   stringsAsFactors = FALSE),
        f, row.names = FALSE)
    }
  )

  observeEvent(input$run_stability, {
    ca  <- cluster_assign_rv()
    mat <- cluster_matrix_rv()
    if (is.null(ca) || is.null(mat)) {
      showNotification("\u274c Generate a heatmap with 'Cut dendrogram' first.",
                       type = "error")
      return()
    }
    # Align the substrate to the partition by name; refuse rather than guess.
    idx <- match(names(ca), rownames(mat))
    if (anyNA(idx)) {
      showNotification("\u274c Cluster membership no longer matches the current matrix.",
                       type = "error")
      return()
    }
    withProgress(message = "Resampling clusters...", value = 0.5, {
      res <- cluster_stability(mat[idx, , drop = FALSE],
                               assign  = ca,
                               k       = length(unique(ca)),
                               linkage = input$linkage %||% "complete",
                               B       = 50L,
                               # Fixed so repeated clicks on the same partition
                               # return the same numbers. Passed explicitly
                               # rather than left to the default so the value
                               # is visible at the call site and lands in the
                               # reproducibility bundle via attr(res, "seed").
                               seed    = 1L)
      cluster_stability_rv(res)
    })
    if (is.null(cluster_stability_rv()))
      showNotification("\u26a0\ufe0f Stability could not be computed for this partition.",
                       type = "warning", duration = 8)
  })

  output$cluster_stability_out <- renderUI({
    res <- cluster_stability_rv()
    if (is.null(res)) return(NULL)
    # Band of a single value. Kept as a local so the verdict shown for the mean
    # and the verdicts checked at the interval edges cannot drift apart.
    band <- function(x) ifelse(is.na(x), "\u2014",
                        ifelse(x > 0.75, "high",
                        ifelse(x < 0.6, "low", "moderate")))
    lo <- res$stability - 2 * res$se
    hi <- res$stability + 2 * res$se
    tagList(
      br(),
      # Wording corrected to describe what is actually computed. The routine
      # draws with replacement and then de-duplicates, which leaves ~63% of the
      # genes (1 - 1/e) -- a subsample, not a bootstrap. Hennig's > 0.75 /
      # < 0.6 rules of thumb are calibrated for resampling WITH replacement, so
      # they are not quoted as verdicts here; the same cut points are kept as
      # relative bands, which is what the numbers can honestly support.
      tags$p(style = "font-size:0.85em; color:#555;",
             "Mean Jaccard \u00b1 2 standard errors over ",
             attr(res, "B") %||% 50, " resamples, each re-clustering roughly 63% of the genes. ",
             "Compare clusters against each other: ", tags$b("higher = more reproducible"),
             ". Treat the bands as indicative, not as a significance test."),
      # The interval is what makes the verdict honest. A cluster whose mean sits
      # near a cut point is not "high" or "moderate" -- the resampling cannot
      # tell, and saying so is more useful than picking whichever side the mean
      # happened to land on. Marked "borderline" whenever +-2 SE spans two bands.
      tags$p(style = "font-size:0.85em; color:#555;",
             "Seeded (seed = ", attr(res, "seed") %||% "none",
             "), so repeating this on the same partition reproduces these numbers. ",
             tags$b("borderline"), " means the interval spans two bands \u2014 the ",
             "resampling cannot separate them at this number of resamples."),
      renderTable({
        data.frame(
          cluster     = res$cluster,
          n_genes     = res$n_genes,
          stability   = ifelse(is.na(res$se),
                               formatC(res$stability, format = "f", digits = 3),
                               sprintf("%.3f \u00b1 %.3f", res$stability, 2 * res$se)),
          resamples   = res$n_resamples,
          reproducibility = ifelse(
            is.na(res$se) | band(lo) == band(hi),
            band(res$stability),
            paste0(band(lo), "\u2013", band(hi), " (borderline)")),
          stringsAsFactors = FALSE, check.names = FALSE
        )
      })
    )
  })

  output$cart_heatmap_panel <- renderUI({
    req(cart_heatmap_genes())
    genes <- cart_heatmap_genes()
    tagList(
      hr(),
      h4("\U0001f525 Custom Gene List Heatmap"),
      tags$p(
        style = "color:#1F3864; font-size:0.9em;",
        icon("circle-info"),
        paste(length(genes), "genes from GO Heatmap cart")
      ),
      fluidRow(
        column(3, sliderInput("cart_ht_h",  "Height (px)", 200, 1400, 500, 50, ticks = FALSE)),
        column(3, sliderInput("cart_ht_fs", "Font size",   6,   24,   12,  1,  ticks = FALSE)),
        column(3, br(), downloadButton("dl_cart_ht_pdf", "\u2b07 PDF",
                                       class = "btn-sm btn-success")),
        column(3, br(), actionButton("close_cart_heatmap", "\u2716 Close",
                                     class = "btn-sm btn-danger", width = "100%"))
      ),
      uiOutput("cart_heatmap_container")
    )
  })

  output$cart_heatmap_container <- renderUI({
    req(cart_heatmap_genes())
    h <- if (is.null(input$cart_ht_h)) 500 else input$cart_ht_h
    plotOutput("cart_heatmap_plot", height = paste0(h, "px"))
  })

  output$cart_heatmap_plot <- renderPlot({
    req(cart_heatmap_genes(), scaled_matrix_rv())
    mat   <- scaled_matrix_rv()          # DISPLAY matrix — color only
    genes <- canon_ids(cart_heatmap_genes())
    keep  <- ids_in(rownames(mat), genes)
    req(sum(keep) > 0)
    sub_mat <- mat[keep, , drop = FALSE]
    
    ord <- col_order_rv()
    if (!is.null(ord)) {
      ord_valid <- intersect(ord, colnames(sub_mat))
      if (length(ord_valid) == ncol(sub_mat))
        sub_mat <- sub_mat[, ord_valid, drop = FALSE]
    }
    
    # Re-cluster this subset on the SAME substrate the main heatmap uses
    # (cluster_matrix_rv, VST/expression-scale) — not on sub_mat, which is
    # display-normalized. Rows/columns are re-aligned to sub_mat by name so
    # the hclust object always matches what Heatmap() below actually draws.
    # Falls back to sub_mat if cluster_matrix_rv is unavailable, so the panel
    # degrades gracefully instead of erroring.
    clust_src <- cluster_matrix_rv()
    sub_clust <- NULL
    if (!is.null(clust_src)) {
      idx <- match(rownames(sub_mat), rownames(clust_src))
      if (!anyNA(idx)) {
        common_cols <- intersect(colnames(clust_src), colnames(sub_mat))
        if (length(common_cols) == ncol(sub_mat)) {
          sub_clust <- clust_src[idx, common_cols, drop = FALSE]
        }
      }
    }
    if (is.null(sub_clust)) sub_clust <- sub_mat  # fallback only
    
    # Same contract as the main heatmap (>= 3 conditions, no constant rows,
    # same linkage). This matters most on the `sub_clust <- sub_mat` fallback
    # above: that path clusters the DISPLAY matrix, where zscore_col can leave
    # rows that are constant on the raw scale and make cor() return NA.
    # $reason is intentionally not surfaced here: renderPlot re-runs on every
    # font-size / height slider tick and would spam notifications.
    cluster_rows_cart <- cluster_rows_or_reason(
      sub_clust, linkage = input$linkage %||% "complete")$hc

    col_fun    <- col_fun_rv()
    gl         <- gene_labels_rv()
    row_labels <- rownames(sub_mat)

    degs_cart   <- deg_genes_rv()
    is_deg_cart <- ids_in(rownames(sub_mat), degs_cart)
    cart_col    <- ifelse(is_deg_cart, "#FF00FF", "black")
    cart_face   <- ifelse(is_deg_cart, "bold",    "plain")

    if (!is.null(gl) && length(gl) > 0) {
      idx   <- ids_match(rownames(sub_mat), names(gl))
      valid <- !is.na(idx)
      if (any(valid)) row_labels[valid] <- as.character(gl[idx[valid]])
    }
    
    fs <- if (is.null(input$cart_ht_fs)) 12 else input$cart_ht_fs

    ht <- Heatmap(
      sub_mat,
      name             = "zscore",
      col              = col_fun,
      cluster_rows     = cluster_rows_cart,
      cluster_columns  = FALSE,
      show_row_names   = TRUE,
      row_labels       = row_labels,
      row_names_gp     = gpar(fontsize = fs, col = cart_col, fontface = cart_face),
      column_names_gp  = gpar(fontsize = fs, fontface = "bold"),
      heatmap_legend_param = list(
        title     = scale_label_rv(),
        title_gp  = gpar(fontsize = fs, fontface = "bold"),
        labels_gp = gpar(fontsize = fs - 1)
      )
    )
    grid::grid.newpage()
    draw(ht)
  })

  output$dl_cart_ht_pdf <- downloadHandler(
    filename = function() paste0("cart_heatmap_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf"),
    content  = function(f) {
      req(cart_heatmap_genes(), scaled_matrix_rv())
      mat   <- scaled_matrix_rv()          # DISPLAY matrix — color only
      genes <- canon_ids(cart_heatmap_genes())
      keep  <- ids_in(rownames(mat), genes)
      sub_mat <- mat[keep, , drop = FALSE]
      
      ord <- col_order_rv()
      if (!is.null(ord)) {
        ord_valid <- intersect(ord, colnames(sub_mat))
        if (length(ord_valid) == ncol(sub_mat))
          sub_mat <- sub_mat[, ord_valid, drop = FALSE]
      }
      
      # Same substrate fix as the live cart heatmap (see the renderPlot
      # handler above) — keeps the downloaded PDF's row order consistent with
      # what is shown on screen.
      clust_src <- cluster_matrix_rv()
      sub_clust <- NULL
      if (!is.null(clust_src)) {
        idx <- match(rownames(sub_mat), rownames(clust_src))
        if (!anyNA(idx)) {
          common_cols <- intersect(colnames(clust_src), colnames(sub_mat))
          if (length(common_cols) == ncol(sub_mat)) {
            sub_clust <- clust_src[idx, common_cols, drop = FALSE]
          }
        }
      }
      if (is.null(sub_clust)) sub_clust <- sub_mat  # fallback only
      
      # Same helper call as the live cart heatmap above, so the downloaded PDF
      # and the on-screen panel can never disagree about row order.
      cluster_rows_cart <- cluster_rows_or_reason(
        sub_clust, linkage = input$linkage %||% "complete")$hc

      col_fun    <- col_fun_rv()
      gl         <- gene_labels_rv()
      row_labels <- rownames(sub_mat)
      if (!is.null(gl) && length(gl) > 0) {
        idx   <- ids_match(rownames(sub_mat), names(gl))
        valid <- !is.na(idx)
        if (any(valid)) row_labels[valid] <- as.character(gl[idx[valid]])
      }
      fs <- if (is.null(input$cart_ht_fs)) 12 else input$cart_ht_fs

      ht <- Heatmap(sub_mat,
                    name            = "zscore",
                    col             = col_fun,
                    cluster_rows    = cluster_rows_cart,
                    cluster_columns = FALSE,
                    show_row_names  = TRUE,
                    row_labels      = row_labels,
                    row_names_gp    = gpar(fontsize = fs),
                    column_names_gp = gpar(fontsize = fs, fontface = "bold"),
                    heatmap_legend_param = list(title = scale_label_rv()))
      pdf(f, width = 8, height = max(4, nrow(sub_mat) * 0.25))
      draw(ht)
      dev.off()
    }
  )

  observeEvent(input$close_cart_heatmap, {
    cart_heatmap_genes(NULL)
  })
}

shinyApp(ui, server)
