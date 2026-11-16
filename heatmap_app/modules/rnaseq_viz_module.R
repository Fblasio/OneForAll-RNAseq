# =============================================================================
# rnaseq_viz_module.R
# -----------------------------------------------------------------------------
# Shiny module for exploratory visualisation of the RNA-seq dataset.
# Currently implements PCA using Variance Stabilising Transformation (VST).
#
# Inputs (reactive arguments from rnaseqDdsServer):
#   dds_rv      — fitted DESeqDataSet
#   metadata_rv — saved metadata table (used to populate the colour-by selector)
#
# The PCA is computed with DESeq2::vst(blind = TRUE) followed by
# DESeq2::plotPCA(..., returnData = TRUE), which returns a data frame
# that is then passed to a custom ggplot2 layer for full visual control.
# Sample labels use ggrepel to avoid overlapping text.
# =============================================================================


# -----------------------------------------------------------------------------
# rnaseqVizControlsUI(): two-row control bar rendered above the plot.
# Split from the main plot UI so it can be placed in a different panel
# if the layout requires it.
# -----------------------------------------------------------------------------
rnaseqVizControlsUI <- function(id) {
ns <- NS(id)

tagList(
fluidRow(
column(
4,
uiOutput(ns("pca_group_ui"))        # dynamic: populated from metadata columns
),
column(
4,
numericInput(ns("pca_ntop"), lbl("Top variable genes", "The plot uses only this many genes, ranked by variance across samples - so the number is part of the result, not a display setting. 500 is standard. If the picture changes a lot between 250 and 2000, treat the structure with caution."), value = 500, min = 100, step = 100)
# number of most variable genes used to compute the PCA
),
column(
4,
numericInput(ns("pca_point_size"), "Point size", value = 3, min = 1, max = 10, step = 0.5)
)
),
fluidRow(
column(
4,
numericInput(ns("pca_text_size"), "Text size", value = 13, min = 6, max = 24, step = 1)
# base font size applied to all ggplot2 text elements
),
column(
4,
checkboxInput(ns("pca_show_labels"), "Show sample labels", value = FALSE)
# when TRUE, adds ggrepel labels with sample names
),
column(
4,
sliderInput(ns("pca_height"), "Plot height (px)", min = 300, max = 1400, value = 700, step = 50)
)
),
textInput(ns("pca_title"), "PCA title", value = "PCA plot")
)
}


# -----------------------------------------------------------------------------
# rnaseqVizPlotUI(): full plot panel combining the controls and the plot output.
# plotOutput height is set dynamically via the pca_height slider.
# -----------------------------------------------------------------------------
rnaseqVizPlotUI <- function(id) {
ns <- NS(id)

tagList(
h3("PCA"),
# Folded away on demand so the figure can sit at the top of the tab: the
# controls are two full rows plus a title field, which otherwise push the plot
# below the fold on a laptop screen. Expanded by default, so nothing moves for
# anyone who was not looking for this. See collapsible_panel() in global.R.
collapsible_panel("\U0001f39b️ PCA options",
                  rnaseqVizControlsUI(id),
                  open = TRUE),
br(),
plotOutput(ns("pca_plot"), height = "auto")
)
}


# -----------------------------------------------------------------------------
# SERVER
# rnaseqVizServer(): computes VST, extracts PCA data, builds ggplot2 figure.
# -----------------------------------------------------------------------------
rnaseqVizServer <- function(id, dds_rv, metadata_rv) {
moduleServer(id, function(input, output, session) {

  # ---------------------------------------------------------------------------
  # output$pca_group_ui: dynamic selectInput populated from the metadata columns.
  # Excludes the "sample" column since that is the row identifier, not a group.
  # ---------------------------------------------------------------------------
  output$pca_group_ui <- renderUI({
    req(metadata_rv())
    cols <- setdiff(colnames(metadata_rv()), "sample")
    if (length(cols) == 0) return(NULL)

    selectInput(
      session$ns("pca_group"),
      "Colour by",
      choices  = cols,
      selected = cols[1],
      width    = "100%"
    )
  })

  # ---------------------------------------------------------------------------
  # pca_data: reactive that computes VST and returns PCA coordinates.
  # blind = TRUE: VST ignores the design formula, recommended for QC plots.
  # returnData = TRUE: returns a data frame instead of a ggplot object,
  #   allowing full custom styling. The data frame contains:
  #     PC1, PC2, name (sample ID), group (selected colour variable)
  #   The "percentVar" attribute holds the fraction of variance per PC.
  # Re-runs whenever dds_rv, pca_group, or pca_ntop change.
  # ---------------------------------------------------------------------------
  pca_data <- reactive({
    req(dds_rv())
    req(input$pca_group)
    dds <- dds_rv()
    req(input$pca_group %in% colnames(SummarizedExperiment::colData(dds)))
    vsd <- tryCatch(
      DESeq2::vst(dds, blind = TRUE, nsub = min(1000, nrow(dds))),
      error = function(e)
        DESeq2::varianceStabilizingTransformation(dds, blind = TRUE)
    )
    
    DESeq2::plotPCA(
      vsd,
      intgroup   = input$pca_group,
      ntop       = input$pca_ntop,
      returnData = TRUE
    )
  })

  # ---------------------------------------------------------------------------
  # output$pca_plot: builds the final ggplot2 PCA figure.
  #
  # Layout decisions:
  #   - 18% padding added to both axes so points never touch the panel border.
  #   - coord_cartesian(clip = "off") prevents labels near edges from being
  #     clipped by the panel boundary.
  #   - plot.margin (1.2 cm top/bottom, 2.5 cm left/right) provides room
  #     for long sample names that extend beyond the plot area.
  #   - theme_minimal with pca_text_size applied uniformly to all text elements
  #     (axis text, axis titles, legend title, legend text).
  #
  # Optional ggrepel labels (when pca_show_labels is TRUE):
  #   - Uses the "name" column returned by plotPCA() as the label.
  #   - size is derived from pca_text_size / 3.2 to match ggplot pt units.
  #   - box.padding and point.padding control minimum spacing between labels
  #     and between labels and points.
  #   - max.overlaps = Inf ensures no label is ever hidden.
  #   - min.segment.length = 0 always draws the connecting segment.
  # ---------------------------------------------------------------------------
  output$pca_plot <- renderPlot({
    req(pca_data())

    pca_df      <- pca_data()
    percent_var <- round(100 * attr(pca_df, "percentVar"))

    # Compute axis ranges with 18% padding to avoid points touching the border
    x_rng <- range(pca_df$PC1, na.rm = TRUE)
    y_rng <- range(pca_df$PC2, na.rm = TRUE)
    x_pad <- max(diff(x_rng) * 0.18, 1)
    y_pad <- max(diff(y_rng) * 0.18, 1)

    p <- ggplot2::ggplot(
      pca_df,
      ggplot2::aes(
        x     = PC1,
        y     = PC2,
        color = .data[[input$pca_group]]
      )
    ) +
      ggplot2::geom_point(size = input$pca_point_size) +
      ggplot2::labs(
        title = input$pca_title,
        x     = paste0("PC1: ", percent_var[1], "% variance"),
        y     = paste0("PC2: ", percent_var[2], "% variance"),
        color = input$pca_group
      ) +
      # clip = "off" prevents ggrepel labels near the panel border from being cut
      ggplot2::coord_cartesian(
        xlim = c(x_rng[1] - x_pad, x_rng[2] + x_pad),
        ylim = c(y_rng[1] - y_pad, y_rng[2] + y_pad),
        clip = "off"
      ) +
      ggplot2::theme_minimal(base_size = input$pca_text_size) +
      ggplot2::theme(
        plot.title    = ggplot2::element_text(face = "bold", hjust = 0.5),
        axis.title    = ggplot2::element_text(size = input$pca_text_size + 1),
        axis.text     = ggplot2::element_text(size = input$pca_text_size),
        legend.title  = ggplot2::element_text(size = input$pca_text_size),
        legend.text   = ggplot2::element_text(size = input$pca_text_size),
        # generous margins prevent long sample labels from being clipped
        plot.margin   = grid::unit(c(1.2, 2.5, 1.2, 2.5), "cm")
      )

    # --- Optional: ggrepel sample labels ------------------------------------
    if (isTRUE(input$pca_show_labels)) {
      p <- p + ggrepel::geom_text_repel(
        ggplot2::aes(label = name),
        size               = input$pca_text_size / 3.2,  # convert pt to ggplot size units
        box.padding        = 0.6,     # padding around label bounding box
        point.padding      = 0.35,    # minimum gap between label and point
        max.overlaps       = Inf,     # never hide a label, even in dense areas
        min.segment.length = 0,       # always draw connecting segment
        show.legend        = FALSE
      )
    }

    p

  # Plot height is driven by the pca_height slider (in pixels)
  }, height = function() {
    input$pca_height
  })
})
}
