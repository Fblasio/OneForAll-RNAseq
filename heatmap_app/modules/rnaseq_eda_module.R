# =============================================================================
# rnaseq_eda_module.R  (goes in modules/, sourced from global.R after
#                       rnaseq_viz_module.R)
# -----------------------------------------------------------------------------
# Complements rnaseq_viz_module.R (PCA, already in the app) with the other
# QC/EDA views reviewers expect as table stakes for an RNA-seq tool:
#   - MDS            (classical multidimensional scaling on sample distances)
#   - Sample corr.   (pairwise correlation heatmap, ComplexHeatmap-based —
#                      same visual language as the app's main heatmap)
#   - Batch diag.    (variance of each PC explained by each metadata factor —
#                      a lightweight PVCA-style diagnostic: lets a user see
#                      "does PC1 track my batch column, or my condition
#                      column?" without a new Bioconductor dependency)
#
# All three reuse the SAME "VST -> top-variable-genes" preprocessing that
# rnaseq_viz_module.R already uses for PCA, so the three views are directly
# comparable (same gene subset feeds PCA, MDS, correlation and batch calcs).
#
# No new package dependencies: DESeq2/SummarizedExperiment (already used for
# the dds), stats (dist/cmdscale/prcomp/lm — all base), ggplot2,
# ComplexHeatmap/circlize (already used by the main heatmap), color_utils.R's
# get_cat_colors() for a consistent annotation palette.
#
# Wiring (see app.r instructions):
#   UI:     rnaseqEdaMdsUI("rnaseq_eda")   / rnaseqEdaCorUI("rnaseq_eda")
#           / rnaseqEdaBatchUI("rnaseq_eda")   — one tabPanel each, alongside
#           the existing PCA tabPanel(rnaseqVizPlotUI("rnaseq_viz")).
#   Server: rnaseqEdaServer("rnaseq_eda", dds_rv = ..., metadata_rv = ...)
# =============================================================================


# -----------------------------------------------------------------------------
# Internal helpers (module-private computation; no Shiny reactivity here so
# they're easy to reason about/unit-test in isolation).
# -----------------------------------------------------------------------------

# VST-transform + subset to the ntop most-variable genes. Mirrors exactly what
# rnaseq_viz_module.R's PCA already does, so MDS/correlation/batch views are
# computed on the SAME gene subset as the existing PCA plot.
.eda_vst_matrix <- function(dds, ntop = 500) {
  vsd <- tryCatch(
    DESeq2::vst(dds, blind = TRUE, nsub = min(1000, nrow(dds))),
    error = function(e) DESeq2::varianceStabilizingTransformation(dds, blind = TRUE)
  )
  mat <- SummarizedExperiment::assay(vsd)
  if (!is.null(ntop) && is.finite(ntop) && ntop > 0 && ntop < nrow(mat)) {
    rv  <- apply(mat, 1, stats::var)
    sel <- order(rv, decreasing = TRUE)[seq_len(ntop)]
    mat <- mat[sel, , drop = FALSE]
  }
  mat
}

# Classical MDS on Euclidean sample-to-sample distances (base stats only —
# no limma::plotMDS dependency). Returns a data.frame with Dim1/Dim2/name,
# and a "percentVar" attribute analogous to DESeq2::plotPCA's, so the MDS
# plot can reuse the exact same rendering logic as the existing PCA plot.
.eda_compute_mds <- function(mat, k = 2) {
  if (is.null(mat) || ncol(mat) < 3) return(NULL)   # need >= 3 samples for MDS
  d   <- stats::dist(t(mat))
  fit <- tryCatch(stats::cmdscale(d, k = min(k, ncol(mat) - 1), eig = TRUE),
                  error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  pts <- as.data.frame(fit$points)
  colnames(pts) <- paste0("Dim", seq_len(ncol(pts)))
  pts$name <- colnames(mat)
  eig <- fit$eig
  percent_var <- if (!is.null(eig)) {
    pos <- pmax(eig, 0)
    if (sum(pos) > 0) round(100 * pos / sum(pos), 1) else rep(NA_real_, length(eig))
  } else {
    rep(NA_real_, ncol(pts) - 1)
  }
  attr(pts, "percentVar") <- percent_var
  pts
}

# Sample-by-sample correlation matrix.
.eda_compute_cor <- function(mat, method = "pearson") {
  stats::cor(mat, method = method)
}

# PCA scores for the first n_pcs components (base prcomp — unlike
# DESeq2::plotPCA, which only ever returns PC1/PC2, this exposes as many PCs
# as the batch-diagnostic view needs).
.eda_compute_pca_scores <- function(mat, n_pcs = 5) {
  if (is.null(mat) || ncol(mat) < 3) return(NULL)
  pca <- stats::prcomp(t(mat), center = TRUE, scale. = FALSE)
  k   <- min(n_pcs, ncol(pca$x))
  scores <- as.data.frame(pca$x[, seq_len(k), drop = FALSE])
  var_explained <- (pca$sdev^2 / sum(pca$sdev^2))[seq_len(k)]
  list(scores = scores, percentVar = round(100 * var_explained, 1))
}

# The batch/covariate diagnostic itself: for every (PC, metadata factor)
# pair, fit PC_score ~ factor and take the R^2. A factor that tracks a PC
# strongly (R^2 near 1) is a candidate confound (if it's a batch/technical
# column) or the expected biological driver (if it's your condition column)
# — the point of the view is to make that distinction visible, not to make
# it for the user.
.eda_pc_factor_r2 <- function(pca_scores, metadata, factor_cols) {
  pcs <- colnames(pca_scores)
  out <- matrix(NA_real_, nrow = length(pcs), ncol = length(factor_cols),
                dimnames = list(pcs, factor_cols))
  for (fc in factor_cols) {
    fvals <- metadata[[fc]]
    if (is.null(fvals) || length(unique(stats::na.omit(fvals))) < 2) next
    for (pc in pcs) {
      fit <- tryCatch(stats::lm(pca_scores[[pc]] ~ as.factor(fvals)),
                      error = function(e) NULL)
      if (!is.null(fit)) {
        out[pc, fc] <- tryCatch(summary(fit)$r.squared, error = function(e) NA_real_)
      }
    }
  }
  out
}

# Keep only metadata columns that have a value for EVERY sample.
# Why this matters for the batch diagnostic: stats::lm() defaults to
# na.action = na.omit, so a column with missing values would silently be fitted
# on a SUBSET of samples and its R^2 shown alongside columns fitted on all
# samples — two numbers that aren't comparable. Rather than present that
# side by side, such columns are excluded and reported to the user.
# Treats NA and empty/whitespace-only strings (common when metadata comes from
# a CSV with blank cells) as missing.
.eda_complete_cols <- function(metadata, cols) {
  if (is.null(metadata) || length(cols) == 0) return(character(0))
  keep <- vapply(cols, function(cc) {
    v <- metadata[[cc]]
    if (is.null(v)) return(FALSE)
    ch <- trimws(as.character(v))
    !any(is.na(v) | is.na(ch) | !nzchar(ch))
  }, logical(1))
  cols[keep]
}

# ── Batch-corrected PCA helpers ─────────────────────────────────────────────
# Detect whether `batch_col` and `preserve_col` are confounded (aliased). If
# every level of one occurs with exactly one level of the other, no method can
# separate their effects — removing the batch would also remove the biology.
# The rank of the joint model matrix is the standard test: rank-deficient means
# the two are not independently estimable.
.eda_is_confounded <- function(meta, batch_col, preserve_col) {
  if (is.null(preserve_col) || identical(preserve_col, "__none__")) return(FALSE)
  mm <- tryCatch(
    stats::model.matrix(~ as.factor(meta[[batch_col]]) + as.factor(meta[[preserve_col]])),
    error = function(e) NULL)
  if (is.null(mm)) return(TRUE)
  qr(mm)$rank < ncol(mm)
}

# Remove the effect of `batch_col` from `mat` while PROTECTING the variation
# due to `preserve_col`.
#
# The `design` argument is the critical part: it tells removeBatchEffect which
# variation to keep. Without it, biological signal correlated with batch gets
# stripped out along with the batch effect. When no factor is preserved we pass
# an intercept-only design, which is limma's default and removes more than most
# users expect — hence the UI nudges toward always choosing a factor to preserve.
.eda_remove_batch <- function(mat, meta, batch_col, preserve_col) {
  batch <- as.factor(meta[[batch_col]])
  design <- if (!is.null(preserve_col) && !identical(preserve_col, "__none__")) {
    stats::model.matrix(~ as.factor(meta[[preserve_col]]))
  } else {
    matrix(1, ncol(mat), 1)
  }
  limma::removeBatchEffect(mat, batch = batch, design = design)
}

# PCA scores -> tidy data.frame for plotting, with percentVar attached.
.eda_pca_plot_df <- function(mat, meta, color_col, panel_label) {
  pca <- .eda_compute_pca_scores(mat, n_pcs = 2)
  if (is.null(pca)) return(NULL)
  pct <- pca$percentVar
  df <- data.frame(
    PC1   = pca$scores[[1]],
    PC2   = pca$scores[[2]],
    name  = colnames(mat),
    group = meta[[color_col]],
    stringsAsFactors = FALSE
  )
  df$panel <- sprintf("%s  (PC1 %.1f%%, PC2 %.1f%%)", panel_label, pct[1], pct[2])
  df
}

# Align a metadata data.frame's rows to the column order of an expression
# matrix by sample name. Every EDA view needs this so a factor's levels line
# up with the correct sample, regardless of the order metadata was saved in.
.eda_align_metadata <- function(mat, metadata) {
  if (is.null(metadata) || !"sample" %in% colnames(metadata)) return(NULL)
  idx <- match(colnames(mat), as.character(metadata$sample))
  if (any(is.na(idx))) return(NULL)
  metadata[idx, , drop = FALSE]
}


# =============================================================================
# MDS
# =============================================================================
rnaseqEdaMdsUI <- function(id) {
  ns <- NS(id)
  tagList(
    h4("Multidimensional scaling (MDS)"),
    helpText("Classical MDS on Euclidean distances between VST-transformed samples — ",
             "a distance-based complement to PCA. Points close together are transcriptionally similar."),
    fluidRow(
      column(4, uiOutput(ns("mds_group_ui"))),
      column(4, numericInput(ns("mds_ntop"), lbl("Top variable genes", "The plot uses only this many genes, ranked by variance across samples - so the number is part of the result, not a display setting. 500 is standard. If the picture changes a lot between 250 and 2000, treat the structure with caution."), value = 500, min = 100, step = 100)),
      column(4, checkboxInput(ns("mds_show_labels"), "Show sample labels", value = FALSE))
    ),
    # ── Personalization — same always-visible layout as the PCA tab
    # (rnaseqVizControlsUI in rnaseq_viz_module.R), not the collapsible
    # panel used in GSEA: this view sits right next to PCA, so its
    # controls should look and behave the same way PCA's already do.
    fluidRow(
      column(4, numericInput(ns("mds_point_size"), "Point size", value = 3, min = 1, max = 10, step = 0.5)),
      column(4, numericInput(ns("mds_text_size"), "Text size", value = 13, min = 6, max = 24, step = 1)),
      column(4, sliderInput(ns("mds_height"), "Plot height (px)", min = 300, max = 1400, value = 600, step = 50))
    ),
    textInput(ns("mds_title"), "MDS title", value = "MDS plot"),
    downloadButton(ns("dl_mds"), "\u2b07 Download plot (PNG)", class = "btn-sm btn-success"),
    br(), br(),
    plotOutput(ns("mds_plot"), height = "auto")
  )
}

# =============================================================================
# Sample correlation heatmap
# =============================================================================
rnaseqEdaCorUI <- function(id) {
  ns <- NS(id)
  tagList(
    h4("Sample correlation"),
    helpText("Pairwise correlation between samples on VST-transformed, top-variable genes. ",
             "Replicates should cluster tightly with high correlation (typically > 0.9)."),
    fluidRow(
      column(3, selectInput(ns("cor_method"), "Method",
                            choices = c("Pearson" = "pearson", "Spearman" = "spearman"),
                            selected = "pearson")),
      column(3, numericInput(ns("cor_ntop"), lbl("Top variable genes", "The plot uses only this many genes, ranked by variance across samples - so the number is part of the result, not a display setting. 500 is standard. If the picture changes a lot between 250 and 2000, treat the structure with caution."), value = 500, min = 100, step = 100)),
      column(3, uiOutput(ns("cor_annot_ui"))),
      column(3, checkboxInput(ns("cor_show_values"), "Show r values", value = TRUE))
    ),
    # ── Personalization — same always-visible layout as PCA's controls ──
    fluidRow(
      column(4, colourpicker::colourInput(ns("cor_col_low"),  "Low",  value = "#3B7DD8")),
      column(4, colourpicker::colourInput(ns("cor_col_mid"),  "Mid",  value = "white")),
      column(4, colourpicker::colourInput(ns("cor_col_high"), "High", value = "#C0392B"))
    ),
    fluidRow(
      column(4, numericInput(ns("cor_font_size"), "Label font size", value = 10, min = 6, max = 20, step = 1)),
      column(4, checkboxInput(ns("cor_cluster_rows"), "Cluster rows", value = TRUE)),
      column(4, checkboxInput(ns("cor_cluster_cols"), "Cluster columns", value = TRUE))
    ),
    textInput(ns("cor_title"), "Heatmap title", value = ""),
    downloadButton(ns("dl_cor"), "\u2b07 Download plot (PNG)", class = "btn-sm btn-success"),
    br(), br(),
    plotOutput(ns("cor_plot"), height = "600px")
  )
}

# =============================================================================
# Batch / covariate diagnostic
# =============================================================================
rnaseqEdaBatchUI <- function(id) {
  ns <- NS(id)
  tagList(
    h4("Batch / covariate diagnostic"),
    helpText("For each principal component and each metadata column, the tile shows R\u00b2 ",
             "from a one-way model of that PC on that factor \u2014 how much of that PC's ",
             "variance the factor explains. A technical column (batch, lane, date) lighting up ",
             "on an early PC is a signal to consider including it in your DESeq2 design formula. ",
             "This is a diagnostic, not a p-value \u2014 always look at sample size per group too."),
    fluidRow(
      column(4, numericInput(ns("batch_ntop"), "Top variable genes", value = 500, min = 100, step = 100)),
      column(4, numericInput(ns("batch_npcs"), "Number of PCs", value = 5, min = 2, max = 15, step = 1))
    ),
    # ── Personalization — same always-visible layout as PCA's controls ──
    fluidRow(
      column(4, colourpicker::colourInput(ns("batch_col_low"),  "Low R\u00b2",  value = "white")),
      column(4, colourpicker::colourInput(ns("batch_col_high"), "High R\u00b2", value = "#C0392B")),
      column(4, numericInput(ns("batch_text_size"), "Text size", value = 13, min = 6, max = 24, step = 1))
    ),
    textInput(ns("batch_title"), "Plot title",
             value = "Variance of each PC explained by each metadata factor"),
    downloadButton(ns("dl_batch"), "\u2b07 Download plot (PNG)", class = "btn-sm btn-success"),
    br(), br(),
    uiOutput(ns("batch_dropped_note")),
    plotOutput(ns("batch_plot"), height = "500px")
  )
}


# =============================================================================
# PCA with batch correction (VISUALIZATION ONLY)
# =============================================================================
rnaseqEdaPcaCorrectedUI <- function(id) {
  ns <- NS(id)
  tagList(
    h4("PCA with batch correction"),
    div(style = "color:#8a6d3b; background:#fcf8e3; border:1px solid #faebcc;
                 border-radius:4px; padding:10px; margin-bottom:12px;",
        strong("Visualization only \u2014 do not use this to run differential expression."),
        br(),
        "This removes a factor's effect from the expression matrix so you can SEE what ",
        "the data looks like without it. To actually account for a batch effect in your ",
        "statistics, add the factor to your DESeq2 design formula instead ",
        "(e.g. ", code("~ batch + condition"), ") and leave the counts untouched. ",
        "Pre-correcting counts and then testing them produces anti-conservative p-values.",
        br(), br(),
        em("Also note: a corrected PCA always looks cleaner \u2014 that is what the correction ",
           "does. It is not evidence that correcting was appropriate.")),
    # The controls fold; the warning box above does NOT. Hiding a caveat behind a
    # click is exactly the wrong thing to make collapsible.
    collapsible_panel(
      "\U0001f39b️ PCA options",
      fluidRow(
        column(4, uiOutput(ns("bc_batch_ui"))),
        column(4, uiOutput(ns("bc_preserve_ui"))),
        column(4, uiOutput(ns("bc_color_ui")))
      ),
      fluidRow(
        column(3, numericInput(ns("bc_ntop"), lbl("Top variable genes", "The plot uses only this many genes, ranked by variance across samples - so the number is part of the result, not a display setting. 500 is standard. If the picture changes a lot between 250 and 2000, treat the structure with caution."), value = 500, min = 100, step = 100)),
        column(3, numericInput(ns("bc_point_size"), "Point size", value = 3, min = 1, max = 10, step = 0.5)),
        column(3, numericInput(ns("bc_text_size"), "Text size", value = 13, min = 6, max = 24, step = 1)),
        column(3, br(), checkboxInput(ns("bc_show_labels"), "Show sample labels", value = FALSE))
      ),
      open = TRUE
    ),
    downloadButton(ns("dl_bc"), "\u2b07 Download plot (PNG)", class = "btn-sm btn-success"),
    br(), br(),
    uiOutput(ns("bc_status")),
    plotOutput(ns("bc_plot"), height = "550px")
  )
}


# =============================================================================
# SERVER
# =============================================================================
rnaseqEdaServer <- function(id, dds_rv, metadata_rv) {
  moduleServer(id, function(input, output, session) {

    metadata_cols <- reactive({
      req(metadata_rv())
      setdiff(colnames(metadata_rv()), "sample")
    })

    # ---- MDS ------------------------------------------------------------
    output$mds_group_ui <- renderUI({
      cols <- metadata_cols()
      if (length(cols) == 0) return(NULL)
      selectInput(session$ns("mds_group"), "Colour by",
                  choices = cols, selected = cols[1], width = "100%")
    })

    mds_data <- reactive({
      req(dds_rv())
      mat <- .eda_vst_matrix(dds_rv(), ntop = input$mds_ntop %||% 500)
      pts <- .eda_compute_mds(mat)
      req(!is.null(pts))
      meta <- .eda_align_metadata(mat, metadata_rv())
      req(!is.null(meta), input$mds_group %in% colnames(meta))
      pts$group <- meta[[input$mds_group]]
      pts
    })

    build_mds_plot <- function() {
      df <- mds_data()
      pv <- attr(df, "percentVar")
      x_rng <- range(df$Dim1, na.rm = TRUE); y_rng <- range(df$Dim2, na.rm = TRUE)
      x_pad <- max(diff(x_rng) * 0.18, 1); y_pad <- max(diff(y_rng) * 0.18, 1)
      text_size <- input$mds_text_size %||% 13

      p <- ggplot2::ggplot(df, ggplot2::aes(x = Dim1, y = Dim2, color = group)) +
        ggplot2::geom_point(size = input$mds_point_size %||% 3) +
        ggplot2::labs(
          title = input$mds_title %||% "MDS plot",
          x = paste0("Dim1", if (!is.na(pv[1])) paste0(" (", pv[1], "%)") else ""),
          y = paste0("Dim2", if (length(pv) >= 2 && !is.na(pv[2])) paste0(" (", pv[2], "%)") else ""),
          color = input$mds_group
        ) +
        ggplot2::coord_cartesian(
          xlim = c(x_rng[1] - x_pad, x_rng[2] + x_pad),
          ylim = c(y_rng[1] - y_pad, y_rng[2] + y_pad), clip = "off"
        ) +
        ggplot2::theme_minimal(base_size = text_size) +
        ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
                       axis.title = ggplot2::element_text(size = text_size + 1),
                       axis.text  = ggplot2::element_text(size = text_size),
                       legend.title = ggplot2::element_text(size = text_size),
                       legend.text  = ggplot2::element_text(size = text_size),
                       plot.margin = grid::unit(c(1.2, 2.5, 1.2, 2.5), "cm"))

      if (isTRUE(input$mds_show_labels)) {
        p <- p + ggrepel::geom_text_repel(ggplot2::aes(label = name), size = text_size / 3.2,
                                          max.overlaps = Inf, show.legend = FALSE)
      }
      p
    }

    output$mds_plot <- renderPlot({
      req(mds_data()); build_mds_plot()
    }, height = function() { input$mds_height %||% 600 })

    output$dl_mds <- downloadHandler(
      filename = function() paste0("MDS_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
      content = function(file) {
        req(mds_data())
        ggplot2::ggsave(file, plot = build_mds_plot(), width = 9, height = 7, dpi = 300)
      }
    )

    # ---- Sample correlation ----------------------------------------------
    output$cor_annot_ui <- renderUI({
      cols <- metadata_cols()
      selectInput(session$ns("cor_annot"), "Annotate by",
                  choices = c("None" = "__none__", cols),
                  selected = "__none__", width = "100%")
    })

    cor_mat <- reactive({
      req(dds_rv())
      mat <- .eda_vst_matrix(dds_rv(), ntop = input$cor_ntop %||% 500)
      .eda_compute_cor(mat, method = input$cor_method %||% "pearson")
    })

    build_cor_heatmap <- function() {
      m <- cor_mat()
      rng <- range(m, na.rm = TRUE)
      col_fun <- circlize::colorRamp2(
        c(rng[1], mean(rng), rng[2]),
        c(input$cor_col_low %||% "#3B7DD8",
          input$cor_col_mid %||% "white",
          input$cor_col_high %||% "#C0392B"))
      fs <- input$cor_font_size %||% 10

      top_annot <- NULL
      annot_choice <- input$cor_annot %||% "__none__"
      if (!identical(annot_choice, "__none__")) {
        meta <- .eda_align_metadata(m, metadata_rv())
        if (!is.null(meta) && annot_choice %in% colnames(meta)) {
          lvls <- unique(as.character(meta[[annot_choice]]))
          pal  <- get_cat_colors(lvls)
          top_annot <- ComplexHeatmap::HeatmapAnnotation(
            df  = data.frame(setNames(list(as.character(meta[[annot_choice]])), annot_choice)),
            col = setNames(list(pal), annot_choice)
          )
        }
      }

      ComplexHeatmap::Heatmap(
        m, name = "r",
        column_title = if (nzchar(input$cor_title %||% "")) input$cor_title else NULL,
        col = col_fun,
        top_annotation = top_annot,
        cluster_rows    = isTRUE(input$cor_cluster_rows %||% TRUE),
        cluster_columns = isTRUE(input$cor_cluster_cols %||% TRUE),
        cell_fun = if (isTRUE(input$cor_show_values)) {
          function(j, i, x, y, width, height, fill)
            grid::grid.text(sprintf("%.2f", m[i, j]), x, y,
                            gp = grid::gpar(fontsize = fs - 1,
                                            col = if (abs(m[i, j]) > 0.6) "white" else "black"))
        } else NULL,
        row_names_gp = grid::gpar(fontsize = fs),
        column_names_gp = grid::gpar(fontsize = fs)
      )
    }

    output$cor_plot <- renderPlot({ req(cor_mat()); ComplexHeatmap::draw(build_cor_heatmap()) })

    output$dl_cor <- downloadHandler(
      filename = function() paste0("sample_correlation_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
      content = function(file) {
        req(cor_mat())
        grDevices::png(file, width = 9, height = 7, units = "in", res = 300)
        ComplexHeatmap::draw(build_cor_heatmap())
        grDevices::dev.off()
      }
    )

    # ---- Batch diagnostic ---------------------------------------------
    batch_r2 <- reactive({
      req(dds_rv(), metadata_rv())
      mat <- .eda_vst_matrix(dds_rv(), ntop = input$batch_ntop %||% 500)
      pca <- .eda_compute_pca_scores(mat, n_pcs = input$batch_npcs %||% 5)
      req(!is.null(pca))
      meta <- .eda_align_metadata(mat, metadata_rv())
      req(!is.null(meta))
      cols <- metadata_cols()
      req(length(cols) > 0)
      # Only columns with a value for every sample: an incomplete column would
      # be fitted on a subset by lm()'s na.omit and its R^2 wouldn't be
      # comparable to the others (see .eda_complete_cols).
      cols_ok <- .eda_complete_cols(meta, cols)
      dropped <- setdiff(cols, cols_ok)
      shiny::validate(shiny::need(
        length(cols_ok) > 0,
        paste0("No metadata column has a value for every sample, so no batch ",
               "diagnostic can be computed. Incomplete column(s): ",
               paste(dropped, collapse = ", "), ".")))
      list(r2 = .eda_pc_factor_r2(pca$scores, meta, cols_ok),
           pct = pca$percentVar, dropped = dropped)
    })

    build_batch_plot <- function() {
      res <- batch_r2()
      r2  <- res$r2
      pct <- res$pct
      df <- as.data.frame(as.table(r2))
      colnames(df) <- c("PC", "Factor", "R2")
      # Order PCs by variance explained (already in prcomp order) — keep as factor levels.
      df$PC <- factor(df$PC, levels = rev(rownames(r2)))
      pc_labels <- setNames(paste0(rownames(r2), " (", pct, "%)"), rownames(r2))
      fs <- input$batch_text_size %||% 13

      ggplot2::ggplot(df, ggplot2::aes(x = Factor, y = PC, fill = R2)) +
        ggplot2::geom_tile(color = "white") +
        ggplot2::geom_text(ggplot2::aes(label = ifelse(is.na(R2), "", sprintf("%.2f", R2))),
                           size = fs / 3.2) +
        ggplot2::scale_fill_gradient(low  = input$batch_col_low  %||% "white",
                                     high = input$batch_col_high %||% "#C0392B",
                                     limits = c(0, 1), na.value = "grey90",
                                     name = expression(R^2)) +
        ggplot2::scale_y_discrete(labels = pc_labels[levels(df$PC)]) +
        ggplot2::labs(title = input$batch_title %||% "", x = NULL, y = NULL) +
        ggplot2::theme_minimal(base_size = fs) +
        ggplot2::theme(
          plot.title  = ggplot2::element_text(face = "bold", hjust = 0.5, size = fs + 2),
          axis.text.x = ggplot2::element_text(angle = 30, hjust = 1, size = fs),
          axis.text.y = ggplot2::element_text(size = fs)
        )
    }

    output$batch_dropped_note <- renderUI({
      res <- tryCatch(batch_r2(), error = function(e) NULL)
      if (is.null(res) || length(res$dropped) == 0) return(NULL)
      div(style = "color:#8a6d3b; background:#fcf8e3; border:1px solid #faebcc;
                   border-radius:4px; padding:8px; margin-bottom:10px;",
          sprintf("Excluded %d metadata column(s) with missing values: %s. Only columns with a value for every sample are shown, so all R\u00b2 values are computed on the same samples and remain comparable.",
                  length(res$dropped), paste(res$dropped, collapse = ", ")))
    })

    output$batch_plot <- renderPlot({ req(batch_r2()); build_batch_plot() })

    output$dl_batch <- downloadHandler(
      filename = function() paste0("batch_diagnostic_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
      content = function(file) {
        req(batch_r2())
        ggplot2::ggsave(file, plot = build_batch_plot(), width = 8, height = 6, dpi = 300)
      }
    )
    # ---- PCA with batch correction ---------------------------------------
    # Only complete metadata columns are offered: an incomplete batch column
    # would leave NA levels that removeBatchEffect can't model (same reasoning
    # as the batch diagnostic — see .eda_complete_cols).
    bc_cols <- reactive({
      meta <- metadata_rv()
      req(!is.null(meta))
      .eda_complete_cols(meta, metadata_cols())
    })

    output$bc_batch_ui <- renderUI({
      cols <- bc_cols()
      if (length(cols) == 0) return(helpText("No complete metadata column available."))
      selectInput(session$ns("bc_batch"),
                  lbl("Remove effect of (batch)", "Subtracts a technical effect (batch, run, harvest day) FOR THIS PLOT ONLY. Never feed corrected values to differential expression - put the batch in the design formula instead. Use this tab to check whether a batch effect exists."),
                  choices = cols, selected = cols[1], width = "100%")
    })
    output$bc_preserve_ui <- renderUI({
      cols <- bc_cols()
      selectInput(session$ns("bc_preserve"),
                  lbl("Preserve (biology)", "The biological variable to protect while the batch is removed - normally treatment or genotype. Leaving it at 'Nothing' can quietly delete the very signal you are looking for. Must be a different column from the one above."),
                  choices = c("Nothing (not recommended)" = "__none__", cols),
                  selected = if (length(cols) >= 2) cols[2] else "__none__",
                  width = "100%")
    })
    output$bc_color_ui <- renderUI({
      cols <- bc_cols()
      if (length(cols) == 0) return(NULL)
      selectInput(session$ns("bc_color"), "Colour by",
                  choices = cols,
                  selected = if (length(cols) >= 2) cols[2] else cols[1],
                  width = "100%")
    })

    bc_data <- reactive({
      req(dds_rv(), metadata_rv(), input$bc_batch, input$bc_color)
      shiny::validate(shiny::need(
        requireNamespace("limma", quietly = TRUE),
        "The 'limma' package is required for batch correction. Install it with: BiocManager::install('limma')"))

      mat  <- .eda_vst_matrix(dds_rv(), ntop = input$bc_ntop %||% 500)
      meta <- .eda_align_metadata(mat, metadata_rv())
      req(!is.null(meta))

      batch_col    <- input$bc_batch
      preserve_col <- input$bc_preserve %||% "__none__"
      color_col    <- input$bc_color

      shiny::validate(shiny::need(
        !identical(batch_col, preserve_col),
        "'Remove effect of' and 'Preserve' must be different columns \u2014 removing and preserving the same factor is contradictory."))
      shiny::validate(shiny::need(
        length(unique(as.character(meta[[batch_col]]))) >= 2,
        sprintf("Column '%s' has only one level, so there is no batch effect to remove.", batch_col)))
      shiny::validate(shiny::need(
        !.eda_is_confounded(meta, batch_col, preserve_col),
        sprintf(paste0("'%s' and '%s' are confounded (each level of one occurs with only one level of the other). ",
                       "Removing '%s' would also remove the biological signal, so the correction is blocked. ",
                       "This is a limitation of the experimental design, not of the software."),
                batch_col, preserve_col, batch_col)))

      corrected <- tryCatch(.eda_remove_batch(mat, meta, batch_col, preserve_col),
                            error = function(e) e)
      shiny::validate(shiny::need(
        !inherits(corrected, "error"),
        paste("Batch correction failed:", if (inherits(corrected, "error")) conditionMessage(corrected) else "")))

      # Same gene set in both panels so the before/after comparison is fair.
      df_before <- .eda_pca_plot_df(mat,       meta, color_col, "Uncorrected")
      df_after  <- .eda_pca_plot_df(corrected, meta, color_col, "Batch-corrected")
      req(!is.null(df_before), !is.null(df_after))

      out <- rbind(df_before, df_after)
      out$panel <- factor(out$panel, levels = c(unique(df_before$panel), unique(df_after$panel)))
      attr(out, "meta") <- list(batch = batch_col, preserve = preserve_col, color = color_col)
      out
    })

    output$bc_status <- renderUI({
      df <- tryCatch(bc_data(), error = function(e) NULL)
      if (is.null(df)) return(NULL)
      m <- attr(df, "meta")
      if (identical(m$preserve, "__none__")) {
        div(style = "color:#a94442; background:#f2dede; border:1px solid #ebccd1;
                     border-radius:4px; padding:8px; margin-bottom:10px;",
            sprintf("No factor is being preserved, so removing '%s' may also strip out biological signal. Choose a factor under 'Preserve (biology)' unless you have a specific reason not to.", m$batch))
      } else {
        helpText(sprintf("Removed the effect of '%s' while preserving '%s'.", m$batch, m$preserve))
      }
    })

    build_bc_plot <- function() {
      df <- bc_data()
      ts <- input$bc_text_size %||% 13
      m  <- attr(df, "meta")

      p <- ggplot2::ggplot(df, ggplot2::aes(x = PC1, y = PC2, color = group)) +
        ggplot2::geom_point(size = input$bc_point_size %||% 3) +
        ggplot2::facet_wrap(~ panel, scales = "free") +
        ggplot2::labs(title = "PCA before vs after batch correction",
                     color = m$color) +
        ggplot2::theme_minimal(base_size = ts) +
        ggplot2::theme(
          plot.title   = ggplot2::element_text(face = "bold", hjust = 0.5),
          strip.text   = ggplot2::element_text(face = "bold", size = ts),
          axis.title   = ggplot2::element_text(size = ts),
          axis.text    = ggplot2::element_text(size = ts - 2),
          legend.title = ggplot2::element_text(size = ts),
          legend.text  = ggplot2::element_text(size = ts)
        )

      if (isTRUE(input$bc_show_labels)) {
        p <- p + ggrepel::geom_text_repel(ggplot2::aes(label = name), size = ts / 3.2,
                                          max.overlaps = Inf, show.legend = FALSE)
      }
      p
    }

    output$bc_plot <- renderPlot({ build_bc_plot() })

    output$dl_bc <- downloadHandler(
      filename = function() paste0("PCA_batch_corrected_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
      content = function(file) {
        ggplot2::ggsave(file, plot = build_bc_plot(), width = 12, height = 6, dpi = 300)
      }
    )
  })
}
