# =============================================================================
# de_cross_comparison_module.R  (goes in modules/, sourced from global.R after
#                                de_comparison_log.R)
# -----------------------------------------------------------------------------
# Cross-comparison meta-analysis over EVERY DE contrast logged this session.
#
# WHY THIS EXISTS / WHAT IT READS:
#   de_comparison_log.R already auto-captures each DE run and returns a
#   reactive `comparisons` (list of entries, each with $label, $result (the full
#   DE table), $alpha, $lfc). Nothing consumed that log for analysis — only for
#   export. This module plugs into that SAME reactive (zero changes to the
#   capture logic) and turns the already-collected data into two views iDEP
#   handles weakly:
#     1. UpSet plot  — how the significant-gene SETS across contrasts overlap
#                      (shared vs contrast-unique DEGs).
#     2. Rank-rank   — for a pair of contrasts, are effect sizes concordant?
#                      (a signed -log10(p) x -log10(p) / log2FC x log2FC scatter,
#                       the lightweight cousin of an RRHO plot.)
#
# NO NEW HARD DEPENDENCY:
#   UpSet uses ComplexHeatmap::make_comb_mat()/UpSet() (already used by the main
#   heatmap). Rank-rank uses ggplot2 (already used everywhere). No UpSetR, no
#   RRHO package required.
#
# Wiring (see app.r instructions):
#   UI:     deCrossComparisonUI("de_cross")
#   Server: deCrossComparisonServer("de_cross", comparisons = de_comp_log$comparisons)
#           where de_comp_log is the value returned by deComparisonLogServer().
# =============================================================================


deCrossComparisonUI <- function(id) {
  ns <- NS(id)
  tagList(
    h4("\U0001f500 Cross-comparison meta-analysis (this session)"),
    helpText("Synthesizes every contrast in the DE comparison log. ",
             "Run at least two contrasts, then explore how their results agree."),
    tabsetPanel(
      # ---- UpSet ---------------------------------------------------------
      tabPanel(
        "DEG overlap (UpSet)",
        br(),
        fluidRow(
          column(4, selectInput(ns("upset_dir"), "Genes to include",
                                choices = c("Up- and down-regulated" = "both",
                                            "Up-regulated only"       = "up",
                                            "Down-regulated only"     = "down"),
                                selected = "both")),
          column(4, numericInput(ns("upset_min_size"), "Min intersection size",
                                 value = 1, min = 1, step = 1)),
          column(4, br(), downloadButton(ns("dl_upset"), "\u2b07 PNG",
                                         class = "btn-sm btn-success"))
        ),
        helpText("Each column is a combination of contrasts; bar height = number ",
                 "of genes significant in exactly that combination."),
        plotOutput(ns("upset_plot"), height = "550px"),
        br(),
        downloadButton(ns("dl_overlap_csv"), "\u2b07 Download membership table (CSV)",
                       class = "btn-sm btn-info")
      ),
      # ---- Rank-rank concordance ----------------------------------------
      tabPanel(
        "Concordance (rank-rank)",
        br(),
        fluidRow(
          column(4, uiOutput(ns("rr_x_ui"))),
          column(4, uiOutput(ns("rr_y_ui"))),
          column(4, selectInput(ns("rr_metric"), "Compare by",
                                choices = c("log2 fold change" = "log2fc",
                                            "signed -log10(p)"  = "signed_p"),
                                selected = "log2fc"))
        ),
        fluidRow(
          column(4, checkboxInput(ns("rr_sig_only"),
                                  "Only genes significant in at least one", value = TRUE)),
          column(4, br(), downloadButton(ns("dl_rr"), "\u2b07 PNG",
                                         class = "btn-sm btn-success"))
        ),
        plotOutput(ns("rr_plot"), height = "550px"),
        uiOutput(ns("rr_stats"))
      )
    )
  )
}


deCrossComparisonServer <- function(id, comparisons) {
  moduleServer(id, function(input, output, session) {

    # ---- shared helpers ----------------------------------------------------

    # Significant-gene vector for one log entry, honoring its own alpha/lfc and
    # the requested direction. Uses the SAME significance rule the log itself
    # uses when it computes n_up / n_down, so the sets are consistent with the
    # counts the user already saw in the log table.
    sig_genes <- function(entry, direction = "both") {
      r <- entry$result
      if (is.null(r) || !all(c("gene_id", "log2FoldChange", "padj") %in% colnames(r)))
        return(character(0))
      alpha <- entry$alpha %||% 0.05
      lfc   <- entry$lfc %||% 0
      keep  <- !is.na(r$padj) & r$padj <= alpha
      if (!is.null(lfc) && lfc > 0) keep <- keep & abs(r$log2FoldChange) >= lfc
      if (identical(direction, "up"))   keep <- keep & r$log2FoldChange > 0
      if (identical(direction, "down")) keep <- keep & r$log2FoldChange < 0
      unique(as.character(r$gene_id[keep]))
    }

    # Unique, human-friendly labels (the log allows duplicate labels if the same
    # contrast is run twice; disambiguate so UpSet/pickers don't collapse them).
    entry_labels <- reactive({
      h <- comparisons()
      if (length(h) == 0) return(character(0))
      raw <- vapply(h, function(e) e$label %||% "contrast", "")
      make.unique(raw, sep = " #")
    })

    # ---- UpSet -------------------------------------------------------------
    sig_sets <- reactive({
      h <- comparisons()
      labs <- entry_labels()
      if (length(h) < 2) return(NULL)
      sets <- lapply(h, sig_genes, direction = input$upset_dir %||% "both")
      names(sets) <- labs
      sets <- sets[vapply(sets, length, 0L) > 0]      # drop empty contrasts
      if (length(sets) < 2) return(NULL)
      sets
    })

    build_upset <- function() {
      sets <- sig_sets()
      shiny::validate(shiny::need(!is.null(sets),
                    "Need at least two contrasts that each have >=1 significant gene."))
      cm <- ComplexHeatmap::make_comb_mat(sets, mode = "distinct")
      min_sz <- input$upset_min_size %||% 1
      cm <- cm[ComplexHeatmap::comb_size(cm) >= min_sz]
      shiny::validate(shiny::need(length(ComplexHeatmap::comb_size(cm)) > 0,
                    "No intersections meet the minimum size. Lower the threshold."))
      ComplexHeatmap::UpSet(
        cm,
        comb_order = order(ComplexHeatmap::comb_size(cm), decreasing = TRUE),
        top_annotation = ComplexHeatmap::upset_top_annotation(cm, add_numbers = TRUE),
        right_annotation = ComplexHeatmap::upset_right_annotation(cm, add_numbers = TRUE)
      )
    }

    output$upset_plot <- renderPlot({ ComplexHeatmap::draw(build_upset()) })

    output$dl_upset <- downloadHandler(
      filename = function() paste0("DEG_upset_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
      content = function(file) {
        grDevices::png(file, width = 10, height = 6, units = "in", res = 300)
        ComplexHeatmap::draw(build_upset())
        grDevices::dev.off()
      }
    )

    # Long membership table: gene x contrast (1/0), for download / downstream use.
    membership_df <- reactive({
      sets <- sig_sets()
      req(!is.null(sets))
      all_genes <- sort(unique(unlist(sets)))
      mat <- vapply(sets, function(s) all_genes %in% s, logical(length(all_genes)))
      out <- data.frame(gene_id = all_genes, mat, check.names = FALSE)
      out$n_contrasts <- rowSums(mat)
      out[order(-out$n_contrasts), ]
    })

    output$dl_overlap_csv <- downloadHandler(
      filename = function() paste0("DEG_membership_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
      content = function(file) {
        df <- tryCatch(membership_df(), error = function(e) NULL)
        write.csv(if (is.null(df)) data.frame() else df, file, row.names = FALSE)
      }
    )

    # ---- Rank-rank concordance --------------------------------------------
    output$rr_x_ui <- renderUI({
      labs <- entry_labels()
      selectInput(session$ns("rr_x"), "X-axis contrast",
                  choices = setNames(seq_along(labs), labs),
                  selected = if (length(labs) >= 1) 1 else NULL)
    })
    output$rr_y_ui <- renderUI({
      labs <- entry_labels()
      selectInput(session$ns("rr_y"), "Y-axis contrast",
                  choices = setNames(seq_along(labs), labs),
                  selected = if (length(labs) >= 2) 2 else NULL)
    })

    # Merge two contrasts on gene_id, producing the chosen metric for each.
    rr_data <- reactive({
      h <- comparisons()
      labs <- entry_labels()
      req(length(h) >= 2, input$rr_x, input$rr_y)
      ix <- as.integer(input$rr_x); iy <- as.integer(input$rr_y)
      req(ix >= 1, iy >= 1, ix <= length(h), iy <= length(h))
      shiny::validate(shiny::need(ix != iy, "Pick two different contrasts."))

      metric_col <- function(r) {
        if (identical(input$rr_metric %||% "log2fc", "signed_p")) {
          p <- r$pvalue
          if (is.null(p)) p <- r$padj
          pos <- p[is.finite(p) & p > 0]
          floor_p <- if (length(pos)) min(pos) else .Machine$double.xmin
          p[!is.finite(p) | p <= 0] <- floor_p
          sign(r$log2FoldChange) * -log10(p)
        } else {
          r$log2FoldChange
        }
      }

      rx <- h[[ix]]$result; ry <- h[[iy]]$result
      req(!is.null(rx), !is.null(ry))
      dx <- data.frame(gene_id = as.character(rx$gene_id), x = metric_col(rx),
                       x_sig = !is.na(rx$padj) & rx$padj <= (h[[ix]]$alpha %||% 0.05))
      dy <- data.frame(gene_id = as.character(ry$gene_id), y = metric_col(ry),
                       y_sig = !is.na(ry$padj) & ry$padj <= (h[[iy]]$alpha %||% 0.05))
      m <- merge(dx, dy, by = "gene_id")
      m <- m[is.finite(m$x) & is.finite(m$y), , drop = FALSE]
      if (isTRUE(input$rr_sig_only)) m <- m[m$x_sig | m$y_sig, , drop = FALSE]
      attr(m, "labs") <- c(labs[ix], labs[iy])
      m
    })

    build_rr_plot <- function() {
      m <- rr_data()
      shiny::validate(shiny::need(nrow(m) >= 3, "Fewer than 3 shared genes to compare with these settings."))
      labs <- attr(m, "labs")
      axis_lab <- if (identical(input$rr_metric %||% "log2fc", "signed_p"))
        "signed -log10(p)" else "log2 fold change"

      ggplot2::ggplot(m, ggplot2::aes(x = x, y = y)) +
        ggplot2::geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
        ggplot2::geom_vline(xintercept = 0, color = "grey70", linewidth = 0.3) +
        ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                             color = "grey50") +
        ggplot2::geom_point(alpha = 0.35, size = 1.3, color = "#1F3864") +
        ggplot2::labs(
          title = "Rank-rank concordance",
          subtitle = sprintf("%d shared genes", nrow(m)),
          x = paste0(labs[1], "\n", axis_lab),
          y = paste0(labs[2], "\n", axis_lab)
        ) +
        ggplot2::theme_bw(base_size = 13) +
        ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
    }

    output$rr_plot <- renderPlot({ build_rr_plot() })

    output$rr_stats <- renderUI({
      m <- tryCatch(rr_data(), error = function(e) NULL)
      if (is.null(m) || nrow(m) < 3) return(NULL)
      rho <- suppressWarnings(stats::cor(m$x, m$y, method = "spearman",
                                         use = "complete.obs"))
      concordant <- mean(sign(m$x) == sign(m$y), na.rm = TRUE)
      helpText(sprintf(
        "Spearman rho = %.3f across %d shared genes. %.0f%% share the same sign of effect.",
        rho, nrow(m), 100 * concordant))
    })

    output$dl_rr <- downloadHandler(
      filename = function() paste0("rank_rank_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
      content = function(file) {
        ggplot2::ggsave(file, plot = build_rr_plot(), width = 8, height = 7, dpi = 300)
      }
    )
  })
}
