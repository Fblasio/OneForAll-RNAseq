# =============================================================================
# rnaseq_de_module.R
# -----------------------------------------------------------------------------
# Shiny module for differential expression analysis using DESeq2.
#
# Inputs (reactive arguments from rnaseqDdsServer):
#   dds_rv      — fitted DESeqDataSet produced by rnaseqDdsServer
#   metadata_rv — saved metadata table (used to populate contrast level selectors)
#
# Workflow:
#   1. User selects a metadata factor and two levels (test vs. reference).
#   2. On button press, DESeq2::results() extracts statistics for that contrast.
#   3. Filters are applied: log2FC threshold, significance, direction.
#   4. Results are shown in an interactive DT table with export buttons.
#
# The log2 fold change is always expressed as log2(test / reference):
#   positive  → gene is higher in the test condition
#   negative  → gene is higher in the reference condition
# =============================================================================


# -----------------------------------------------------------------------------
# UI
# rnaseqDeUI(): sidebar controls for contrast definition and result filtering.
# -----------------------------------------------------------------------------
rnaseqDeUI <- function(id, start_step = 5) {
ns <- NS(id)

# Continues the numbering that rnaseqDdsUI() started; see the note there on why
# the offset is passed in rather than hard-coded.
s_run    <- start_step         # 5
s_filter <- start_step + 1L    # 6

tagList(
step_header(s_run, "Run a contrast"),

# --- Dynamic contrast selectors (rendered by the server) -------------------
# Populated once metadata is available. Contains:
#   contrast_factor  — which metadata column defines the groups
#   contrast_num     — test level (numerator of log2FC)
#   contrast_den     — reference level (denominator of log2FC)

# Result type is shown only when the fitted model contains interaction terms.
uiOutput(ns("result_mode_ui")),

# Dynamic result controls.
# Simple contrasts use factor/test/reference selectors.
# Interaction mode exposes fitted interaction coefficients only.
uiOutput(ns("contrast_ui")),
uiOutput(ns("interaction_ui")),
uiOutput(ns("interaction_reference_info")),

# --- Model settings (only take effect on Run) --------------------------------
# alpha_cutoff and apply_shrinkage are consumed inside de_results_all(), which
# is an eventReactive on run_de. Editing them without clicking Run changes
# nothing, so they are grouped here rather than among the filters below, which
# DO apply instantly. See the alpha-drift guard in de_results_filtered().
tags$strong("Model settings"),
helpText("Used when the model is evaluated — click Run after changing them."),
numericInput(ns("alpha_cutoff"), "Adjusted p-value cutoff", value = 0.05, min = 0, max = 1, step = 0.01),
helpText("Also the target FDR for DESeq2's independent filtering, so the padj column itself depends on this value."),
checkboxInput(ns("apply_shrinkage"),
# Deliberately says what the helpText below does NOT: that significance is
# untouched, and what the two settings are actually FOR. The two are meant to
# complement each other rather than repeat.
lbl("Apply LFC shrinkage (ashr)", "Pulls the fold changes of low-count, noisy genes toward zero so log2FC is comparable across genes. It does NOT change p-values or padj - only the effect size. Leave it on for volcano plots, ranking and GSEA."),
value = TRUE),
helpText("Shrinks log2FC for low-count/high-variance genes toward zero — recommended for ranking and the volcano plot. Uncheck to see raw DESeq2 MLE log2FC instead."),

hr(),

# --- Result filters (applied instantly to the last Run) ----------------------
step_header(s_filter, "Filter results", optional = TRUE),
helpText("Applied to the table immediately; no re-run needed."),
numericInput(ns("lfc_cutoff"),   "Absolute log2FC cutoff",  value = 0,    min = 0, step = 0.1),
helpText("Applied to the log2FC column as shown, i.e. to the shrunken estimate when LFC shrinkage is on."),

# sig_only: restrict to genes with padj <= alpha_cutoff
checkboxInput(ns("sig_only"),   "Show only significant genes",    value = TRUE),
# up_only / down_only: mutually exclusive direction filters
checkboxInput(ns("up_only"),    "Show only upregulated genes",    value = FALSE),
checkboxInput(ns("down_only"),  "Show only downregulated genes",  value = FALSE),

h4("\U0001f4dc Session Log"),
appLogUI("applog_rnaseq"),

# ── Normalised counts ──────────────────────────────────────────────────────
# Exported because plotting a gene across groups is the first thing asked of a
# DE result, and until now that matrix could only be reached from outside the
# application. Deliberately NOT added to the reproducibility bundle: these are
# the user's data on a derived scale, and the bundle carries no data.
downloadButton(ns("dl_norm_counts"), "Normalised counts (CSV)",
               class = "btn-outline-secondary", style = "width:100%;"),
helpText("counts(dds, normalized = TRUE) for every gene that passed the pre-filter, on the design currently fitted."),

hr(),
# Pinned to the bottom of the viewport while the sidebar scrolls, so re-running
# never requires scrolling back past the contrast selectors, the model settings
# and the filters above it. Same mechanism as "Generate Heatmap": the
# .cheski-sticky-run class is declared once in app.r's tags$head and is
# position:sticky, so it needs no JS and renders in place when the sidebar is
# shorter than the window.
#
# Placed LAST on purpose. position:sticky pins an element within its parent's
# box, so anything rendered after it would scroll UNDER the bar -- which is why
# the Session Log moved above rather than the button staying where it was.
#
# NOTE for anyone adding a wrapper here: sticky is cancelled by any ancestor
# with overflow hidden/auto, so the sidebar must stay free of those.
div(class = "cheski-sticky-run",
    actionButton(ns("run_de"), "Run differential analysis",
                 icon  = icon("play"),
                 class = "btn-success btn-lg",
                 width = "100%"))
)
}


# -----------------------------------------------------------------------------
# SERVER
# rnaseqDeServer(): handles contrast selection, DE computation, and table output.
# -----------------------------------------------------------------------------
# species_rv: OPTIONAL reactive returning the active organism. Supplied by app.r
# so the staleness guard in de_sig_genes() can actually fire -- see the comment
# there. NULL-defaulted, so any other caller keeps the previous behaviour.
rnaseqDeServer <- function(id, dds_rv, metadata_rv, species_rv = NULL) {
moduleServer(id, function(input, output, session) {
  freezeReactiveValue(input, "result_mode")
  freezeReactiveValue(input, "contrast_factor")
  freezeReactiveValue(input, "contrast_num")
  freezeReactiveValue(input, "contrast_den")
  freezeReactiveValue(input, "interaction_name")
  # ---------------------------------------------------------------------------
  # output$contrast_ui: renders the factor selector and a placeholder for the
  # level selectors. Depends on the saved metadata being available.
  # Excludes the "sample" column since it is the row ID, not a grouping factor.
  # ---------------------------------------------------------------------------
  # Render simple contrast selectors only when the standard contrast mode is active.
  
  # Detect whether the fitted DESeq2 object exposes interaction coefficients.
  # The decision is based on resultsNames(dds), which reflects the fitted model.
  has_interaction_terms <- reactive({
    req(dds_rv())
    
    coef_names <- DESeq2::resultsNames(dds_rv())
    interaction_names <- coef_names[
      !coef_names %in% "Intercept" &
        !grepl("_vs_", coef_names) &
        grepl("[.:]", coef_names)
    ]
    
    length(interaction_names) > 0
  })
  

  # Render the result mode selector only when interaction coefficients exist.
  # For one-factor designs, the selector is hidden because only simple contrasts apply.
  output$result_mode_ui <- renderUI({
    req(dds_rv())
    
    if (!isTRUE(has_interaction_terms())) {
      return(NULL)
    }
    
    selectInput(
      session$ns("result_mode"),
      "Result type",
      choices = c(
        "Simple contrast" = "contrast",
        "Interaction term" = "interaction"
      ),
      selected = "contrast",
      width = "100%"
    )
  })
  
  # Reset DE input state whenever a new fitted DESeq2 object is received.
  # The DE controls are marked as not ready until the new choices are rebuilt.
  observeEvent(dds_rv(), {
    req(dds_rv())
    
    meta <- as.data.frame(SummarizedExperiment::colData(dds_rv()))
    design_vars <- all.vars(DESeq2::design(dds_rv()))
    valid_factors <- intersect(design_vars, colnames(meta))
    valid_factors <- valid_factors[sapply(valid_factors, function(x) {
      vals <- unique(as.character(meta[[x]]))
      vals <- vals[!is.na(vals) & nzchar(vals)]
      length(vals) >= 2
    })]
    
    default_factor <- if (length(valid_factors) >= 1) valid_factors[1] else NULL
    default_levels <- if (!is.null(default_factor)) {
      vals <- unique(as.character(meta[[default_factor]]))
      vals[!is.na(vals) & nzchar(vals)]
    } else {
      character(0)
    }
    
    coef_names <- DESeq2::resultsNames(dds_rv())
    interaction_names <- coef_names[
      !coef_names %in% "Intercept" &
        !grepl("_vs_", coef_names) &
        grepl("[.:]", coef_names)
    ]
    
    if (isTRUE(has_interaction_terms())) {
      updateSelectInput(
        session, "result_mode",
        choices = c(
          "Main effect / simple contrast" = "contrast",
          "Interaction term" = "interaction"
        ),
        selected = "contrast"
      )
    } else {
      updateSelectInput(
        session, "result_mode",
        choices = c(
          "Main effect / simple contrast" = "contrast"
        ),
        selected = "contrast"
      )
      
      updateSelectInput(
        session, "interaction_name",
        choices = character(0),
        selected = character(0)
      )
    }
    
    if (!is.null(default_factor)) {
      updateSelectInput(
        session, "contrast_factor",
        choices = valid_factors,
        selected = default_factor
      )
    }
    
    if (length(default_levels) >= 2) {
      updateSelectInput(
        session, "contrast_num",
        choices = default_levels,
        selected = default_levels[1]
      )
      updateSelectInput(
        session, "contrast_den",
        choices = default_levels,
        selected = default_levels[2]
      )
    }
    
    if (length(interaction_names) >= 1) {
      updateSelectInput(
        session, "interaction_name",
        choices = interaction_names,
        selected = interaction_names[1]
      )
    }
  }, ignoreInit = FALSE)
  
  # Render simple contrast controls when no interaction is selected.
  # Factor names are taken directly from colData(dds) to guarantee exact match
  # with the names used internally by DESeq2::results().
  output$contrast_ui <- renderUI({
    req(dds_rv())
    
    meta <- as.data.frame(SummarizedExperiment::colData(dds_rv()))
    design_vars <- all.vars(DESeq2::design(dds_rv()))
    
    cols <- intersect(design_vars, colnames(meta))
    cols <- cols[sapply(cols, function(x) {
      vals <- unique(as.character(meta[[x]]))
      vals <- vals[!is.na(vals) & nzchar(vals)]
      length(vals) >= 2
    })]
    
    if (length(cols) == 0) return(NULL)
    
    current_mode <- input$result_mode %||% "contrast"
    show_interaction_mode <- isTRUE(has_interaction_terms())
    
    if (show_interaction_mode && identical(current_mode, "interaction")) {
      return(NULL)
    }
    
    tagList(
      selectInput(
        session$ns("contrast_factor"),
        "Factor",
        choices = cols,
        selected = cols[1],
        width = "100%"
      ),
      uiOutput(session$ns("contrast_levels_ui"))
    )
  })

  # ---------------------------------------------------------------------------
  # output$contrast_levels_ui: renders the Treatment and Control level selectors.
  # Re-renders whenever contrast_factor changes to reflect the available levels
  # for the newly selected column.
  # NA and empty strings are excluded from the level list.
  # Requires at least 2 distinct levels to proceed.
  # ---------------------------------------------------------------------------
  # Render the Treatment and Control level selectors using colData(dds) levels.
  # This ensures level names match exactly what DESeq2::results() expects.
  output$contrast_levels_ui <- renderUI({
    req(dds_rv(), input$contrast_factor)
    
    meta <- as.data.frame(SummarizedExperiment::colData(dds_rv()))
    req(input$contrast_factor %in% colnames(meta))
    
    lvls <- levels(meta[[input$contrast_factor]])
    if (is.null(lvls)) lvls <- unique(as.character(meta[[input$contrast_factor]]))
    lvls <- lvls[!is.na(lvls) & nzchar(lvls)]
    req(length(lvls) >= 2)
    
    # Default Control = the design's reference level (lvls[1]); Treatment = next.
    # The input IDs stay contrast_num / contrast_den: they are read by
    # de_results_all(), by contrast_params() and by the reproducibility bundle,
    # and DESeq2's own vocabulary for the pair is numerator/denominator. Only
    # the LABELS change, to the words this audience actually uses.
    tagList(
      selectInput(session$ns("contrast_num"),
                  lbl("Treatment level", "The condition being tested - the numerator of the fold change. log2FC > 0 means higher in Treatment than in Control. Swapping Treatment and Control flips every sign in the table, and nothing will warn you."),
                  choices = lvls, selected = lvls[2], width = "100%"),
      selectInput(session$ns("contrast_den"),
                  lbl("Control level", "The baseline being compared against - normally control, mock or wild-type. 'Up-regulated' always means up relative to this. Defaults to the control level the model was fitted with."),
                  choices = lvls, selected = lvls[1], width = "100%")
    )
  })

  # ---------------------------------------------------------------------------
  # de_results_filtered: eventReactive triggered by the "Run" button.
  # ignoreInit = TRUE: only on explicit button press.
  #
  # Steps:
  #   1. Build contrast vector: c(factor, test_level, reference_level).
  #   2. Call DESeq2::results() with that contrast.
  #        log2FC = log2(mean_test / mean_reference)
  #   3. Convert to data frame; move gene IDs from rownames to a column.
  #   4. Drop genes with NA log2FoldChange (typically low-count genes
  #      that were not tested due to the independent filtering step).
  #   5. Apply absolute log2FC threshold (default: |log2FC| >= 0).
  #   6. Optionally filter to significant genes (padj <= alpha_cutoff).
  #   7. Optionally restrict to upregulated or downregulated genes.
  #        Note: if both up_only and down_only are checked, neither applies.
  #   8. Sort by padj ascending (NAs pushed to the bottom).
  # ---------------------------------------------------------------------------
  # Render interaction controls only when the fitted model contains interactions
  # and the interaction mode is explicitly selected.
  
  output$interaction_ui <- renderUI({
    req(dds_rv())
    
    current_mode <- input$result_mode %||% "contrast"
    if (!isTRUE(has_interaction_terms()) || !identical(current_mode, "interaction")) {
      return(NULL)
    }
    
    coef_names <- DESeq2::resultsNames(dds_rv())
    interaction_names <- coef_names[
      !coef_names %in% "Intercept" &
        !grepl("_vs_", coef_names) &
        grepl("[.:]", coef_names)
    ]
    
    if (length(interaction_names) == 0) {
      return(NULL)
    }
    
    tagList(
      selectInput(
        session$ns("interaction_name"),
        "Interaction coefficient",
        choices = interaction_names,
        selected = interaction_names[1],
        width = "100%"
      )
    )
  })
 
  # Show reference levels only when interaction mode is active.
  # This avoids displaying interaction-specific information for one-factor designs.
  output$interaction_reference_info <- renderUI({
    req(dds_rv())
    
    current_mode <- input$result_mode %||% "contrast"
    if (!isTRUE(has_interaction_terms()) || !identical(current_mode, "interaction")) {
      return(NULL)
    }
    
    meta_df <- as.data.frame(SummarizedExperiment::colData(dds_rv()))
    design_vars <- all.vars(DESeq2::design(dds_rv()))
    design_vars <- intersect(design_vars, colnames(meta_df))
    if (length(design_vars) == 0) return(NULL)
    
    ref_lines <- lapply(design_vars, function(var_name) {
      var_data <- meta_df[[var_name]]
      
      if (!is.factor(var_data)) {
        return(tags$li(sprintf("%s: not a factor", var_name)))
      }
      
      tags$li(sprintf("%s: control = %s", var_name, levels(var_data)[1]))
    })

    div(
      class = "alert alert-info",
      tags$b("Active control levels"),
      tags$ul(ref_lines),
      helpText(
        "To change interaction interpretation, change the control level of a factor and rebuild the DESeq2 object."
      )
    )
  })
  
  # Compute DE results according to the selected result mode.
  # When no interaction selector is available, the mode defaults to "contrast".
  # Compute the full DE result table once per button click.
  # This object keeps all tested genes and serves as the source for the volcano plot.
  de_results_all <- eventReactive(input$run_de, {
    
    req(dds_rv())
    current_mode <- input$result_mode %||% "contrast"
    
    # Force simple contrast mode when the current model has no interaction terms.
    if (!isTRUE(has_interaction_terms())) {
      current_mode <- "contrast"
    }
    
    message("=== DE trigger ===")
    message("current_mode : ", current_mode)
    message("contrast_factor: ", input$contrast_factor %||% "NULL")
    message("contrast_num : ", input$contrast_num %||% "NULL")
    message("contrast_den : ", input$contrast_den %||% "NULL")
    
    if (identical(current_mode, "interaction")) {
      req(input$interaction_name)
      res <- DESeq2::results(dds_rv(), name = input$interaction_name, alpha = input$alpha_cutoff %||% 0.05)
    } else {
      req(input$contrast_factor, input$contrast_num, input$contrast_den)
      
      meta <- as.data.frame(SummarizedExperiment::colData(dds_rv()))
      req(input$contrast_factor %in% colnames(meta))
      
      valid_levels <- unique(as.character(meta[[input$contrast_factor]]))
      valid_levels <- valid_levels[!is.na(valid_levels) & nzchar(valid_levels)]
      
      req(input$contrast_num %in% valid_levels)
      req(input$contrast_den %in% valid_levels)
      req(input$contrast_num != input$contrast_den)
      
      contrast_vec <- c(input$contrast_factor, input$contrast_num, input$contrast_den)
      
      res <- tryCatch(
        DESeq2::results(dds_rv(), contrast = contrast_vec, alpha = input$alpha_cutoff %||% 0.05),
        error = function(e) {
          shiny::showNotification(
            paste(
              "DESeq2 could not extract the selected contrast:",
              conditionMessage(e)
            ),
            type = "error",
            duration = 10
          )
          NULL
        }
      )
      
      if (is.null(res)) {
        return(NULL)
      }
    }
    
    if (is.null(res)) {
      return(NULL)
    }
    
    # Wald statistic (log2FC / lfcSE) on the unshrunken MLE. lfcShrink() drops
    # this column, so it is captured here before the shrinkage block below.
    stat_mle <- if ("stat" %in% colnames(res)) {
      setNames(res$stat, rownames(res))
    } else NULL

    # Make shrinkage to the results
    # Track whether shrinkage was actually applied to this result, so
    # downstream consumers (table caption, volcano title, comparison log)
    # can say which kind of log2FC they're showing.
    shrinkage_applied <- FALSE
    
    # LFC shrinkage — user-controlled. Only runs when the checkbox is on.
    if (isTRUE(input$apply_shrinkage %||% TRUE)) {
      res <- tryCatch({
        # ashr is the longest step in this tab on a full transcriptome and had
        # no indicator at all. with_progress_messages() (log_utils.R) surfaces
        # lfcShrink()'s own message while it runs; the call itself is unchanged.
        shrunk <- with_progress_messages(
          if (identical(current_mode, "interaction")) {
            DESeq2::lfcShrink(dds_rv(), coef = input$interaction_name,
                              type = "ashr", res = res)
          } else {
            DESeq2::lfcShrink(dds_rv(), contrast = contrast_vec,
                              type = "ashr", res = res)
          },
          "Shrinking log2 fold changes (ashr)")
        shrinkage_applied <- TRUE
        shrunk
      }, error = function(e) {
        showNotification(paste("LFC shrinkage skipped:", conditionMessage(e)),
                         type = "warning", duration = 6)
        res   # fall back to unshrunken MLE
      })
    }
    
    # Convert result object to a regular data frame and expose gene IDs explicitly.
    res_df <- as.data.frame(res)
    
    # rownames(dds) are already canonical tomato IDs at this point, because
    # the counts matrix was canonicalized in import_raw_counts()/
    # import_salmon_counts() before the DESeqDataSet was built.
    res_df$gene_id <- rownames(res_df)
    # Restore the MLE Wald statistic removed by lfcShrink(). Matched by name
    # rather than relying on row order.
    if (!is.null(stat_mle)) res_df$stat <- unname(stat_mle[res_df$gene_id])
    attr(res_df, "shrinkage_applied") <- shrinkage_applied
    rownames(res_df) <- NULL

    res_df <- res_df[, c("gene_id", setdiff(colnames(res_df), "gene_id")), drop = FALSE]
    
    # Remove rows where log2 fold change was not estimated.
    # These genes are not informative for ranking or plotting.
    res_df <- res_df[!is.na(res_df$log2FoldChange), , drop = FALSE]
    
    
    # Sort by adjusted p-value, sending missing values to the bottom.
    out <- res_df[order(res_df$padj, na.last = TRUE), , drop = FALSE]
    attr(out, "shrinkage_applied") <- shrinkage_applied
    # Record the alpha this padj column was actually computed with. DESeq2's
    # independent filtering picks the baseMean quantile that maximises the
    # number of genes with padj < alpha, so padj is specific to this value.
    attr(out, "alpha_used") <- input$alpha_cutoff %||% 0.05
    # Record which species this contrast was computed under, so the staleness
    # check in de_sig_genes() can ask the result directly instead of guessing
    # from the shape of the gene IDs. Same mechanism as the two attributes above.
    attr(out, "species") <- ACTIVE_SPECIES
    out
  }, ignoreInit = TRUE)
  
  # Derive the filtered DE table from the complete result table.
  # This filtered view is used for the interactive table and DEG-based downstream steps.
  de_results_filtered <- reactive({
    res_df <- de_results_all()
    req(!is.null(res_df))
    
    # Alpha-drift guard. de_results_all() only fires on run_de, but this
    # reactive re-reads input$alpha_cutoff every time it changes. Without this,
    # typing a new alpha silently re-filters a padj column computed under the
    # OLD alpha's independent-filtering threshold: the table updates, looks
    # correct, and is wrong. Warn instead of recomputing, so no lfcShrink()
    # runs on a keystroke.
    alpha_used <- attr(res_df, "alpha_used")
    alpha_now  <- suppressWarnings(as.numeric(input$alpha_cutoff %||% 0.05))
    if (!is.null(alpha_used) && !is.na(alpha_used) && !is.na(alpha_now) &&
        !isTRUE(all.equal(as.numeric(alpha_used), alpha_now))) {
      showNotification(
        sprintf(paste0("\u26a0\ufe0f padj was computed with alpha = %s (independent ",
                       "filtering is tuned to that target). You are filtering at %s ",
                       "\u2014 click 'Run differential analysis' to recompute."),
                format(alpha_used), format(alpha_now)),
        type = "warning", duration = NULL, id = "alpha_drift")
    } else {
      removeNotification("alpha_drift")
    }
    
    # Filter thresholds, coerced to finite values before ANY subsetting below.
    # Clearing a numericInput sends NA, and `%||%` does not catch that: it only
    # replaces NULL / length-0 / "" (global.R:42-46). Subsetting by an all-NA
    # logical index does NOT return zero rows -- R returns one all-NA row per
    # input row, so the table silently filled with <NA> instead of erroring or
    # emptying. Falls back to each widget's own default (0 / 0.05) so the filter
    # matches what the box shows once it is repopulated. Same is.na-guard
    # pattern already used for top_var_n and universe_threshold in app.r.
    # Deliberately NOT reusing `alpha_now` from the drift guard above: that
    # guard must keep seeing the raw NA and stay silent, rather than compare a
    # substituted 0.05 against alpha_used and warn about a drift the user never
    # caused.
    lfc_cut <- suppressWarnings(as.numeric(input$lfc_cutoff %||% 0))
    if (!is.finite(lfc_cut)) lfc_cut <- 0
    alpha_cut <- suppressWarnings(as.numeric(input$alpha_cutoff %||% 0.05))
    if (!is.finite(alpha_cut)) alpha_cut <- 0.05

    # Apply minimum absolute log2 fold change threshold.
    # NOTE: this is a post-hoc filter on the log2FC column as displayed, i.e. on
    # the ashr-shrunken estimate when shrinkage is on. padj still tests
    # H0: log2FC = 0, so a gene passing both is one that is significantly
    # non-zero AND whose (shrunken) estimate clears the cutoff. Shrinkage is
    # what makes this defensible: it already pulls in the inflated log2FC of
    # low-count genes that a filter on the raw MLE would wrongly admit.
    res_df <- res_df[abs(res_df$log2FoldChange) >= lfc_cut, , drop = FALSE]

    # Keep only statistically significant genes when requested.
    if (isTRUE(input$sig_only)) {
      res_df <- res_df[!is.na(res_df$padj) & res_df$padj <= alpha_cut, , drop = FALSE]
    }

    # Apply one-direction filter only when exactly one direction is selected.
    if (isTRUE(input$up_only) && !isTRUE(input$down_only)) {
      res_df <- res_df[res_df$log2FoldChange > lfc_cut, , drop = FALSE]
    }

    if (isTRUE(input$down_only) && !isTRUE(input$up_only)) {
      res_df <- res_df[res_df$log2FoldChange < (-lfc_cut), , drop = FALSE]
    }
    
    # Order by adjusted p-value, sending missing values to the bottom.
    res_df[order(res_df$padj, na.last = TRUE), , drop = FALSE]
  })

  # ---------------------------------------------------------------------------
  # output$de_table: interactive DT table of the filtered DE results.
  # Features:
  #   - Column-level filters in the header row (filter = "top").
  #   - Export buttons: copy to clipboard, CSV, Excel.
  #   - Horizontal scrolling for wide tables.
  #   - 25 rows per page as a reasonable default for gene tables.
  # ---------------------------------------------------------------------------
  # Builds species-aware gene links for the DEG table.
  .make_deg_gene_link <- function(gid) {
    # gid is already canonical (Solyc##g######.<version>); only the
    # SGN link identifier (without version) needs to be derived here.
    gid <- trimws(as.character(gid))
    if (identical(ACTIVE_SPECIES, "tomato")) {
      gid_link <- stringr::str_extract(gid, "Solyc[0-9]{2}g[0-9]+")
      if (is.na(gid_link)) gid_link <- gid
      url <- paste0(GENE_DB_URL, gid_link)
      paste0("<a href='", url, "' target='_blank' rel='noopener noreferrer'>", gid, "</a>")
    } else {
      paste0("<a href='https://bar.utoronto.ca/thalemine/portal.do?externalids=",
             toupper(gid), "' target='_blank' rel='noopener noreferrer'>", gid, "</a>")
    }
  }
  
  
  output$de_table <- DT::renderDT({
    res_tbl <- de_results_filtered()
    req(!is.null(res_tbl))
    
    restbl_display <- res_tbl
    
    restbl_display$gene_id <- ifelse(
      is.na(restbl_display$gene_id) | !nzchar(restbl_display$gene_id),
      restbl_display$gene_id,
      vapply(restbl_display$gene_id, .make_deg_gene_link, character(1))
    )
    
    shrink_note <- if (isTRUE(attr(res_tbl, "shrinkage_applied"))) {
      "log2FC values are ashr-shrunken."
    } else {
      "log2FC values are raw DESeq2 MLE (no shrinkage applied)."
    }
    message("DE table: ", shrink_note)
    
    DT::datatable(
      restbl_display,
      filter = "top",
      rownames = FALSE,
      escape = FALSE,
      extensions = c("Buttons"),
      options = list(
        scrollX = TRUE,
        pageLength = 25,
        dom = "Bfrtip",
        buttons = list(
          list(
            extend = "copy",
            title = NULL,
            header = FALSE,
            exportOptions = list(
              columns = 0,
              modifier = list(
                page = "all",
                selected = NULL
              )
            )
          ),
          list(
            extend = "csv",
            exportOptions = list(
              modifier = list(
                page = "all",
                selected = NULL
              )
            )
          ),
          list(
            extend = "excel",
            exportOptions = list(
              modifier = list(
                page = "all",
                selected = NULL
              )
            )
          )
        )
      )
    )
  }, server = FALSE)
  
  # Reactive exposing only the gene_id column of the current filtered DE table.
  # Used by selectionHistoryServer to let the user save named gene lists.
  de_sig_genes <- reactive({
    tbl <- de_results_filtered()
    if (is.null(tbl) || nrow(tbl) == 0) return(NULL)
    
    # Staleness guard, asked of the RESULT rather than of the gene IDs.
    # Switching species calls rnaseq_handles$reset() (app.r:375) and nulls the
    # dds, but de_results_all() is an eventReactive and keeps serving the
    # previous species' result -- enriching that against the new annotation
    # would be wrong. Only the result can say which species produced it; a gene
    # ID cannot, which is why the ID-shape test this replaces was unsound.
    #
    # What it replaces: a per-gene regex (^AT[0-9MC]G[0-9]+ /
    # ^Solyc[0-9]{2}g[0-9]+) that DROPPED every row not matching it. That
    # removed legitimate rows of a perfectly correct matrix -- transposable-
    # element genes (AT1TE...), novel StringTie loci (MSTRG...), spike-ins,
    # miRNAs -- and, for tomato, every version-less or uppercase ID, because
    # canonicalize_tomato_ids_impl() only rewrites IDs carrying a version
    # suffix and the tomato pattern was case-sensitive. When nothing matched it
    # returned NULL with only a console message(), so the GO/KEGG buttons did
    # nothing at all with no user-facing explanation. Genes are no longer
    # dropped here; IDs the annotation cannot map are reported by the
    # enrichment step itself, where the user can actually see the count.
    # The current organism is read through a REACTIVE, not from the global.
    # ACTIVE_SPECIES is a plain variable that app.r rewrites with
    # assign(..., .GlobalEnv), so rewriting it invalidates nothing: this reactive
    # kept serving its cached value and the guard below never ran. Measured before
    # the fix -- after switching to tomato, de_sig_genes() still handed out 17
    # Arabidopsis IDs and no notification appeared, so GO/KEGG then enriched
    # Arabidopsis IDs against the tomato OrgDb and reported "no terms found"
    # rather than an ID mismatch.
    # species_rv is app.r's own reactiveVal (app.r:712), already updated by the
    # apply_species observer (app.r:750), so depending on it here is what makes
    # the switch invalidate this reactive. Falls back to the global when the
    # argument is absent, which preserves the old behaviour for any other caller.
    active_species <- if (is.function(species_rv)) species_rv() else ACTIVE_SPECIES

    res_species <- attr(de_results_all(), "species")
    if (!is.null(res_species) && !identical(res_species, active_species)) {
      showNotification(
        paste0("⚠️ These DE results were computed for ", res_species,
               ", but the active species is now ", active_species,
               ". Re-run the differential analysis before using them for enrichment."),
        type = "warning", duration = 12, id = "de_species_stale")
      return(NULL)
    }

    # Drop only unusable entries (NA / blank), never well-formed IDs. Kept as
    # an explicit filter rather than relying on an upstream pattern match, so
    # an empty or missing gene_id can never reach canon_ids()/enrichGO().
    genes_out <- as.character(tbl$gene_id)
    genes_out <- genes_out[!is.na(genes_out) & nzchar(trimws(genes_out))]
    if (length(genes_out) == 0) return(NULL)

    genes_out
  })
  
  # ── Normalised count matrix ──────────────────────────────────────────────
  # counts(dds, normalized = TRUE) divides each column by its size factor, so
  # the values stay on a count scale and stay comparable between samples. This
  # is the matrix a per-gene expression plot is drawn from; the VST matrix the
  # heatmap and the PCA use answers a different question and is not what is
  # written here. The header records the design, because the size factors are
  # estimated from the object that design was fitted on.
  output$dl_norm_counts <- downloadHandler(
    filename = function() {
      paste0("normalised_counts_", format(Sys.Date(), "%Y%m%d"), ".csv")
    },
    content = function(file) {
      dds <- dds_rv()
      # Must be fully qualified: global.R loads jsonlite AFTER shiny, and
      # jsonlite::validate() (argument `txt`) masks shiny::validate(), which
      # fails with "is.character(txt) is not TRUE". Same trap as gsea_module.R.
      shiny::validate(shiny::need(
        !is.null(dds), "Fit a model before downloading counts."))
      m <- DESeq2::counts(dds, normalized = TRUE)
      con <- file(file, open = "wt", encoding = "UTF-8")
      on.exit(close(con), add = TRUE)
      writeLines(c(
        "# Size-factor normalised counts, counts(dds, normalized = TRUE).",
        paste0("# Design: ",
               tryCatch(paste(deparse(DESeq2::design(dds)), collapse = ""),
                        error = function(e) "unknown")),
        paste0("# Genes: ", nrow(m), "   Samples: ", ncol(m)),
        "# Size factors are estimated from the genes that passed the pre-filter,",
        "# so the same raw matrix filtered differently gives different values."
      ), con)
      utils::write.csv(data.frame(gene_id = rownames(m), m, check.names = FALSE),
                       con, row.names = FALSE)
    }
  )

  list(
    de_results = reactive(de_results_filtered()),
    de_results_all = reactive(de_results_all()),
    de_sig_genes = de_sig_genes,
    run_trigger = reactive(input$run_de),
    contrast_params = reactive(list(
      mode        = input$result_mode %||% "contrast",
      factor      = input$contrast_factor,
      test        = input$contrast_num,
      reference   = input$contrast_den,
      interaction = input$interaction_name,
      alpha       = input$alpha_cutoff %||% 0.05,
      lfc         = input$lfc_cutoff %||% 0,
      shrinkage   = isTRUE(input$apply_shrinkage %||% TRUE),
      design      = tryCatch(paste(deparse(DESeq2::design(dds_rv())), collapse = ""),
                             error = function(e) NA_character_),
      # ── Result-defining settings that were recorded nowhere ────────────────
      # sig_only / up_only / down_only are applied in de_results_filtered() and
      # therefore decide de_sig_genes() -- the literal gene vector handed to GO
      # and KEGG. Without them the bundle cannot explain where that list came
      # from, even though it stores the list itself.
      sig_only    = isTRUE(input$sig_only),
      up_only     = isTRUE(input$up_only),
      down_only   = isTRUE(input$down_only),
      # Reference level per design factor, DERIVED from the fitted object
      # rather than from the ref_level_* inputs. levels()[1] is by definition
      # what DESeq2 contrasted against, so this cannot drift from what was
      # actually fitted the way a re-read input can.
      reference_levels = tryCatch({
        md <- as.data.frame(SummarizedExperiment::colData(dds_rv()))
        vars <- intersect(all.vars(DESeq2::design(dds_rv())), colnames(md))
        stats::setNames(
          lapply(vars, function(v) if (is.factor(md[[v]])) levels(md[[v]])[1] else NA_character_),
          vars)
      }, error = function(e) NULL),
      # The count pre-filter, read back from the object it was applied to.
      prefilter = tryCatch(S4Vectors::metadata(dds_rv())$cheski_prefilter,
                           error = function(e) NULL)
    ))
  )
})
}
