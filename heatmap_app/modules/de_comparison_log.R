# =============================================================================
# de_comparison_log.R
# -----------------------------------------------------------------------------
# Session-level log of EVERY differential-expression comparison run during a
# session. Each time the user runs a DE contrast, one entry is appended
# automatically (no manual "save" step). Exports the whole session as:
#   - an R script that re-runs each contrast with DESeq2
#   - a JSON manifest (parameters + per-comparison summary)
#   - a ZIP of one CSV per comparison (full result tables)
#
# It stores results passively; it does NOT re-run or alter any DESeq2 logic.
# =============================================================================

deComparisonLogUI <- function(id) {
  ns <- NS(id)
  tagList(
    h4("\U0001f5c3\ufe0f DE comparison log (this session)"),
    helpText("Every comparison you run is recorded here automatically."),
    DT::DTOutput(ns("log_table")),
    br(),
    fluidRow(
      column(6, actionButton(ns("delete_sel"), "\U0001f5d1\ufe0f Delete selected",
                             class = "btn-danger btn-sm", width = "100%")),
      column(6, actionButton(ns("clear_all"), "\U0001f9f9 Clear all",
                             class = "btn-warning btn-sm", width = "100%"))
    ),
    br(),
    tags$b("Export whole session:"),
    fluidRow(
      column(4, downloadButton(ns("dl_r"),    "\u2b07 R script", class = "btn-sm btn-success", width = "100%")),
      column(4, downloadButton(ns("dl_json"), "\u2b07 JSON",     class = "btn-sm btn-info",    width = "100%")),
      column(4, downloadButton(ns("dl_zip"),  "\u2b07 CSVs (zip)", class = "btn-sm btn-primary", width = "100%"))
    )
  )
}

# server args are the reactives already exposed by rnaseqDeServer + the raw
# contrast inputs, passed in from app.r so this module reads no globals.
deComparisonLogServer <- function(id,
                                  run_trigger,        # reactive/eventReactive that fires on Run
                                  de_results_all,     # reactive: full result data.frame
                                  contrast_params,     # reactive: named list of the contrast settings
                                  session_info = NULL) {  # optional reactive: named list, species/pkg versions 
  moduleServer(id, function(input, output, session) {
    
    comparisons <- reactiveVal(list())
    
    # ── Auto-capture on each run ────────────────────────────────────────────
    # observeEvent on run_trigger(): whenever a new DE result is produced,
    # snapshot params + results + summary and append.
    observeEvent(run_trigger(), {
      res <- tryCatch(de_results_all(), error = function(e) NULL)
      if (is.null(res) || nrow(res) == 0) return()
      p <- tryCatch(contrast_params(), error = function(e) list())
      
      # abs() >= lfc is what the exported table was filtered by; with lfc = 0
      # (or absent) it is always TRUE, so an unfiltered run counts as before.
      n_up   <- sum(res$log2FoldChange > 0 & abs(res$log2FoldChange) >= (p$lfc %||% 0) & !is.na(res$padj) & res$padj <= (p$alpha %||% 0.05), na.rm = TRUE)
      n_down <- sum(res$log2FoldChange < 0 & abs(res$log2FoldChange) >= (p$lfc %||% 0) & !is.na(res$padj) & res$padj <= (p$alpha %||% 0.05), na.rm = TRUE)
      
      label <- if (identical(p$mode, "interaction")) {
        paste0("interaction: ", p$interaction %||% "?")
      } else {
        paste0(p$factor %||% "?", ": ", p$test %||% "?", " vs ", p$reference %||% "?")
      }
      
      entry <- list(
        # The modulo has to happen BEFORE as.integer(), not after. Milliseconds
        # since the epoch is ~1.79e12, far past .Machine$integer.max, so
        # as.integer() saw the raw value, returned NA and warned
        # ("NAs introduced by coercion to integer range") on EVERY DE run --
        # and NA %% anything is still NA, so every comparison was stamped id=NA
        # and serialised as null in the reproducibility manifest.
        id        = as.integer((as.numeric(Sys.time()) * 1000) %% .Machine$integer.max),
        label     = label,
        mode      = p$mode %||% "contrast",
        factor    = p$factor %||% NA_character_,
        test      = p$test %||% NA_character_,
        reference = p$reference %||% NA_character_,
        interaction = p$interaction %||% NA_character_,
        alpha     = p$alpha %||% 0.05,
        lfc       = p$lfc %||% 0,
        shrinkage = isTRUE(p$shrinkage),
        design    = p$design %||% NA_character_,
        # ── Result-defining settings repro_module.R already knows how to use ──
        # prefilter: repro_module.R:250 reads hist[[i]]$prefilter and emits it as
        # EXECUTABLE code ("keep <- rowSums(counts(dds) >= 10) >= 6"). Without it
        # the generated rerun_DE.R fell back to a comment, and the script could
        # not reproduce its own padj -- the pre-filter decides which genes exist,
        # which decides the baseMean quantile DESeq2's independent filtering
        # picks, which decides EVERY padj.
        # The other four decide de_sig_genes(), i.e. the literal gene vector
        # handed to GO and KEGG, so without them the bundle stored that vector
        # without recording which filters produced it.
        # contrast_params() (rnaseq_de_module.R:758-779) has supplied all five all
        # along; only this list was not passing them on. Kept unwrapped (no %||%)
        # because repro_module.R tests them with is.list()/isTRUE() and needs to
        # be able to tell "not recorded" (NULL) from "recorded as FALSE".
        prefilter        = p$prefilter,
        reference_levels = p$reference_levels,
        sig_only         = isTRUE(p$sig_only),
        up_only          = isTRUE(p$up_only),
        down_only        = isTRUE(p$down_only),
        env       = if (is.function(session_info)) tryCatch(session_info(), error = function(e) NULL) else NULL,
        n_total   = nrow(res),
        n_up      = n_up,
        n_down    = n_down,
        time      = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
        result    = res
      )
      comparisons(c(comparisons(), list(entry)))
      showNotification(paste0("\U0001f5c3\ufe0f Comparison logged: ", label,
                              " (", n_up, " up / ", n_down, " down)"),
                       type = "message", duration = 4)
    }, ignoreInit = TRUE)
    
    # ── Table ───────────────────────────────────────────────────────────────
    output$log_table <- DT::renderDT({
      h <- comparisons()
      if (length(h) == 0)
        return(DT::datatable(data.frame(Message = "No comparisons run yet."),
                             rownames = FALSE, options = list(dom = "t")))
      df <- data.frame(
        `#`        = seq_along(h),
        Comparison = vapply(h, `[[`, "", "label"),
        Up         = vapply(h, `[[`, 0L, "n_up"),
        Down       = vapply(h, `[[`, 0L, "n_down"),
        Total      = vapply(h, `[[`, 0L, "n_total"),
        Time       = vapply(h, `[[`, "", "time"),
        check.names = FALSE
      )
      DT::datatable(df, rownames = FALSE, selection = "single",
                    options = list(dom = "tp", pageLength = 10))
    })
    
    observeEvent(input$delete_sel, {
      sel <- input$log_table_rows_selected
      req(sel)
      h <- comparisons()
      comparisons(h[-sel])
    })
    observeEvent(input$clear_all, { comparisons(list()) })
    
    # ── Helpers for export ──────────────────────────────────────────────────
    .safe <- function(x) gsub("[^A-Za-z0-9]+", "_", x)
    
    # ── R script export ─────────────────────────────────────────────────────
    output$dl_r <- downloadHandler(
      filename = function() paste0("de_session_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".R"),
      content  = function(file) {
        h <- comparisons()
        env <- if (is.function(session_info)) tryCatch(session_info(), error = function(e) NULL) else NULL
        env_lines <- if (!is.null(env)) {
          c("## --- Session environment ---",
            vapply(names(env), function(nm) sprintf("## %s: %s", nm, env[[nm]]), ""),
            "## ---------------------------")
        } else character(0)
        
        lines <- c(
          "## DE session script (auto-generated by One for All)",
          paste0("## Generated: ", format(Sys.time())),
          env_lines,
          "## Assumes a fitted DESeqDataSet named `dds` is already in your workspace.",
          "library(DESeq2)", "library(ashr)", "",
          if (length(h) == 0) "## (no comparisons were logged)" else NULL
        )
        for (e in h) {
          obj <- paste0("res_", .safe(e$label))
          if (identical(e$mode, "interaction")) {
            call_line <- sprintf('%s <- results(dds, name = "%s", alpha = %s)',
                                 obj, e$interaction, format(e$alpha))
          } else {
            call_line <- sprintf('%s <- results(dds, contrast = c("%s","%s","%s"), alpha = %s)',
                                 obj, e$factor, e$test, e$reference, format(e$alpha))
          }
          
          # Re-create shrinkage ONLY if it was actually applied to this run —
          # otherwise the regenerated script would silently change the values.
          shrink_line <- if (isTRUE(e$shrinkage)) {
            if (identical(e$mode, "interaction")) {
              sprintf('%s <- lfcShrink(dds, coef = "%s", type = "ashr", res = %s)',
                      obj, e$interaction, obj)
            } else {
              sprintf('%s <- lfcShrink(dds, contrast = c("%s","%s","%s"), type = "ashr", res = %s)',
                      obj, e$factor, e$test, e$reference, obj)
            }
          } else NULL
          
          lines <- c(lines, "",
                     paste0("## ", e$label, "  |  ", e$time,
                            "  |  up=", e$n_up, " down=", e$n_down, " total=", e$n_total,
                            "  |  shrinkage=", isTRUE(e$shrinkage)),
                     paste0("## design: ", e$design %||% "NA"),
                     call_line,
                     shrink_line,
                     if (!is.null(e$lfc) && e$lfc > 0)
                       sprintf('%s <- subset(%s, !is.na(padj) & padj <= %s & abs(log2FoldChange) >= %s)',
                               obj, obj, format(e$alpha), format(e$lfc))
                     else
                       sprintf('%s <- subset(%s, !is.na(padj) & padj <= %s)', obj, obj, format(e$alpha)),
                     sprintf('write.csv(as.data.frame(%s), "%s.csv")', obj, obj)
          )
        }
        writeLines(lines, file)
      }
    )
    
    # ── JSON export ─────────────────────────────────────────────────────────
    output$dl_json <- downloadHandler(
      filename = function() paste0("de_session_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".json"),
      content  = function(file) {
        h <- comparisons()
        manifest <- lapply(h, function(e) {
          e2 <- e
          # keep the JSON light: replace the full table with the gene-level rows
          # so parameters + results are both present but structured.
          e2$result <- NULL
          e2$genes  <- {
            r <- e$result
            keep <- !is.na(r$padj) & r$padj <= (e$alpha %||% 0.05)
            if (!is.null(e$lfc) && e$lfc > 0) keep <- keep & abs(r$log2FoldChange) >= e$lfc
            r2 <- r[keep, , drop = FALSE]
            lapply(seq_len(nrow(r2)), function(i) list(
              gene_id        = r2$gene_id[i],
              log2FoldChange = round(r2$log2FoldChange[i], 4),
              padj           = signif(r2$padj[i], 4)
            ))
          }
          e2
        })
        jsonlite::write_json(list(
          generated = format(Sys.time()),
          session_info = if (is.function(session_info)) tryCatch(session_info(), error = function(e) NULL) else NULL,
          n_comparisons = length(h),
          comparisons = manifest
        ), file, auto_unbox = TRUE, pretty = TRUE, null = "null")
      }
    )
    
    # ── ZIP of full CSV tables ──────────────────────────────────────────────
    output$dl_zip <- downloadHandler(
      filename = function() paste0("de_session_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".zip"),
      content  = function(file) {
        h <- comparisons()
        tmp <- file.path(tempdir(), paste0("de_zip_", as.integer(Sys.time())))
        dir.create(tmp, showWarnings = FALSE)
        # session_info.txt — annotation/package versions + full sessionInfo(),
        # so this export stands on its own without the full analysis bundle.
        si_path <- file.path(tmp, "session_info.txt")
        con <- file(si_path, open = "wt"); on.exit(close(con), add = TRUE)
        writeLines("=== One fo All DE session: environment ===", con)
        writeLines(paste0("Generated: ", format(Sys.time())), con)
        writeLines("", con)
        env <- if (is.function(session_info)) tryCatch(session_info(), error = function(e) NULL) else NULL
        if (!is.null(env)) {
          writeLines("--- Annotation / package versions ---", con)
          for (nm in names(env)) writeLines(sprintf("%s: %s", nm, env[[nm]]), con)
          writeLines("", con)
        }
        writeLines("--- sessionInfo() ---", con)
        writeLines(utils::capture.output(utils::sessionInfo()), con)
        paths <- character(0)
        for (i in seq_along(h)) {
          e  <- h[[i]]
          fn <- file.path(tmp, sprintf("%02d_%s.csv", i, .safe(e$label)))
          write.csv(e$result, fn, row.names = FALSE)
          paths <- c(paths, fn, si_path)
        }
        # also drop the JSON manifest inside the zip
        mf <- file.path(tmp, "manifest.json")
        jsonlite::write_json(list(
          generated = format(Sys.time()),
          comparisons = lapply(h, function(e) list(
            label = e$label, factor = e$factor, test = e$test,
            reference = e$reference, interaction = e$interaction,
            alpha = e$alpha, lfc = e$lfc, design = e$design,
            n_up = e$n_up, n_down = e$n_down, n_total = e$n_total, time = e$time))
        ), mf, auto_unbox = TRUE, pretty = TRUE)
        paths <- c(paths, mf)
        old <- setwd(tmp); on.exit(setwd(old), add = TRUE)
        utils::zip(zipfile = file, files = basename(paths))
      }
    )
    
    return(list(comparisons = comparisons))
  })
}