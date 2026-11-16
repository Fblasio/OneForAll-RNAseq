# =============================================================================
# repro_module.R
# -----------------------------------------------------------------------------
# Reproducibility layer: one download button producing a .zip that describes
# and regenerates everything done in a session.
#
# The bundle holds EVERY run of the session, not the most recent one of each
# kind. N heatmaps produce N matrix/script pairs; N GO, KEGG, GSEA and DE runs
# each produce N blocks in their respective scripts. Refused and failed runs are
# recorded too, so a bundle reporting N runs came from N attempts.
#
# Contents
#   manifest.json      every parameter of every run, in scope
#   session_info.txt   sessionInfo() plus annotation package versions
#   gene_lists/*.txt   the exact gene vectors used, one file per run where a run
#                      has its own (GO/KEGG queries and backgrounds, GSEA
#                      rankings), plus the session-level selections
#   rerun_DE.R         every DE contrast, grouped by design formula
#   rerun_GO.R         every GO run, with its table and its figure
#   rerun_KEGG.R       every KEGG run, with its table and its figure
#   rerun_GSEA.R       every GSEA run (tables only; see the note below)
#   rerun_heatmap*.R   one script per generated heatmap
#   heatmap*matrix.csv the exact plotted matrix for each, so the figure
#                      regenerates rather than merely being described
#   README.txt         what the bundle is, and its stated limits
#
# File naming: with a single heatmap the pair is heatmap_matrix.csv /
# rerun_heatmap.R; with several it is heatmap_NN_matrix.csv /
# rerun_heatmap_NN.R, numbered in generation order.
#
# What the bundle does NOT reproduce, stated here so the omissions are visible:
#   * rerun_heatmap*.R draws the heat body and the colour scale only. Annotation
#     bars, gene-symbol labels and the cluster split are recorded in
#     manifest.json but not drawn.
#   * rerun_GSEA.R regenerates tables, not the NES plots. The GSEA panel drives a
#     different plotting function from the GO/KEGG one and does not share the
#     style-snapshot mechanism.
#   * rerun_DE.R needs the user's own counts matrix and metadata; raw counts are
#     not bundled. The count pre-filter that was applied IS emitted.
#
# Every collector is optional. One that is not supplied simply produces no
# corresponding output, and the generated scripts say so rather than failing.
#
# Wiring:
#   UI:     reproDownloadUI("repro")            # a single download button
#   Server: reproServer("repro", collectors = list(...))
#           where each collector is a function() returning a plain list/vector.
# =============================================================================

reproDownloadUI <- function(id, label = "\u2b07 Download analysis bundle (.zip)") {
  ns <- NS(id)
  downloadButton(ns("dl_bundle"), label, class = "btn-sm btn-info", width = "100%")
}

# collectors: named list of zero-arg functions, all optional. Recognised names:
#
#   heatmap_config  list of per-heatmap parameter lists, one per generated figure
#   heatmap_matrix  list of list(matrix = <numeric, genes x conditions>,
#                                config = <the same frozen config>), one per
#                   figure. Produces the CSV + rerun script pair(s).
#   de_params       the CURRENT DE contrast. Fallback only, used when
#                   de_history is absent.
#   de_history      EVERY DE contrast of the session (the comparison log's
#                   `comparisons` reactive). Preferred over de_params.
#   go_params       GO params for the current run, split into heatmap/de panels
#   go_history      EVERY GO run's params; each element carries $gene_ids and
#                   optionally $universe and $plot
#   kegg_params     EVERY KEGG run's params (same shape as go_history)
#   gsea_params     one run_gsea() params list, or a list of them. Fallback only.
#   gsea_history    EVERY GSEA run, each list(params=, ranked_ids=, ranked_scores=)
#   gene_lists      named list of character vectors, written one file per name
#   annotation_info named character: org.db, KEGG code, ITAG version, ID-map md5
#
# The *_history collectors are what make the bundle cover a whole session rather
# than its last action, and are preferred wherever both forms exist.
#
# The style snapshot: a GO or KEGG history element may carry $plot, the frozen
# appearance of that run's figure (see enrich_plot_snapshot() in
# utils/enrichment_functions.R). When present the emitted script redraws the
# figure; when absent it emits a note and writes the table only, so histories
# recorded before the snapshot existed still produce a valid script.
reproServer <- function(id, collectors = list()) {
  moduleServer(id, function(input, output, session) {

    safe_call <- function(nm) {
      f <- collectors[[nm]]
      if (is.null(f)) return(NULL)
      tryCatch(f(), error = function(e) NULL)
    }

    # Normalize gsea_params into a named list-of-runs, whether the caller gave
    # one run (a flat params list with a $database field) or several (a list of
    # such lists). Returns list(key = params, ...) so each run gets its own block.
    normalize_gsea <- function(gp) {
      if (is.null(gp)) return(list())
      # single run: a flat list that itself has a 'database' element
      if (!is.null(gp$database)) return(list(gsea_1 = gp))
      # already a list of runs
      out <- list()
      for (i in seq_along(gp)) {
        nm <- names(gp)[i]
        if (is.null(nm) || !nzchar(nm)) nm <- paste0("gsea_", i)
        out[[nm]] <- gp[[i]]
      }
      out
    }

    # Replace bulky per-run fields with their size.
    #
    # The history collectors carry the actual data as well as the parameters:
    # go/kegg entries hold `gene_ids` + `universe`, gsea entries hold the whole
    # ranked vector, and de entries hold the complete result data frame. All of
    # that is already written to gene_lists/*.txt and to the rerun scripts, so
    # repeating it in manifest.json would duplicate megabytes and make the file
    # unreadable -- a single DE result table would dwarf everything else.
    # Sizes are kept (n_gene_ids, n_universe, n_result, ...) so the manifest
    # still states the shape of every run.
    #
    # Only the five names below are stripped, which is why the $plot style
    # snapshot survives into the manifest.
    .slim <- function(h) {
      if (is.null(h) || length(h) == 0) return(h)
      bulky <- c("gene_ids", "universe", "ranked_ids", "ranked_scores", "result")
      lapply(h, function(e) {
        if (!is.list(e)) return(e)
        for (nm in intersect(bulky, names(e))) {
          e[[paste0("n_", nm)]] <-
            if (is.data.frame(e[[nm]])) nrow(e[[nm]]) else length(e[[nm]])
          e[[nm]] <- NULL
        }
        # gsea_history nests its parameters one level down.
        if (is.list(e$params)) e$params <- .slim(list(e$params))[[1]]
        e
      })
    }

    # The key set below is an allow-list. A collector not named here is absent
    # from manifest.json by design because it is serialised elsewhere:
    # `heatmap_matrix` becomes the heatmap CSV(s), `gene_lists` becomes
    # gene_lists/*.txt, and only their sizes appear here.
    #
    # The three *_history keys must stay in this list. Without them the manifest
    # would describe GO, KEGG, GSEA and DE by a single run each while the rerun
    # scripts reproduced all of them -- two accounts of the same session that
    # disagree on how much of it happened.
    build_manifest <- function() {
      list(
        app            = "One for All",
        generated_at   = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
        active_species = if (exists("ACTIVE_SPECIES")) ACTIVE_SPECIES else NA,
        heatmap_config = safe_call("heatmap_config"),
        de_params      = safe_call("de_params"),
        de_history     = .slim(safe_call("de_history")),
        go_params      = safe_call("go_params"),
        go_history     = .slim(safe_call("go_history")),
        kegg_params    = .slim(safe_call("kegg_params")),
        gsea_params    = normalize_gsea(safe_call("gsea_params")),
        gsea_history   = .slim(safe_call("gsea_history")),
        annotation     = as.list(safe_call("annotation_info")),
        gene_list_sizes = {
          gl <- safe_call("gene_lists")
          if (is.null(gl)) NULL else lapply(gl, length)
        }
      )
    }

    write_session_info <- function(path) {
      con <- file(path, open = "wt")
      on.exit(close(con))
      writeLines("=== One for All reproducibility: session & annotation versions ===", con)
      writeLines(paste0("Generated: ", format(Sys.time())), con)
      writeLines("", con)
      ann <- safe_call("annotation_info")
      if (!is.null(ann)) {
        writeLines("--- Annotation ---", con)
        for (nm in names(ann)) writeLines(sprintf("%s: %s", nm, ann[[nm]]), con)
        writeLines("", con)
      }
      writeLines("--- sessionInfo() ---", con)
      writeLines(utils::capture.output(utils::sessionInfo()), con)
    }

    # Generate a runnable script reproducing EVERY DE contrast this session.
    #
    # Source of truth is the `de_history` collector (the DE comparison log's
    # `comparisons` reactive) — the only place the app accumulates every run.
    # Falls back to the single-run `de_params` collector if history isn't wired,
    # so this degrades gracefully rather than emitting nothing.
    #
    # Contrasts are GROUPED BY DESIGN FORMULA: a session can rebuild the dds
    # with a different design, and every contrast under one design shares a
    # single DESeq() fit — which is both faithful to what the app did and much
    # faster to re-run than refitting per contrast.
    write_rerun_script <- function(path) {
      hist <- safe_call("de_history")
      if (is.null(hist) || length(hist) == 0) {
        de <- safe_call("de_params")
        # is.null(de) is NOT the right test. The de_params collector is built as
        # list(design = input[[...]], factor = input[[...]], ...), and a list of
        # NULLs is a list of length 6, not NULL -- so before any contrast has been
        # run this fallback produced one "record" whose every field was empty, the
        # length-0 guard below never fired, and the bundle shipped a runnable-
        # looking script with design "~ condition" and
        # results(dds, contrast = c("", "", "")). That script errors if run and
        # describes an analysis that never happened, while manifest.json correctly
        # reported no DE history -- two files in one bundle disagreeing.
        # A usable fallback needs exactly the three values results(contrast = )
        # takes; anything less is not a contrast and is dropped.
        have <- unlist(de[c("factor", "test", "reference")])
        hist <- if (length(have) == 3 && all(nzchar(as.character(have)))) list(de) else list()
      }

      con <- file(path, open = "wt")
      on.exit(close(con))
      writeLines(c(
        "# Auto-generated by One for All - reproduces EVERY differential expression",
        "# contrast run in this session.",
        "# You must provide the same counts matrix + metadata used in the app.",
        "# This script intentionally does NOT hide any parameter.",
        "",
        "library(DESeq2)",
        "",
        "# counts: genes x samples integer matrix; coldata: data.frame with sample metadata",
        "# counts  <- as.matrix(read.csv('counts.csv', row.names = 1))",
        "# coldata <- read.csv('metadata.csv', row.names = 1)",
        ""), con)

      if (length(hist) == 0) {
        writeLines("# No DE contrast was run in this session.", con)
        return(invisible(NULL))
      }

      designs <- vapply(hist, function(e) as.character(e$design %||% "~ condition"), "")
      idx_all <- seq_along(hist)

      for (d in unique(designs)) {
        idx <- idx_all[designs == d]

        # Emit the recorded pre-filter as executable code rather than as an
        # instruction to remember it. min_count decides which genes exist, which
        # decides the baseMean quantile DESeq2's independent filtering picks,
        # which decides every padj in the table -- so a script that merely asked
        # the reader to "apply the same pre-filtering" would not reproduce the
        # result it claims to.
        #
        # rnaseq_dds_module.R stores the filter on the object as
        # S4Vectors::metadata(dds)$cheski_prefilter and rnaseq_de_module.R
        # carries it through contrast_params$prefilter. Read per design group,
        # because min_group_size derives from the first design variable.
        pf <- hist[[idx[1]]]$prefilter
        filter_lines <- if (is.list(pf) && !is.null(pf$min_count) &&
                            !is.null(pf$min_group_size)) {
          c(sprintf(paste0("# pre-filter recorded from the app run: at least %s count(s) ",
                           "in >= %s samples"),
                    format(pf$min_count), format(pf$min_group_size)),
            sprintf("# (the smallest group of '%s'). In the app this kept %s of %s genes.",
                    pf$grouping_var %||% "?",
                    format(pf$n_genes_after %||% NA), format(pf$n_genes_before %||% NA)),
            sprintf("keep <- rowSums(counts(dds) >= %s) >= %s",
                    format(pf$min_count), format(pf$min_group_size)),
            "dds  <- dds[keep, ]")
        } else {
          "# NOTE: no pre-filter was recorded for this design; apply the same one you used in the app."
        }

        writeLines(c(
          sprintf("# ===== design: %s =====", d),
          sprintf("design_formula <- %s", deparse(d)),
          "dds <- DESeqDataSetFromMatrix(round(counts), coldata, as.formula(design_formula))",
          filter_lines,
          "dds <- DESeq(dds)",
          ""), con)

        for (i in idx) {
          e     <- hist[[i]]
          alpha <- e$alpha %||% 0.05
          is_int <- identical(e$mode, "interaction")
          # the comparison log stores the interaction coef under `interaction`;
          # the single-run de_params fallback may use `interaction_name`.
          coef_nm <- e$interaction %||% e$interaction_name %||% ""

          res_line <- if (is_int) {
            sprintf("res <- results(dds, name = %s, alpha = %s)",
                    deparse(coef_nm), alpha)
          } else {
            sprintf("res <- results(dds, contrast = c(%s, %s, %s), alpha = %s)",
                    deparse(e$factor %||% ""), deparse(e$test %||% ""),
                    deparse(e$reference %||% ""), alpha)
          }

          # Reproduce the ashr shrinkage step when the app applied it —
          # matches rnaseq_de_module.R's lfcShrink() call exactly.
          shrink_line <- if (isTRUE(e$shrinkage)) {
            if (is_int) {
              sprintf("res <- lfcShrink(dds, coef = %s, type = 'ashr', res = res)",
                      deparse(coef_nm))
            } else {
              sprintf("res <- lfcShrink(dds, contrast = c(%s, %s, %s), type = 'ashr', res = res)",
                      deparse(e$factor %||% ""), deparse(e$test %||% ""),
                      deparse(e$reference %||% ""))
            }
          } else {
            "# (no LFC shrinkage was applied to this contrast in the app)"
          }

          lfc_line <- if (!is.null(e$lfc) && is.finite(suppressWarnings(as.numeric(e$lfc))) &&
                          as.numeric(e$lfc) > 0) {
            sprintf("# app also applied an |log2FC| >= %s filter when calling genes significant",
                    e$lfc)
          } else NULL

          writeLines(c(
            sprintf("# --- comparison %02d: %s ---", i, e$label %||% "contrast"),
            res_line,
            shrink_line,
            "res <- as.data.frame(res)",
            "res$gene_id <- rownames(res)",
            "res <- res[!is.na(res$log2FoldChange), ]",
            "res <- res[order(res$padj), ]",
            lfc_line,
            sprintf("write.csv(res, 'DE_results_comparison_%02d.csv', row.names = FALSE)", i),
            ""), con)
        }
      }
    }

    # ── The figure, not just the table ───────────────────────────────────────
    # rerun_GO.R / rerun_KEGG.R emit the FIGURE as well as the table. Every
    # style control is recorded per run (see attach_plot_snapshot() in
    # utils/enrichment_functions.R), so the plot a table was shown as is
    # reproducible rather than merely described.
    #
    # make_enrichment_plot() is app code, not a CRAN package, so a standalone
    # script cannot library() it. It is therefore INLINED into the head of each
    # rerun script that needs it, rather than shipped as a separate file the
    # scripts source().
    #
    # Inlining duplicates ~115 lines between rerun_GO.R and rerun_KEGG.R, which
    # would be the wrong trade in hand-maintained source -- but these are
    # GENERATED artefacts that nobody edits, and the duplication buys two things
    # that matter more here:
    #   * each script is self-contained. Copy rerun_GO.R out of the bundle, mail
    #     it to a reviewer, and it still draws its figures.
    #   * no silent failure mode. The previous `if (file.exists(...)) source()`
    #     guard skipped the load when the helper was missing and the script then
    #     died several lines later on "could not find function".
    #
    # srcref is stripped before deparse() so the emitted source comes from the
    # AST rather than the original file text: with keep.source = TRUE (the
    # RStudio default) deparse() would otherwise reproduce the comments, and
    # their UTF-8 would have to survive writeLines() into a bundle that is
    # opened on an unknown locale.
    .fn_source <- function(f) {
      f2 <- f
      attr(f2, "srcref") <- NULL
      deparse(f2, width.cutoff = 500L)
    }

    # The plotting preamble, or character(0) when it is not needed (no run in
    # this history recorded a style) or not available.
    plot_helper_lines <- function(hist) {
      if (!exists("make_enrichment_plot", mode = "function")) return(character(0))
      if (is.null(hist) || length(hist) == 0) return(character(0))
      if (!any(vapply(hist, function(p) is.list(p) && is.list(p$plot), logical(1))))
        return(character(0))
      src <- .fn_source(make_enrichment_plot)
      c(
        "# ---------------------------------------------------------------------",
        "# The app's own plotting function, inlined so this script is standalone.",
        "# Copied verbatim from utils/enrichment_functions.R; do not edit.",
        "# ---------------------------------------------------------------------",
        "library(ggplot2)",
        "",
        "`%||%` <- function(a, b) {",
        "  if (is.null(a) || length(a) == 0) return(b)",
        "  if (is.character(a) && length(a) == 1L && !nzchar(a)) return(b)",
        "  a",
        "}",
        "",
        paste("make_enrichment_plot <-", src[1]),
        src[-1],
        "# --------------------------- end of inlined helper -------------------",
        ""
      )
    }

    # The emitted make_enrichment_plot() call for one recorded run. `p$plot` is
    # the snapshot frozen when Run was pressed; a run recorded before this
    # feature existed simply has no $plot and gets no figure block, so old
    # bundles and old sessions degrade to the previous behaviour.
    emit_plot_block <- function(p, df_expr, out_png) {
      s <- p$plot
      if (!is.list(s)) {
        return("# (no plot style was recorded for this run - table only)")
      }
      c(
        sprintf("plot_df <- %s", df_expr),
        sprintf(paste0(
          "p <- make_enrichment_plot(plot_df,\n",
          "        title = %s, top_n = %s, pinned = %s, font_size = %s,\n",
          "        style = %s, x_metric = %s,\n",
          "        color_mode = %s, color_metric = %s,\n",
          "        color_low = %s, color_high = %s, flat_color = %s,\n",
          "        size_by_count = %s, theme_style = %s, label_wrap = %s,\n",
          "        sort_desc = %s)"),
          deparse(s$title %||% ""), s$top_n %||% 20,
          paste(deparse(as.character(s$pinned %||% character(0))), collapse = ""),
          s$font_size %||% 13,
          deparse(s$style %||% "dot"), deparse(s$x_metric %||% "GeneRatio"),
          deparse(s$color_mode %||% "gradient"), deparse(s$color_metric %||% "p.adjust"),
          deparse(s$color_low %||% "#3B7DD8"), deparse(s$color_high %||% "#C0392B"),
          deparse(s$flat_color %||% "#3B7DD8"),
          if (isTRUE(s$size_by_count)) "TRUE" else "FALSE",
          deparse(s$theme_style %||% "bw"), s$label_wrap %||% 45,
          if (isTRUE(s$sort_desc)) "TRUE" else "FALSE"),
        sprintf("ggsave('%s', p, width = 10, height = 8, dpi = 150)", out_png)
      )
    }

    # Generate a runnable script for the GO enrichment step(s) actually run.
    # gene_dir: where per-run gene files are written (gene_lists/ in the bundle).
    write_rerun_GO <- function(path, gene_dir) {
      hist <- safe_call("go_history")
      con <- file(path, open = "wt")
      on.exit(close(con))
      writeLines(c(
        "# Auto-generated by One for All - reproduces EVERY GO enrichment run this session.",
        "# Requires the gene_lists/*.txt files included in this same bundle.",
        "# Each run also redraws its figure with the plot style that was on screen",
        "# when Run GO was pressed. The plotting function is inlined below, so this",
        "# script needs nothing from the bundle except gene_lists/.",
        "",
        "library(clusterProfiler)",
        "",
        plot_helper_lines(hist)
      ), con)

      # p: one run's params list (must carry $gene_ids and optionally $universe).
      emit_block <- function(p, idx) {
        if (is.null(p)) return(invisible(NULL))
        # A failed run is recorded so the attempt is visible, but emitting
        # enrichGO() for it would produce a script that errors on purpose.
        # Note it and move on -- the count of blocks then matches the count of
        # SUCCESSFUL runs, while the comment preserves the audit trail.
        if (isTRUE(p$failed)) {
          writeLines(c(
            sprintf("# --- run %d: %s --- FAILED, not reproducible ---",
                    idx, p$context %||% "GO"),
            sprintf("#     error: %s", p$error %||% "unknown"),
            ""), con)
          return(invisible(NULL))
        }
        label     <- p$context %||% sprintf("GO_run_%d", idx)
        gene_file <- sprintf("go_run_%02d_genes", idx)
        uni_file  <- sprintf("go_run_%02d_universe", idx)

        # write this run's exact gene vector (and universe, if any) to disk
        genes <- p$gene_ids %||% character(0)
        if (length(genes) > 0)
          writeLines(as.character(genes), file.path(gene_dir, paste0(gene_file, ".txt")))
        has_uni <- !is.null(p$universe) && length(p$universe) > 0
        if (has_uni)
          writeLines(as.character(p$universe), file.path(gene_dir, paste0(uni_file, ".txt")))

        writeLines(c(
          sprintf("# --- run %d: %s ---", idx, label),
          sprintf("library(%s)", p$org_db_pkg %||% "org.At.tair.db"),
          sprintf("gene_ids <- readLines('gene_lists/%s.txt')", gene_file),
          if (has_uni) sprintf("universe <- readLines('gene_lists/%s.txt')", uni_file)
          else         "universe <- NULL",
          sprintf(
            # minGSSize/maxGSSize come from the recorded run, not from literals.
            # Hardcoding them here instead would let the script test a
            # different set of terms from the app, adjust p over a different
            # number of them, and produce a table that does not match the
            # figure it claims to reproduce.
            "ego <- enrichGO(gene = gene_ids, OrgDb = %s, keyType = %s, ont = %s,\n                pAdjustMethod = %s, pvalueCutoff = %s, qvalueCutoff = %s,\n                universe = universe, minGSSize = %s, maxGSSize = %s)",
            p$org_db_pkg %||% "org.At.tair.db", deparse(p$key_type %||% "TAIR"),
            deparse(p$ontology %||% p$process %||% "BP"), deparse(p$padj_method %||% "BH"),
            p$pval_cutoff %||% 0.05, p$qval_cutoff %||% 0.05,
            p$min_gs_size %||% 3, p$max_gs_size %||% 2000
          ),
          sprintf(
            "ego <- tryCatch(clusterProfiler::simplify(ego, cutoff = %s, by = 'p.adjust', select_fun = min), error = function(e) ego)",
            p$simplify_cutoff %||% 0.7
          ),
          sprintf("write.csv(as.data.frame(ego@result), 'GO_results_run_%02d.csv')", idx),
          emit_plot_block(p, "as.data.frame(ego)",
                          sprintf("GO_plot_run_%02d.png", idx)),
          ""
        ), con)
      }

      if (is.null(hist) || length(hist) == 0) {
        writeLines("# No GO enrichment was run in this session.", con)
        return(invisible(NULL))
      }
      for (i in seq_along(hist)) emit_block(hist[[i]], i)
    }

    # Generate a runnable script for the KEGG step(s) actually run.
    #
    # Source of truth is the `kegg_params` collector, which carries every run
    # from both KEGG panels (heatmap selection and DEG list).
    #
    # ID SPACE -- the one thing that makes this script different from rerun_GO.R:
    # run_kegg_enrichment() records `gene_ids` in the app's own namespace (TAIR,
    # or canonical Solyc), NOT the IDs handed to enrichKEGG. For Arabidopsis
    # they are the same. For every other species KEGG keys on Entrez and the app
    # converts through the side map named in `id_map_source` (whose md5 is in
    # manifest.json -> annotation), so the conversion is emitted explicitly
    # rather than pretending the recorded IDs are what enrichKEGG saw.
    write_rerun_KEGG <- function(path, gene_dir) {
      hist <- safe_call("kegg_params")
      con <- file(path, open = "wt")
      on.exit(close(con))
      writeLines(c(
        "# Auto-generated by One for All - reproduces EVERY KEGG enrichment run this session.",
        "# Requires the gene_lists/*.txt files included in this same bundle.",
        "# Each run also redraws its figure with the plot style that was on screen",
        "# when Run KEGG was pressed. The plotting function is inlined below, so this",
        "# script needs nothing from the bundle except gene_lists/.",
        "",
        "library(clusterProfiler)",
        "",
        plot_helper_lines(hist)
      ), con)

      if (is.null(hist) || length(hist) == 0) {
        writeLines("# No KEGG enrichment was run in this session.", con)
        return(invisible(NULL))
      }

      for (i in seq_along(hist)) {
        p <- hist[[i]]
        if (is.null(p)) next

        # Same convention as rerun_GO.R: a failed run is documented, not emitted
        # as code that is known to error.
        if (isTRUE(p$failed)) {
          writeLines(c(
            sprintf("# --- run %d --- FAILED, not reproducible ---", i),
            sprintf("#     error: %s", p$error %||% "unknown"),
            ""), con)
          next
        }

        gene_file <- sprintf("kegg_run_%02d_genes", i)
        uni_file  <- sprintf("kegg_run_%02d_universe", i)
        genes <- p$gene_ids %||% character(0)
        if (length(genes) > 0)
          writeLines(as.character(genes), file.path(gene_dir, paste0(gene_file, ".txt")))
        has_uni <- !is.null(p$universe) && length(p$universe) > 0
        if (has_uni)
          writeLines(as.character(p$universe), file.path(gene_dir, paste0(uni_file, ".txt")))

        id_src <- p$id_map_source %||% NA_character_
        map_lines <- if (identical(p$species, "arabidopsis")) {
          "# Arabidopsis: KEGG 'ath' keys on TAIR IDs, so no conversion is needed."
        } else if (!is.na(id_src)) {
          c(sprintf("# %s: KEGG '%s' keys on Entrez/NCBI IDs. The app mapped the IDs",
                    p$species %||% "this species", p$kegg_code %||% "?"),
            sprintf("# below through %s (md5 in manifest.json -> annotation).", id_src),
            sprintf("id_map   <- read.csv('%s', colClasses = 'character')", id_src),
            "gene_ids <- unique(id_map$ENTREZID[id_map[[1]] %in% gene_ids])",
            "gene_ids <- gene_ids[nzchar(gene_ids)]",
            if (has_uni) c(
              "universe <- unique(id_map$ENTREZID[id_map[[1]] %in% universe])",
              "universe <- universe[nzchar(universe)]") else NULL)
        } else {
          c(sprintf("# %s: KEGG '%s' keys on Entrez/NCBI IDs, but no ID map was",
                    p$species %||% "this species", p$kegg_code %||% "?"),
            "# recorded for this run -- convert gene_ids to Entrez before running.")
        }

        writeLines(c(
          sprintf("# --- run %d ---", i),
          sprintf("# %s of %s query gene(s) had a KEGG-compatible ID in the app.",
                  format(p$n_genes_mapped %||% p$n_genes_input %||% NA),
                  format(p$n_genes_input %||% NA)),
          sprintf("gene_ids <- readLines('gene_lists/%s.txt')", gene_file),
          if (has_uni) sprintf("universe <- readLines('gene_lists/%s.txt')", uni_file)
          else         "universe <- NULL",
          map_lines,
          sprintf(
            # Bounds come from the recorded run, not from literals, for the same
            # reason rerun_GO.R stopped hardcoding minGSSize.
            "res <- enrichKEGG(gene = gene_ids, organism = %s, universe = universe,\n                  minGSSize = %s, maxGSSize = %s, pAdjustMethod = %s,\n                  pvalueCutoff = %s, qvalueCutoff = %s, use_internal_data = FALSE)",
            deparse(p$kegg_code %||% "ath"), p$min_gs_size %||% 4, p$max_gs_size %||% 1000,
            deparse(p$padj_method %||% "BH"), p$pval_cutoff %||% 0.05, p$qval_cutoff %||% 0.05),
          # as.data.frame(), not res@result: the slot holds every tested pathway
          # and only the S3 method applies the cutoffs. Same reason the app now
          # converts through the accessor -- see utils/enrichment_functions.R.
          sprintf("write.csv(as.data.frame(res), 'KEGG_results_run_%02d.csv')", i),
          emit_plot_block(p, "as.data.frame(res)",
                          sprintf("KEGG_plot_run_%02d.png", i)),
          ""
        ), con)
      }
    }

    # Generate a runnable script for the GSEA step(s) actually run.
    #
    # Mirrors run_gsea() in gsea_functions.R: it rebuilds the named ranked
    # vector from the two gene_lists files, then calls gseGO()/gseKEGG() with
    # the same eps = 0 and seed = TRUE the app uses, so the permutation-based
    # p-values reproduce exactly rather than merely resembling the app's output.
    #
    # Tables only: the GSEA panel plots through make_gsea_summary_plot(), which
    # takes a different argument set from the GO/KEGG plotter and is not covered
    # by the style-snapshot mechanism, so no figure is emitted here.
    #
    # gene_dir: where per-run ranked-vector files are written (gene_lists/).
    write_rerun_GSEA <- function(path, gsea_runs, gene_dir) {
      con <- file(path, open = "wt")
      on.exit(close(con))
      writeLines(c(
        "# Auto-generated by One for All - reproduces EVERY GSEA (fgsea) run this session.",
        "# Requires the gene_lists/gsea_run_*_ranked_ids.txt and *_ranked_scores.txt",
        "# files included in this same bundle.",
        "",
        "library(clusterProfiler)",
        ""
      ), con)

      if (length(gsea_runs) == 0) {
        writeLines("# No GSEA was run in this session.", con)
        return(invisible(NULL))
      }

      for (i in seq_along(gsea_runs)) {
        run <- gsea_runs[[i]]
        if (is.null(run)) next
        p   <- run$params %||% run          # tolerate params-only entries too
        db  <- p$database %||% "GO"
        key <- sprintf("gsea_run_%02d", i)

        # write this run's exact ranked vector to its own pair of files
        ids    <- run$ranked_ids    %||% character(0)
        scores <- run$ranked_scores %||% numeric(0)
        if (length(ids) > 0) {
          writeLines(as.character(ids),    file.path(gene_dir, paste0(key, "_ranked_ids.txt")))
          writeLines(as.character(scores), file.path(gene_dir, paste0(key, "_ranked_scores.txt")))
        }

        common <- c(
          sprintf("# --- run %d : %s GSEA ---", i, db),
          sprintf("library(%s)", p$org_db_pkg %||% "org.At.tair.db"),
          sprintf("ranked_ids    <- readLines('gene_lists/%s_ranked_ids.txt')", key),
          sprintf("ranked_scores <- as.numeric(readLines('gene_lists/%s_ranked_scores.txt'))", key),
          "rank_vec <- setNames(ranked_scores, ranked_ids)",
          "rank_vec <- sort(rank_vec, decreasing = TRUE)  # gseGO/gseKEGG require this",
          ""
        )

        if (identical(db, "GO")) {
          body <- sprintf(
            "res <- gseGO(geneList = rank_vec, OrgDb = %s, keyType = %s, ont = %s,\n             minGSSize = %s, maxGSSize = %s, pvalueCutoff = %s,\n             pAdjustMethod = %s, eps = 0, seed = TRUE, verbose = FALSE)",
            p$org_db_pkg %||% "org.At.tair.db", deparse(p$key_type %||% "TAIR"),
            deparse(p$ontology %||% "BP"), p$min_size %||% 10, p$max_size %||% 500,
            p$pval_cutoff %||% 0.05, deparse(p$padj_method %||% "BH"))
        } else {
          body <- sprintf(
            "# NOTE: for KEGG, ranked_ids were already remapped to KEGG's ID space in-app.\nres <- gseKEGG(geneList = rank_vec, organism = %s,\n               minGSSize = %s, maxGSSize = %s, pvalueCutoff = %s,\n               pAdjustMethod = %s, eps = 0, seed = TRUE, verbose = FALSE,\n               use_internal_data = FALSE)",
            deparse(p$kegg_code %||% "ath"), p$min_size %||% 10, p$max_size %||% 500,
            p$pval_cutoff %||% 0.05, deparse(p$padj_method %||% "BH"))
        }

        writeLines(c(common, body,
          sprintf("write.csv(as.data.frame(res@result), 'GSEA_results_run_%02d.csv')", i),
          ""), con)
      }
    }

    # Write the exact plotted matrix plus a script that redraws each heatmap.
    #
    # The `heatmap_matrix` collector may hand over either a single
    # list(matrix, config) or a list of such pairs, one per heatmap generated in
    # the session. Both forms are accepted: one pair produces the unnumbered
    # heatmap_matrix.csv / rerun_heatmap.R, and N pairs produce N numbered sets.
    #
    # Returns the NUMBER of matrices written, not a logical, because the README
    # branch below needs the count to describe the layout it produced.
    write_heatmap_reproduction <- function(dir) {
      hm <- safe_call("heatmap_matrix")
      if (is.null(hm)) return(0L)
      runs <- if (!is.null(hm$matrix)) list(hm) else hm
      n <- 0L
      for (i in seq_along(runs)) {
        if (!is.null(runs[[i]]$matrix)) {
          n <- n + 1L
          # Single run keeps the historical file names so nothing that already
          # reads a bundle has to learn a new layout.
          csv_name <- if (length(runs) == 1L) "heatmap_matrix.csv"
                      else sprintf("heatmap_%02d_matrix.csv", i)
          scr_name <- if (length(runs) == 1L) "rerun_heatmap.R"
                      else sprintf("rerun_heatmap_%02d.R", i)
          .write_one_heatmap(runs[[i]], file.path(dir, csv_name),
                             file.path(dir, scr_name), csv_name)
        }
      }
      n
    }

    # csv_name is passed separately from csv_path because the generated script
    # reads the CSV by BASENAME (both files land side by side in the bundle),
    # while the writer needs the full path.
    .write_one_heatmap <- function(hm, csv_path, script_path, csv_name) {
      m <- hm$matrix
      cfg <- hm$config %||% list()

      # matrix -> CSV (genes as first column so it round-trips with row.names)
      utils::write.csv(as.data.frame(m), csv_path, row.names = TRUE)

      cols <- cfg$colors %||% list(low = "#3B7DD8", mid = "white", high = "#C0392B")
      zlim <- cfg$zlim %||% NA
      lo   <- cols$low  %||% "#3B7DD8"
      mid  <- cols$mid  %||% "white"
      hi   <- cols$high %||% "#C0392B"
      lnk  <- cfg$linkage %||% "complete"
      # Same branch the app takes when building col_fun: "none"/"log2" are
      # all-positive magnitudes with no meaningful zero-centre and get a robust
      # 2%/98% quantile scale; every other method is centred on 0 and keeps the
      # symmetric zlim from the slider. Emitting one fixed ramp for both, as
      # this script used to, mis-coloured every "none"/"log2" figure.
      quantile_scale <- identical(cfg$normalization, "none") ||
                        identical(cfg$normalization, "log2")
      # Output image named after its own CSV, so several heatmaps in one bundle
      # cannot overwrite each other's PNG when the scripts are run in place.
      png_name <- paste0(sub(".csv", "", csv_name, fixed = TRUE), "_reproduced.png")
      con <- file(script_path, open = "wt"); on.exit(close(con))
      writeLines(c(
        "# Auto-generated by One for All - redraws the heatmap from the plotted matrix",
        sprintf("# saved alongside this script (%s).", csv_name),
        "#",
        "# REPRODUCED EXACTLY: the numeric matrix, and the colour scale including",
        "#   the app's method-dependent choice of anchor points.",
        "#",
        "# NOT REPRODUCED EXACTLY: the row order.",
        "#   The app clusters rows on the SAMPLE-level VST matrix (Pearson",
        "#   correlation distance, linkage below) -- not on the group-mean display",
        "#   matrix saved here. That substrate is not part of this bundle, so the",
        "#   dendrogram below is built from the display matrix instead: same",
        "#   distance and linkage, different observations, so expect a similar but",
        "#   not identical row order. For ward.D2 the app additionally uses",
        "#   sqrt(2*(1-r)), the Euclidean embedding of correlation distance, which",
        "#   ComplexHeatmap does not apply.",
        "#",
        "#   Earlier versions of this script used ComplexHeatmap's DEFAULT clustering",
        "#   (Euclidean, on the display matrix) and row_km = k (k-means). Neither is",
        "#   what the app does -- it cuts the hierarchical tree with cutree() -- so",
        "#   both have been removed rather than left to look authoritative.",
        "",
        "library(ComplexHeatmap)",
        "library(circlize)",
        "",
        sprintf("mat <- as.matrix(read.csv('%s', row.names = 1, check.names = FALSE))", csv_name),
        "",
        if (quantile_scale) {
          c("zmin <- as.numeric(quantile(mat, 0.02, na.rm = TRUE))",
            "zmax <- as.numeric(quantile(mat, 0.98, na.rm = TRUE))",
            "if (!is.finite(zmin) || !is.finite(zmax) || zmin == zmax) {",
            "  zmin <- min(mat, na.rm = TRUE); zmax <- max(mat, na.rm = TRUE)",
            "}",
            sprintf("col_fun <- colorRamp2(c(zmin, (zmin + zmax)/2, zmax), c('%s','%s','%s'))",
                    lo, mid, hi))
        } else {
          # The midpoint is NOT hardcoded to 0. It is 0 only when 0 lies inside
          # the recorded window; otherwise it is the centre of that window.
          # This mirrors heatmap_color_breaks() (utils/color_utils.R) line for
          # line -- writing `c(zmin, 0, zmax)` here would hand colorRamp2()
          # non-increasing breaks for a window like 0.5..2 and redraw a figure
          # the app never produced.
          c(if (length(zlim) == 2 && all(is.finite(zlim))) {
              sprintf("zmin <- %s; zmax <- %s", format(zlim[1]), format(zlim[2]))
            } else {
              "zmax <- max(abs(mat), na.rm = TRUE); zmin <- -zmax"
            },
            "if (!is.finite(zmin) || !is.finite(zmax)) { zmin <- -1; zmax <- 1 }",
            "if (zmin > zmax) { .t <- zmin; zmin <- zmax; zmax <- .t }",
            "if (zmin == zmax) { zmin <- -1; zmax <- 1 }",
            "zmid <- if (zmin < 0 && zmax > 0) 0 else (zmin + zmax) / 2",
            sprintf("col_fun <- colorRamp2(c(zmin, zmid, zmax), c('%s','%s','%s'))",
                    lo, mid, hi))
        },
        "",
        sprintf("ht <- Heatmap(mat, name = 'z',\n              col = col_fun,\n              cluster_rows = %s,\n              clustering_distance_rows = 'pearson',\n              clustering_method_rows = '%s',\n              cluster_columns = %s,\n              clustering_distance_columns = 'spearman',\n              clustering_method_columns = 'complete',\n              show_row_names = nrow(mat) <= 100)",
                if (isTRUE(cfg$cluster_rows)) "TRUE" else "FALSE",
                lnk,
                if (isTRUE(cfg$cluster_cols)) "TRUE" else "FALSE"),
        "",
        sprintf("png('%s', width = 1200, height = 1400, res = 150)", png_name),
        "draw(ht)",
        "dev.off()",
        sprintf("cat('Wrote %s\\n')", png_name)
      ), con)
      TRUE
    }

    output$dl_bundle <- downloadHandler(
      filename = function()
        paste0("OneForAll_analysis_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".zip"),
      content = function(file) {
        tmp <- file.path(tempdir(), paste0("cheski_bundle_", as.integer(Sys.time())))
        dir.create(tmp, showWarnings = FALSE)
        dir.create(file.path(tmp, "gene_lists"), showWarnings = FALSE)
        on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

        # manifest.json
        manifest <- build_manifest()
        writeLines(
          jsonlite::toJSON(manifest, pretty = TRUE, auto_unbox = TRUE, null = "null"),
          file.path(tmp, "manifest.json"))

        # session_info.txt
        write_session_info(file.path(tmp, "session_info.txt"))

        # session_log.txt (if utils/log_utils.R is sourced and events were logged)
        if (exists("get_session_log", mode = "function")) {
          log_lines <- get_session_log(session)
          if (length(log_lines) > 0) {
            writeLines(log_lines, file.path(tmp, "session_log.txt"))
          }
        }

        # gene lists
        gl <- safe_call("gene_lists")
        if (!is.null(gl)) {
          for (nm in names(gl)) {
            v <- gl[[nm]]
            if (!is.null(v) && length(v) > 0)
              writeLines(as.character(v),
                         file.path(tmp, "gene_lists", paste0(nm, ".txt")))
          }
        }

        # rerun scripts
        gene_dir <- file.path(tmp, "gene_lists")
        write_rerun_script(file.path(tmp, "rerun_DE.R"))
        # The plotting helper is inlined into each of these two by
        # plot_helper_lines(), so there is no separate file to write and no
        # ordering dependency between them.
        write_rerun_GO(file.path(tmp, "rerun_GO.R"), gene_dir)
        write_rerun_KEGG(file.path(tmp, "rerun_KEGG.R"), gene_dir)

        # Prefer the full session history; fall back to the single current run
        # (normalize_gsea) if no history collector is wired.
        gsea_runs <- safe_call("gsea_history")
        if (is.null(gsea_runs) || length(gsea_runs) == 0)
          gsea_runs <- normalize_gsea(safe_call("gsea_params"))
        write_rerun_GSEA(file.path(tmp, "rerun_GSEA.R"), gsea_runs, gene_dir)

        n_heatmaps    <- write_heatmap_reproduction(tmp)
        wrote_heatmap <- n_heatmaps > 0L

        # README
        readme <- c(
          "One for All analysis bundle",
          "===========================",
          "manifest.json    - all parameters (heatmap config + DE + GO + KEGG + GSEA) in scope.",
          "session_info.txt - R/package/annotation versions for reproducibility.",
          "session_log.txt  - timeline of uploads/warnings/errors this session (if any).",
          "gene_lists/      - exact gene vectors used (selected, DEGs, cart, DE-significant,",
          "                   GO input/universe vectors, and GSEA ranked ids+scores).",
          "rerun_DE.R       - runnable script to reproduce the DE step with your data,",
          "                   including the count pre-filter recorded from each run.",
          "rerun_GO.R       - runnable script to reproduce the GO enrichment step(s).",
          "rerun_KEGG.R     - runnable script to reproduce the KEGG enrichment step(s).",
          "                   For species whose KEGG code keys on Entrez it also emits the",
          "                   ID-map conversion the app performed (see annotation.id_map_*).",
          "rerun_GSEA.R     - runnable script to reproduce the GSEA step(s) from the ranked",
          "                   vectors in gene_lists/ (same eps/seed as the app).",
          "",
          "SCOPE: the GO and KEGG runs recorded here are those made from the heatmap",
          "enrichment panel and from the DEG (RNA-seq) panel. GO runs made from the",
          "standalone 'GO Analysis' feature are deliberately NOT recorded, so an empty",
          "or short go_history means 'nothing was run in those two panels', not",
          "'no GO was ever run'. Each entry names the panel it came from in its",
          "`context` field, and the runs are listed in the order they happened.",
          "",
          "rerun_GO.R and rerun_KEGG.R each also redraw every run's FIGURE",
          "(GO_plot_run_NN.png / KEGG_plot_run_NN.png) with the plot type, axis, colours,",
          "theme and label wrapping that were on screen when THAT run was launched - not",
          "whatever the style panel said at download time. The plotting function is inlined",
          "at the top of each script, so either file can be lifted out of this bundle and",
          "run on its own (it still needs gene_lists/).")
        if (wrote_heatmap) {
          readme <- c(readme,
          if (n_heatmaps == 1L) {
            c("rerun_heatmap.R  - runnable script that redraws the heatmap from heatmap_matrix.csv.",
              "heatmap_matrix.csv - the exact numeric matrix that was plotted.")
          } else {
            c(sprintf("rerun_heatmap_NN.R / heatmap_NN_matrix.csv - one PAIR per heatmap generated"),
              sprintf("                   this session (%d in total, numbered in the order you", n_heatmaps),
              "                   generated them). Each script redraws its own matrix with the",
              "                   settings that were active when THAT heatmap was built.")
          },
          "                   Reproduces the plotted matrix and the colour scale exactly.",
          "                   Row ORDER is close but not identical: the app clusters on the",
          "                   sample-level VST matrix, which is not included in this bundle,",
          "                   so the script clusters the display matrix with the same distance",
          "                   and linkage instead. The script's own header explains this.")
        } else {
          readme <- c(readme,
          "",
          "Note: no heatmap matrix was captured this session; the heatmap appearance is",
          "still described as configuration in manifest.json.")
        }
        writeLines(readme, file.path(tmp, "README.txt"))

        # zip it
        old <- setwd(tmp); on.exit(setwd(old), add = TRUE)
        zip::zip(zipfile = file, files = list.files(tmp, recursive = TRUE))
      }
    )
  })
}
