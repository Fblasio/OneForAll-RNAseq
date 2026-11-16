# ── in-memory cache ────────────────────────────────────────
.enrich_cache <- new.env(hash = TRUE, parent = emptyenv())

.cache_key <- function(...) digest::digest(list(...))

clear_enrichment_cache <- function() {
  rm(list = ls(envir = .enrich_cache), envir = .enrich_cache)
  message("Enrichment cache cleared.")
}



# ── GO ─────────────────────────────────────────────────────
run_go_enrichment <- function(gene_ids,
                              process        = "BP",
                              padj_method    = "BH",
                              pval_cutoff    = 0.05,
                              qval_cutoff    = 0.05,
                              simplify_cutoff = 0.7,
                              use_cache      = TRUE,
                              universe       = NULL,
                              context = "default") {

  key <- .cache_key(sort(gene_ids), process, padj_method,
                    pval_cutoff, qval_cutoff, simplify_cutoff,
                    sort(universe), context, ACTIVE_SPECIES)

  if (use_cache && exists(key, envir = .enrich_cache)) {
    message("GO: returning cached result")
    return(get(key, envir = .enrich_cache))
  }

  geneids_upper <- canon_ids(gene_ids)
  
  run_params <- list(
    gene_ids        = geneids_upper,
    universe        = if (!is.null(universe)) canon_ids(universe) else NULL,
    process         = process,
    padj_method     = padj_method,
    pval_cutoff     = pval_cutoff,
    qval_cutoff     = qval_cutoff,
    simplify_cutoff = simplify_cutoff,
    org_db_pkg      = ORG_DB_PKG,
    key_type        = KEY_TYPE,
    species         = ACTIVE_SPECIES
  )
  
  # message("=== GO DEBUG ===")
  # message("ACTIVE_SPECIES: ", ACTIVE_SPECIES)
  # message("KEY_TYPE: ", KEY_TYPE)
  # message("ORG_DB_PKG: ", ORG_DB_PKG)
  # message("N input genes: ", length(gene_ids))
  # message("N normalized genes: ", length(geneids_upper))
  # message("Head input genes: ", paste(utils::head(gene_ids, 10), collapse = ", "))
  # message("Head normalized genes: ", paste(utils::head(geneids_upper, 10), collapse = ", "))
  # 
  # if (!is.null(universe)) {
  #   universe_norm <- normalize_gene_ids(universe)
  #   message("Universe length raw: ", length(universe))
  #   message("Universe length normalized: ", length(universe_norm))
  #   message("Head universe raw: ", paste(utils::head(universe, 10), collapse = ", "))
  #   message("Head universe normalized: ", paste(utils::head(universe_norm, 10), collapse = ", "))
  # } else {
  #   message("Universe: NULL")
  # }

  go_res <- enrichGO(
    gene     = geneids_upper,
    OrgDb    = ORG_DB,
    keyType  = KEY_TYPE,
    ont           = process,
    pAdjustMethod = padj_method,
    pvalueCutoff  = pval_cutoff,
    qvalueCutoff  = qval_cutoff,
    universe = if (!is.null(universe)) {
      canon_ids(universe)
    } else NULL,
    readable      = FALSE,
    minGSSize     = 3,
    maxGSSize     = 2000,
    pool          = ifelse(process == "ALL", TRUE, FALSE)
  )

  if (is.null(go_res) || nrow(go_res@result) == 0) {
    result <- empty_result("No GO terms found below the chosen p/q-value cutoffs.")
    result$params <- run_params
    if (use_cache) assign(key, result, envir = .enrich_cache)
    return(result)
  }

  ego_final <- if (process != "ALL") {
    tryCatch(
      clusterProfiler::simplify(go_res, cutoff = simplify_cutoff,
                                by = "p.adjust", select_fun = min),
      error = function(e) go_res
    )
  } else go_res

  ego_final@result$geneID <- vapply(
    ego_final@result$geneID,
    function(x) {
      ids   <- canon_ids(trimws(strsplit(x, "/")[[1]]))
      valid <- ids %in% geneids_upper
      paste(ids[valid], collapse = "/") 
    },
    character(1)
  )
  ego_final@result$Count <- vapply(
    ego_final@result$geneID,
    function(x) if (nzchar(x)) length(strsplit(x, "/")[[1]]) else 0L,
    integer(1)
  )
  ego_final@result <- ego_final@result[ego_final@result$Count > 0, ]

  ego_final <- tryCatch({
    gr <- as.character(ego_final@result$GeneRatio)
    br <- as.character(ego_final@result$BgRatio)
    fe <- DOSE::parse_ratio(gr) / DOSE::parse_ratio(br)
    fe[!is.finite(fe)] <- NA
    ego_final@result$FoldEnrichment <- fe
    ego_final
  }, error = function(e) {
    message("FoldEnrichment computation failed: ", e$message)
    ego_final
  })

  result <- list(table = as.data.frame(ego_final@result), object = ego_final, params = run_params)
  if (use_cache) assign(key, result, envir = .enrich_cache)
  result
}

# ── KEGG ───────────────────────────────────────────────────
run_kegg_enrichment <- function(gene_ids,
                                padj_method = "BH",
                                pval_cutoff = 0.05,
                                qval_cutoff = 0.05,
                                use_cache   = TRUE,
                                universe    = NULL) {
  # `universe` joins the cache key: the same gene list against a different
  # background is a different result, and omitting it here would serve a
  # stale genome-background table after the background changes.
  key <- .cache_key(sort(gene_ids), "kegg", padj_method, pval_cutoff, qval_cutoff,
                    sort(universe), ACTIVE_SPECIES)

  if (use_cache && exists(key, envir = .enrich_cache)) {
    message("KEGG: returning cached result")
    return(get(key, envir = .enrich_cache))
  }
  
  # Arabidopsis KEGG uses TAIR IDs directly; other species need Entrez/NCBI IDs.
  if (ACTIVE_SPECIES == "arabidopsis") {
    kegg_input <- normalize_gene_ids(gene_ids)
  } else {
    ids_canon <- canon_ids(gene_ids)
    entrez_df <- convert_to_entrez(ids_canon)
    if (is.null(entrez_df) || nrow(entrez_df) == 0) {
      result <- empty_result(paste0(
        "KEGG needs Entrez/NCBI IDs, but none of the ", length(unique(ids_canon)),
        " input gene ID(s) could be mapped for this species. KEGG enrichment skipped."))
      if (use_cache) assign(key, result, envir = .enrich_cache)
      return(result)
    }
    kegg_input <- entrez_df$ENTREZID
  }

  # Background must live in the SAME ID space as `kegg_input`, so it goes
  # through the identical conversion. NULL -> enrichKEGG falls back to every
  # gene in the TERM2GENE table (whole-genome background).
  kegg_universe <- NULL
  if (!is.null(universe) && length(universe) > 0) {
    if (ACTIVE_SPECIES == "arabidopsis") {
      kegg_universe <- unique(normalize_gene_ids(universe))
    } else {
      uni_df <- convert_to_entrez(canon_ids(universe))
      if (!is.null(uni_df) && nrow(uni_df) > 0) {
        kegg_universe <- unique(uni_df$ENTREZID)
      }
    }
    if (is.null(kegg_universe) || length(kegg_universe) == 0) {
      kegg_universe <- NULL
      message("KEGG: no background gene could be mapped to KEGG ID space \u2014 ",
              "falling back to the whole-genome background.")
    } else {
      # Foreground must be a subset of the background, or enrichKEGG silently
      # drops the difference and GeneRatio stops matching the input list.
      kegg_universe <- union(kegg_universe, kegg_input)
      message(sprintf("KEGG: using a custom background of %d gene(s).",
                      length(kegg_universe)))
    }
  }

  kegg_res <- enrichKEGG(
    gene     = kegg_input,
    organism = KEGG_CODE,
    universe = kegg_universe,
    minGSSize         = 4,
    maxGSSize         = 1000,
    pAdjustMethod     = padj_method,
    pvalueCutoff      = pval_cutoff,
    qvalueCutoff      = qval_cutoff,
    use_internal_data = FALSE
  )

  if (is.null(kegg_res) || nrow(kegg_res@result) == 0) {
    result <- empty_result(sprintf(
      "No KEGG pathways found below the chosen cutoffs (%d gene ID(s) used as KEGG input).",
      length(unique(kegg_input))))
    if (use_cache) assign(key, result, envir = .enrich_cache)
    return(result)
  }

  result <- list(table = as.data.frame(kegg_res@result), object = kegg_res)
  if (use_cache) assign(key, result, envir = .enrich_cache)
  result
}

# ── Custom dotplot ─────────────────────────────────────────
# pinned: character vector of Description values that must always appear in the plot.
# Remaining slots (top_n - n_pinned) are filled with top-ranked free terms by GeneRatio.
make_enrich_dotplot <- function(result_df, title = "", top_n = 20, pinned = character(0), font_size = 13) {
  req_cols <- c("Description", "GeneRatio", "p.adjust", "Count")
  missing  <- setdiff(req_cols, colnames(result_df))
  if (length(missing) > 0)
    stop("Missing columns: ", paste(missing, collapse = ", "))

  df <- result_df
  df$Description  <- as.character(df$Description)
  df$GeneRatio    <- as.character(df$GeneRatio)
  df$p.adjust     <- as.numeric(df$p.adjust)
  df$Count        <- as.numeric(df$Count)

  df$GeneRatioNum <- sapply(df$GeneRatio, function(x) {
    parts <- strsplit(x, "/")[[1]]
    if (length(parts) == 2) as.numeric(parts[1]) / as.numeric(parts[2])
    else NA_real_
  })

  df <- df[!is.na(df$GeneRatioNum) & !is.na(df$p.adjust), ]
  df <- df[order(df$GeneRatioNum, decreasing = TRUE), ]
  # Apply pinned logic: pinned terms always shown; remaining slots filled with top free terms
  if (length(pinned) > 0) {
    is_pinned <- df$Description %in% pinned
    df_pinned <- df[is_pinned,  ]
    df_free   <- df[!is_pinned, ]
    n_free    <- max(0L, top_n - nrow(df_pinned))
    df        <- rbind(df_free[seq_len(min(n_free, nrow(df_free))), ], df_pinned)
    df        <- df[order(df$GeneRatioNum, decreasing = TRUE), ]
  } else {
    df <- head(df, top_n)
  }
  df$Description <- factor(df$Description, levels = rev(df$Description))

  ggplot(df, aes(x = GeneRatioNum, y = Description,
                 size = Count, color = p.adjust)) +
    geom_point() +
    scale_color_gradient(low = "#C0392B", high = "#00008b", name = "p.adjust") +
    scale_size_continuous(name = "Gene count", range = c(3, 10)) +
    xlab("Gene Ratio") + ylab(NULL) +
    ggtitle(title) +
    theme_bw(base_size = font_size) +
    theme(axis.text.y  = element_text(size = font_size),
          axis.text.x  = element_text(size = font_size),
          axis.title.x = element_text(size = font_size + 2, face = "bold"),
          plot.title   = element_text(size = font_size + 4, face = "bold"),
          legend.text  = element_text(size = font_size - 1),
          legend.title = element_text(size = font_size,     face = "bold"))
}

# ── Flexible, publication-oriented enrichment plot ──────────────────────────
# make_enrichment_plot(): superset of make_enrich_dotplot()/.make_dotplot()
# (both kept as-is elsewhere for backward compatibility), with plot-type,
# x-axis metric, and color choices exposed as arguments instead of being
# hard-coded.
#
# Written against the {Description, GeneRatio, p.adjust, Count} columns that
# both run_go_enrichment() and run_kegg_enrichment() already produce, plus
# the optional FoldEnrichment column run_go_enrichment() adds. This makes it
# usable from GO today and from KEGG later without any changes here — see
# the FoldEnrichment fallback below.
#
# style        : "dot" (classic dotplot), "lollipop" (point + stem from 0),
#                or "bar" (horizontal bar chart).
# x_metric     : "GeneRatio" (k/n, always available) or "FoldEnrichment"
#                ("Enrichment Factor" in the UI; observed/expected ratio).
#                Falls back to GeneRatio with a message() if the column is
#                missing or entirely NA (e.g. a result table that predates
#                this feature, or a source that never computed it).
# color_mode   : "gradient" (color/fill mapped to a significance metric) or
#                "flat" (single user-chosen color, no color legend).
# color_metric : "p.adjust" or "neglog10_padj" (-log10(p.adjust); larger =
#                more significant, sometimes preferred so "more intense
#                color" reads as "more significant").
# size_by_count: TRUE maps point/lollipop-head size to Count; FALSE uses a
#                fixed size (ignored for style = "bar", whose length already
#                encodes the x metric).
# label_wrap   : character width to wrap long term Descriptions onto
#                multiple lines (0/NULL disables wrapping).
# theme_style  : "bw", "minimal", or "classic" — swaps the ggplot2 base theme.
make_enrichment_plot <- function(result_df,
                                  title         = "",
                                  top_n         = 20,
                                  pinned        = character(0),
                                  font_size     = 13,
                                  style         = c("dot", "lollipop", "bar"),
                                  x_metric      = c("GeneRatio", "FoldEnrichment"),
                                  color_mode    = c("gradient", "flat"),
                                  color_metric  = c("p.adjust", "neglog10_padj"),
                                  color_low     = "#3B7DD8",
                                  color_high    = "#C0392B",
                                  flat_color    = "#3B7DD8",
                                  size_by_count = TRUE,
                                  point_size    = 4,
                                  label_wrap    = 45,
                                  theme_style   = c("bw", "minimal", "classic"),
                                  x_lab         = NULL,
                                  sort_desc     = TRUE) {

  style        <- match.arg(style)
  x_metric     <- match.arg(x_metric)
  color_mode   <- match.arg(color_mode)
  color_metric <- match.arg(color_metric)
  theme_style  <- match.arg(theme_style)

  req_cols <- c("Description", "GeneRatio", "p.adjust", "Count")
  missing  <- setdiff(req_cols, colnames(result_df))
  if (length(missing) > 0)
    stop("Missing columns: ", paste(missing, collapse = ", "))

  df <- result_df
  df$Description <- as.character(df$Description)
  df$GeneRatio    <- as.character(df$GeneRatio)
  df$p.adjust     <- as.numeric(df$p.adjust)
  df$Count        <- as.numeric(df$Count)

  df$GeneRatioNum <- sapply(df$GeneRatio, function(x) {
    parts <- strsplit(x, "/")[[1]]
    if (length(parts) == 2) as.numeric(parts[1]) / as.numeric(parts[2]) else NA_real_
  })

  # FoldEnrichment ("Enrichment Factor"): use it only if present and usable;
  # otherwise transparently fall back to GeneRatio so callers never error
  # out just because a source table doesn't compute it (e.g. future KEGG use).
  has_fold <- "FoldEnrichment" %in% colnames(df) &&
    !all(is.na(suppressWarnings(as.numeric(df$FoldEnrichment))))
  if (identical(x_metric, "FoldEnrichment") && !has_fold) {
    message("make_enrichment_plot: FoldEnrichment not available in this ",
            "result table; falling back to GeneRatio.")
    x_metric <- "GeneRatio"
  }
  df$XMetric <- if (identical(x_metric, "FoldEnrichment")) {
    as.numeric(df$FoldEnrichment)
  } else {
    df$GeneRatioNum
  }
  x_axis_label <- x_lab %||% if (identical(x_metric, "FoldEnrichment")) {
    "Enrichment Factor (Fold Enrichment)"
  } else {
    "Gene Ratio"
  }

  df$ColorMetric <- if (identical(color_metric, "neglog10_padj")) {
    -log10(df$p.adjust)
  } else {
    df$p.adjust
  }
  color_legend_name <- if (identical(color_metric, "neglog10_padj")) {
    "-log10(p.adjust)"
  } else {
    "p.adjust"
  }

  df <- df[!is.na(df$XMetric) & !is.na(df$p.adjust), ]
  if (nrow(df) == 0) {
    return(
      ggplot() +
        annotate("text", x = 0.5, y = 0.5,
                 label = "No terms available to plot with the current settings.",
                 size = 6, color = "gray40") +
        theme_void()
    )
  }
  df <- df[order(df$XMetric, decreasing = sort_desc), ]

  # Pinned logic: pinned terms always shown; remaining slots filled with the
  # top-ranked free terms by the CHOSEN x metric (not hard-coded to GeneRatio).
  if (length(pinned) > 0) {
    is_pinned <- df$Description %in% pinned
    df_pinned <- df[is_pinned,  ]
    df_free   <- df[!is_pinned, ]
    n_free    <- max(0L, top_n - nrow(df_pinned))
    df        <- rbind(df_free[seq_len(min(n_free, nrow(df_free))), ], df_pinned)
    df        <- df[order(df$XMetric, decreasing = sort_desc), ]
  } else {
    df <- head(df, top_n)
  }

  if (!is.null(label_wrap) && label_wrap > 0) {
    df$Description <- vapply(df$Description,
                              function(s) paste(strwrap(s, width = label_wrap), collapse = "\n"),
                              character(1))
  }
  # Order retained top-to-bottom as sorted above (rev() so ggplot's bottom-up
  # discrete axis reads top-to-bottom the way the data frame is ordered).
  df$Description <- factor(df$Description, levels = rev(unique(df$Description)))

  p <- ggplot(df, aes(x = XMetric, y = Description))

  if (identical(style, "lollipop")) {
    p <- p + geom_segment(aes(x = 0, xend = XMetric, y = Description, yend = Description),
                           color = "grey60", linewidth = 0.6)
  }

  if (identical(style, "bar")) {
    if (identical(color_mode, "gradient")) {
      p <- p + geom_col(aes(fill = ColorMetric), width = 0.7) +
        scale_fill_gradient(low = color_low, high = color_high, name = color_legend_name)
    } else {
      p <- p + geom_col(fill = flat_color, width = 0.7)
    }
  } else {
    # Four explicit, unambiguous branches instead of trying to merge partial
    # aes()/NULL-size objects — safer than it looks, and easy to verify.
    if (identical(color_mode, "gradient") && isTRUE(size_by_count)) {
      p <- p + geom_point(aes(size = Count, color = ColorMetric)) +
        scale_color_gradient(low = color_low, high = color_high, name = color_legend_name) +
        scale_size_continuous(name = "Gene count", range = c(3, 10))
    } else if (identical(color_mode, "gradient") && !isTRUE(size_by_count)) {
      p <- p + geom_point(aes(color = ColorMetric), size = point_size) +
        scale_color_gradient(low = color_low, high = color_high, name = color_legend_name)
    } else if (!identical(color_mode, "gradient") && isTRUE(size_by_count)) {
      p <- p + geom_point(aes(size = Count), color = flat_color) +
        scale_size_continuous(name = "Gene count", range = c(3, 10))
    } else {
      p <- p + geom_point(color = flat_color, size = point_size)
    }
  }

  p <- p +
    xlab(x_axis_label) + ylab(NULL) +
    ggtitle(title)

  p <- p + switch(theme_style,
                   bw      = theme_bw(base_size = font_size),
                   minimal = theme_minimal(base_size = font_size),
                   classic = theme_classic(base_size = font_size))

  p + theme(axis.text.y  = element_text(size = font_size),
            axis.text.x  = element_text(size = font_size),
            axis.title.x = element_text(size = font_size + 2, face = "bold"),
            plot.title   = element_text(size = font_size + 4, face = "bold"),
            legend.text  = element_text(size = font_size - 1),
            legend.title = element_text(size = font_size,     face = "bold"))
}
