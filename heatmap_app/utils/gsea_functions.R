# =============================================================================
# gsea_functions.R  (goes in utils/, sourced from global.R after
#                    enrichment_functions.R / species_config.R / trust_utils.R)
# -----------------------------------------------------------------------------
# Rank-based enrichment (GSEA) on the FULL differential-expression result, as a
# complement to the over-representation (ORA) analysis in enrichment_functions.R.
#
# Why GSEA and not just ORA:
#   ORA (enrichGO/enrichKEGG) needs a hard DEG cut-off and throws away the rest
#   of the genome. GSEA uses the whole ranked list, so coordinated but
#   individually-subthreshold shifts in a pathway are still detected.
#
# Engine: clusterProfiler::gseGO() / gseKEGG(), both of which run fgsea
# (multilevel) internally. Using them (instead of fgsea directly) means we reuse
# the exact ORG_DB / KEY_TYPE / KEGG_CODE machinery already used by
# run_go_enrichment()/run_kegg_enrichment(), so GO gene sets, ID handling and
# species behaviour stay identical across ORA and GSEA.
#
# Returns the same {table, object, params} shape the GO/KEGG modules expect,
# and uses empty_result() (trust_utils.R) so a blank panel always says WHY.
#
# Depends (already loaded via global.R): clusterProfiler, DOSE. For the classic
# running-score plot the module uses enrichplot::gseaplot2() (an enrichplot
# function; enrichplot is a clusterProfiler dependency, so it is already
# installed — no new package required).
# =============================================================================


# -----------------------------------------------------------------------------
# build_gsea_rank(): turn a DE result data.frame into a named, decreasing,
# de-duplicated numeric vector suitable for gseGO()/gseKEGG().
#
# de_df   : a data.frame like the one returned by rnaseqDeServer()$de_results_all
#           — must contain gene_id and log2FoldChange; pvalue is used by the
#           default metric and stat by the "stat" metric when present.
# metric  :
#   "signed_p" (DEFAULT) sign(log2FC) * -log10(pvalue)
#                        Always available — works even when LFC shrinkage is ON
#                        and DESeq2's `stat` column is therefore absent.
#   "stat"               DESeq2 Wald statistic. Only usable when shrinkage is
#                        OFF (ashr drops `stat`); falls back to "signed_p" with
#                        a message if the column is missing.
#   "log2fc"             raw (or shrunken) log2 fold change.
#
# Returns NULL if nothing usable can be built (caller shows empty_result()).
# -----------------------------------------------------------------------------
build_gsea_rank <- function(de_df, metric = c("stat", "log2fc", "signed_p")) {
  metric <- match.arg(metric)

  if (is.null(de_df) || !is.data.frame(de_df) || nrow(de_df) == 0) return(NULL)
  if (!all(c("gene_id", "log2FoldChange") %in% colnames(de_df)))    return(NULL)

  lfc <- suppressWarnings(as.numeric(de_df$log2FoldChange))

  # Resolve which metric is ACTUALLY used into `resolved_metric`, then compute
  # `score` from a single, explicit if/else-if/else chain keyed on that
  # resolved value. No branch here relies on falling through from a block
  # above — each case sets `score` in its own arm, so adding/reordering a
  # future case can't accidentally leave `score` unset or stale.
  resolved_metric <- metric
  if (identical(metric, "stat") &&
      !("stat" %in% colnames(de_df) && any(is.finite(de_df$stat)))) {
    message("build_gsea_rank: 'stat' not available (LFC shrinkage on?) — ",
            "falling back to signed -log10(p).")
    resolved_metric <- "signed_p"
  }
  if (identical(resolved_metric, "signed_p") && !"pvalue" %in% colnames(de_df)) {
    message("build_gsea_rank: no 'pvalue' column — ranking by log2FC instead.")
    resolved_metric <- "log2fc"
  }

  if (identical(resolved_metric, "stat")) {
    score <- suppressWarnings(as.numeric(de_df$stat))

  } else if (identical(resolved_metric, "signed_p")) {
    p <- suppressWarnings(as.numeric(de_df$pvalue))
    # Floor ONLY exact/underflowed zero p-values (real, very-significant
    # genes) so -log10(p) stays finite. NA p-values are DESeq2's own
    # outlier flag (Cook's distance filtering) — those genes are
    # deliberately left as NA here so the `keep` step below drops them
    # from the ranking entirely, instead of turning them into a fake
    # top-of-list score.
    zero_mask <- is.finite(p) & p <= 0
    if (any(zero_mask)) {
      pos     <- p[is.finite(p) & p > 0]
      floor_p <- if (length(pos)) min(pos) else .Machine$double.xmin
      p[zero_mask] <- floor_p
    }
    score <- sign(lfc) * -log10(p)   # NA p -> NA score -> dropped by `keep`

  } else {   # resolved_metric == "log2fc"
    score <- lfc
  }

  ids  <- canon_ids(de_df$gene_id)
  keep <- is.finite(score) & !is.na(ids) & nzchar(ids)
  score <- score[keep]; ids <- ids[keep]
  if (length(score) < 2) return(NULL)

  # Canonicalization can map several raw IDs to one canonical ID (e.g. tomato
  # version suffixes). fgsea needs unique names — keep the largest-magnitude
  # score per gene, which is the standard collapse for ranked lists.
  if (any(duplicated(ids))) {
    ord   <- order(abs(score), decreasing = TRUE)
    ids   <- ids[ord]; score <- score[ord]
    unq   <- !duplicated(ids)
    ids   <- ids[unq]; score <- score[unq]
  }

  names(score) <- ids
  sort(score, decreasing = TRUE)      # gseGO/gseKEGG require decreasing order
}


# -----------------------------------------------------------------------------
# .remap_rank_for_kegg(): KEGG needs IDs in its own space. Mirrors the ID logic
# in run_kegg_enrichment(): Arabidopsis uses TAIR-style IDs directly; other
# species are converted to Entrez.
#
# Returns list(rank, n_query, n_mapped), or NULL when nothing usable survives.
# It returns the counts rather than just the vector because the caller cannot
# recover them afterwards: length(rank) counts ENTREZ RECORDS, not genes, and
# those differ in both directions (one gene can carry several Entrez IDs, and
# several genes can share one). Only this function still holds the map that
# relates the two, so the gene-level count has to be taken here.
# -----------------------------------------------------------------------------
.remap_rank_for_kegg <- function(rank_vec) {
  n_query <- length(unique(names(rank_vec)))
  if (identical(ACTIVE_SPECIES, "arabidopsis")) {
    v <- rank_vec
    names(v) <- normalize_gene_ids(names(v))
    n_mapped <- length(unique(names(v)))          # no conversion: nothing is lost
  } else {
    map <- convert_to_entrez(names(rank_vec))     # cols: <KEY_TYPE>, ENTREZID
    if (is.null(map) || nrow(map) == 0) return(NULL)
    # Column 1 is the KEY_TYPE column in both branches of convert_to_entrez(),
    # so this counts GENES that mapped. Taken before the filter below, which is
    # a no-op under the current data contract (load_id_map_impl() drops blank
    # ENTREZIDs at load, build_gsea_rank() drops non-finite scores).
    n_mapped <- length(unique(map[[1]]))
    idx <- match(map[[1]], names(rank_vec))
    v   <- rank_vec[idx]
    names(v) <- map$ENTREZID
  }
  v <- v[!is.na(names(v)) & nzchar(names(v)) & is.finite(v)]
  if (length(v) < 2) return(NULL)
  if (any(duplicated(names(v)))) {                # many-to-one -> max |score|
    ord <- order(abs(v), decreasing = TRUE)
    v   <- v[ord]; v <- v[!duplicated(names(v))]
  }
  list(rank = sort(v, decreasing = TRUE), n_query = n_query, n_mapped = n_mapped)
}


# -----------------------------------------------------------------------------
# run_gsea(): the GSEA counterpart to run_go_enrichment()/run_kegg_enrichment().
#
# rank_vec    : named decreasing numeric vector from build_gsea_rank().
# database    : "GO" (gseGO) or "KEGG" (gseKEGG).
# ontology    : "BP"/"MF"/"CC"/"ALL" (GO only; ignored for KEGG).
# rank_metric : label of the metric used, stored in params for reproducibility.
#
# Returns list(table, object, params) or empty_result(reason).
# -----------------------------------------------------------------------------
run_gsea <- function(rank_vec,
                     database    = c("GO", "KEGG"),
                     ontology    = "BP",
                     padj_method = "BH",
                     pval_cutoff = 0.05,
                     min_size    = 10,
                     max_size    = 500,
                     rank_metric = NA_character_,
                     use_cache   = TRUE) {

  database <- match.arg(database)

  if (is.null(rank_vec) || length(rank_vec) < 2) {
    return(empty_result(
      "GSEA needs a ranked list of genes. Run a differential-expression contrast first."))
  }
  rank_vec <- sort(rank_vec, decreasing = TRUE)   # enforce required ordering

  key <- .cache_key("gsea", database, ontology, padj_method, pval_cutoff,
                    min_size, max_size, names(rank_vec), unname(rank_vec),
                    ACTIVE_SPECIES)
  if (use_cache && exists(key, envir = .enrich_cache)) {
    message("GSEA: returning cached result")
    return(get(key, envir = .enrich_cache))
  }

  params <- list(
    database    = database,
    ontology    = if (identical(database, "GO")) ontology else NA_character_,
    padj_method = padj_method,
    pval_cutoff = pval_cutoff,
    min_size    = min_size,
    max_size    = max_size,
    rank_metric = rank_metric,
    n_genes     = length(rank_vec),
    org_db_pkg  = ORG_DB_PKG,
    key_type    = KEY_TYPE,
    kegg_code   = KEGG_CODE,
    species     = ACTIVE_SPECIES
  )
  gsea_warnings <- character(0)

  gsea_obj <- tryCatch({
    if (identical(database, "GO")) {
      # ── ID-recognition check, mirroring the KEGG branch below ───────────────
      # GO GSEA was the only enrichment path in the app that said nothing about
      # ID mapping. run_go_enrichment() reports query/background reconciliation,
      # run_kegg_enrichment() reports Entrez loss, and the KEGG branch a few
      # lines down reports its remap loss -- but if the ranked IDs were not keys
      # of ORG_DB at all, gseGO() simply found no gene sets and the panel
      # reported "No GO gene sets reached significance", which blames the biology
      # for what is actually an annotation mismatch.
      #
      # Reachable in several ordinary ways: StringTie/novel loci (MSTRG...),
      # transposable-element genes, spike-ins, a hand-pasted list, a rebuilt
      # OrgDb, or a species switch that left a stale ranked list behind (the
      # ranked vector is not invalidated by ACTIVE_SPECIES changing).
      #
      # Counted against keys(ORG_DB, KEY_TYPE) -- the same question the KEGG
      # branch asks with convert_to_entrez(): "is this ID in the ID space the
      # engine keys on?". A gene can be a valid key and still carry no GO term;
      # that is a separate, normal shortfall and is deliberately NOT claimed
      # here. tryCatch so an annotation that cannot list its keys degrades to the
      # previous behaviour rather than erroring.
      go_keys <- tryCatch(AnnotationDbi::keys(ORG_DB, keytype = KEY_TYPE),
                          error = function(e) NULL)
      if (!is.null(go_keys)) {
        n_ranked <- length(rank_vec)
        n_known  <- length(intersect(names(rank_vec), go_keys))
        params$n_genes_go_mapped <- n_known
        if (n_known == 0) {
          result <- empty_result(sprintf(paste0(
            "GO GSEA not run: none of the %d ranked gene(s) are %s keys in %s. ",
            "The ranked list and the active annotation are in different ID ",
            "spaces -- check the organism, or re-run the differential analysis ",
            "if you switched species after computing it."),
            n_ranked, KEY_TYPE, ORG_DB_PKG))
          result$params   <- params
          result$warnings <- gsea_warnings
          if (use_cache) assign(key, result, envir = .enrich_cache)
          return(result)
        }
        n_unknown   <- n_ranked - n_known
        pct_unknown <- if (n_ranked > 0) 100 * n_unknown / n_ranked else 0
        if (n_unknown > 0) {
          gsea_warnings <- c(gsea_warnings, sprintf(paste0(
            "%d of %d ranked gene(s) (%.1f%%) are not %s keys in %s and cannot ",
            "belong to any GO set. The enrichment score is computed over the %d ",
            "gene(s) the annotation recognises."),
            n_unknown, n_ranked, pct_unknown, KEY_TYPE, ORG_DB_PKG, n_known))
          if (pct_unknown > 20) {
            gsea_warnings <- c(gsea_warnings, paste0(
              "More than 20% of the ranking is unrecognised. GSEA assumes the ranked list ",
              "represents the whole measured transcriptome; a ranking missing this share of ",
              "it is no longer that list. Check that the gene IDs match the annotation."))
          }
          message(sprintf("GO GSEA: ID recognition = %d of %d ranked gene(s).",
                          n_known, n_ranked))
        }
      }
      clusterProfiler::gseGO(
        geneList      = rank_vec,
        OrgDb         = ORG_DB,
        keyType       = KEY_TYPE,
        ont           = ontology,
        minGSSize     = min_size,
        maxGSSize     = max_size,
        pvalueCutoff  = pval_cutoff,
        pAdjustMethod = padj_method,
        eps           = 0,        # accurate small p-values (fgsea multilevel)
        seed          = TRUE,     # reproducible permutations
        verbose       = FALSE
      )
    } else {
      remap <- .remap_rank_for_kegg(rank_vec)
      if (is.null(remap)) {
        result <- empty_result(paste0(
          "KEGG GSEA: none of the ", length(rank_vec),
          " ranked genes could be mapped to the ID type KEGG expects for ",
          "this species; KEGG GSEA skipped."))
        result$params   <- params
        result$warnings <- gsea_warnings
        if (use_cache) assign(key, result, envir = .enrich_cache)
        return(result)
      }
      kegg_rank <- remap$rank

      # ── Partial ID-mapping loss ─────────────────────────────────────────────
      # Only TOTAL failure was reported above. Partial loss was reported nowhere,
      # and params$n_genes records the PRE-remap ranked list, so a KEGG run whose
      # ranking lost most of its genes was indistinguishable in both the panel
      # and the reproducibility bundle from one that lost none.
      # Recorded as a separate field rather than by overwriting n_genes: that
      # field means "size of the ranked list" for GO runs and anything already
      # reading it must keep that meaning.
      params$n_genes_kegg_mapped <- remap$n_mapped
      n_unmapped   <- remap$n_query - remap$n_mapped
      pct_unmapped <- if (remap$n_query > 0) 100 * n_unmapped / remap$n_query else 0
      if (n_unmapped > 0) {
        map_src_gsea <- if (is.na(ID_MAP_SOURCE)) "the OrgDb" else ID_MAP_SOURCE
        gsea_warnings <- c(gsea_warnings, sprintf(paste0(
          "%d of %d ranked gene(s) (%.1f%%) have no Entrez/NCBI ID for %s and were ",
          "dropped from the ranking before KEGG GSEA. The enrichment score is ",
          "computed over the %d gene(s) that mapped. This is a coverage limit of ",
          "the ID source (%s), not a filter you chose."),
          n_unmapped, remap$n_query, pct_unmapped, SPECIES_COMMON,
          remap$n_mapped, map_src_gsea))
        if (pct_unmapped > 20) {
          gsea_warnings <- c(gsea_warnings, paste0(
            "More than 20% of the ranking has no KEGG-compatible ID. GSEA assumes ",
            "the ranked list represents the whole measured transcriptome; a ranking ",
            "missing this share of it is no longer that list, and the genes lost are ",
            "typically the less-characterised ones. Consider GO GSEA, which keeps ",
            "the full ranking."))
        }
        message(sprintf("KEGG GSEA: ID mapping = %d of %d ranked gene(s).",
                        remap$n_mapped, remap$n_query))
      }

      clusterProfiler::gseKEGG(
        geneList          = kegg_rank,
        organism          = KEGG_CODE,
        minGSSize         = min_size,
        maxGSSize         = max_size,
        pvalueCutoff      = pval_cutoff,
        pAdjustMethod     = padj_method,
        eps               = 0,
        seed              = TRUE,
        verbose           = FALSE,
        use_internal_data = FALSE
      )
    }
  }, error = function(e) {
    message("GSEA error: ", conditionMessage(e))
    NULL
  })

  if (is.null(gsea_obj) || nrow(gsea_obj@result) == 0) {
    result <- empty_result(sprintf(
      "No %s gene sets reached significance by GSEA (p.adjust <= %s) with set-size in [%d, %d].",
      database, pval_cutoff, min_size, max_size))
    result$params   <- params
    result$warnings <- gsea_warnings
    if (use_cache) assign(key, result, envir = .enrich_cache)
    return(result)
  }

  # GO core_enrichment stays in the user's own ID space (readable enough);
  # a NES-sorted table is the most useful default order.
  tbl <- as.data.frame(gsea_obj@result)
  tbl <- tbl[order(-abs(tbl$NES)), , drop = FALSE]

  result <- list(table = tbl, object = gsea_obj, params = params,
                 warnings = gsea_warnings)
  if (use_cache) assign(key, result, envir = .enrich_cache)
  result
}


# -----------------------------------------------------------------------------
# make_gsea_summary_plot(): compact NES summary across enriched sets. This is
# the GSEA analogue of make_enrichment_plot(); GSEA tables have NES/setSize
# instead of GeneRatio/Count, so it needs its own (small) plotter.
#
# style        : "bar" (NES bars) or "dot" (NES vs term, sized by setSize).
# color_by     : "significance" (gradient on p.adjust / -log10(p.adjust)),
#                "direction"    (flat colors split by NES sign — up vs down),
#                "flat"         (single user-chosen color, no legend).
# color_metric : "p.adjust" or "neglog10_padj" — only used when color_by = "significance".
# title, x_lab : custom text; x_lab defaults to the standard NES axis label
#                when NULL/"" so existing callers keep their old output.
# label_wrap   : character width to wrap long term Descriptions onto multiple
#                lines (0/NULL disables wrapping) — same convention as
#                make_enrichment_plot()'s label_wrap in enrichment_functions.R.
# font_family  : "sans"/"serif"/"mono" — base R graphics device fonts, no
#                extra font-management package (extrafont/showtext) required.
# theme_style  : "bw", "minimal", or "classic" ggplot2 base theme.
# -----------------------------------------------------------------------------
make_gsea_summary_plot <- function(result_df,
                                   top_n        = 20,
                                   font_size    = 13,
                                   font_family  = c("sans", "serif", "mono"),
                                   style        = c("bar", "dot"),
                                   color_by     = c("significance", "direction", "flat"),
                                   color_metric = c("p.adjust", "neglog10_padj"),
                                   color_low    = "#3B7DD8",
                                   color_high   = "#C0392B",
                                   up_color     = "#C0392B",
                                   down_color   = "#3B7DD8",
                                   flat_color   = "#3B7DD8",
                                   theme_style  = c("bw", "minimal", "classic"),
                                   title        = NULL,
                                   x_lab        = NULL,
                                   label_wrap   = 45) {
  style        <- match.arg(style)
  color_by     <- match.arg(color_by)
  color_metric <- match.arg(color_metric)
  font_family  <- match.arg(font_family)
  theme_style  <- match.arg(theme_style)

  blank <- function(msg) ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = msg, size = 6, color = "gray40") +
    theme_void()

  req_cols <- c("Description", "NES", "p.adjust", "setSize")
  if (is.null(result_df) || !all(req_cols %in% colnames(result_df)))
    return(blank("No GSEA results to plot with the current settings."))

  df <- result_df
  df$NES      <- suppressWarnings(as.numeric(df$NES))
  df$p.adjust <- suppressWarnings(as.numeric(df$p.adjust))
  df$setSize  <- suppressWarnings(as.numeric(df$setSize))
  df <- df[is.finite(df$NES) & is.finite(df$p.adjust), , drop = FALSE]
  if (nrow(df) == 0) return(blank("No GSEA results to plot."))

  df <- head(df[order(abs(df$NES), decreasing = TRUE), , drop = FALSE], top_n)
  df <- df[order(df$NES), , drop = FALSE]

  if (!is.null(label_wrap) && label_wrap > 0) {
    df$Description <- vapply(as.character(df$Description),
                             function(s) paste(strwrap(s, width = label_wrap), collapse = "\n"),
                             character(1))
  }
  df$Description <- factor(df$Description, levels = unique(df$Description))

  df$ColorMetric <- if (identical(color_metric, "neglog10_padj"))
    -log10(df$p.adjust) else df$p.adjust
  legend_name <- if (identical(color_metric, "neglog10_padj"))
    "-log10(p.adjust)" else "p.adjust"
  df$Direction <- ifelse(df$NES >= 0, "Up (positive NES)", "Down (negative NES)")

  p <- ggplot(df, aes(y = Description))

  if (identical(style, "bar")) {
    if (identical(color_by, "significance")) {
      p <- p + geom_col(aes(x = NES, fill = ColorMetric), width = 0.7) +
        scale_fill_gradient(low = color_low, high = color_high, name = legend_name)
    } else if (identical(color_by, "direction")) {
      p <- p + geom_col(aes(x = NES, fill = Direction), width = 0.7) +
        scale_fill_manual(values = c("Up (positive NES)" = up_color,
                                     "Down (negative NES)" = down_color),
                          name = "Direction")
    } else {
      p <- p + geom_col(aes(x = NES), fill = flat_color, width = 0.7)
    }
  } else {
    if (identical(color_by, "significance")) {
      p <- p + geom_point(aes(x = NES, size = setSize, color = ColorMetric)) +
        scale_color_gradient(low = color_low, high = color_high, name = legend_name) +
        scale_size_continuous(name = "Set size", range = c(3, 10))
    } else if (identical(color_by, "direction")) {
      p <- p + geom_point(aes(x = NES, size = setSize, color = Direction)) +
        scale_color_manual(values = c("Up (positive NES)" = up_color,
                                      "Down (negative NES)" = down_color),
                           name = "Direction") +
        scale_size_continuous(name = "Set size", range = c(3, 10))
    } else {
      p <- p + geom_point(aes(x = NES, size = setSize), color = flat_color) +
        scale_size_continuous(name = "Set size", range = c(3, 10))
    }
  }

  x_axis_label <- if (!is.null(x_lab) && nzchar(x_lab)) x_lab else "Normalized Enrichment Score (NES)"
  plot_title   <- if (!is.null(title) && nzchar(title)) title else ""

  p <- p +
    geom_vline(xintercept = 0, linewidth = 0.4, color = "grey40") +
    xlab(x_axis_label) + ylab(NULL) + ggtitle(plot_title)

  p <- p + switch(theme_style,
                  bw      = theme_bw(base_size = font_size, base_family = font_family),
                  minimal = theme_minimal(base_size = font_size, base_family = font_family),
                  classic = theme_classic(base_size = font_size, base_family = font_family))

  p + theme(axis.text.y  = element_text(size = font_size),
            axis.text.x  = element_text(size = font_size),
            axis.title.x = element_text(size = font_size + 1, face = "bold"),
            plot.title   = element_text(size = font_size + 3, face = "bold", hjust = 0.5),
            legend.text  = element_text(size = font_size - 1),
            legend.title = element_text(size = font_size, face = "bold"))
}
