# ── Normalization choices (shown in the UI selectInput) ───────────────────────
normalization_choices <- c(
  "Z-score per gene (row)"   = "zscore_row",
  "Z-score per sample (col)" = "zscore_col",
  "Log2 (counts+1)"          = "log2",
  "Log2 Z-score"             = "log2zscore",
  "Log2FC vs reference"      = "log2fc",
  "None (raw values)"        = "none"
)

# ── Helper: average raw counts across replicates per condition ────────────────
# mat    : numeric matrix (genes × samples) with replicate columns
# groups : named list; each element is a character vector of column names
#          belonging to one condition. Names become the output column names.
# Returns: numeric matrix (genes × conditions)
average_replicates <- function(mat, groups) {
  out <- do.call(cbind, lapply(names(groups), function(cond) {
    cols <- intersect(groups[[cond]], colnames(mat))
    if (length(cols) == 0) stop("No columns matched for condition: ", cond)
    rowMeans(mat[, cols, drop = FALSE], na.rm = TRUE)
  }))
  colnames(out) <- names(groups)
  out
}

# ── Log2FC normalisation vs reference condition ───────────────────────────────
# mat           : numeric matrix (genes × conditions), already averaged per condition
# ref_condition : character vector of column names to treat as reference
# Each gene is divided by the mean expression across all reference columns.
# Returns       : log2FC matrix (same dimensions as mat)
log2fc_vs_reference <- function(mat, ref_condition, already_log_scale = FALSE) {
  ref_cols <- intersect(ref_condition, colnames(mat))
  if (length(ref_cols) == 0)
    stop("No reference conditions found in the matrix columns.")
  ref_mean <- rowMeans(mat[, ref_cols, drop = FALSE], na.rm = TRUE)
  if (isTRUE(already_log_scale)) {
    # Input already on a log2-like scale (e.g. DESeq2 VST). A ratio of logs is
    # meaningless; the correct relative measure on a log scale is a difference:
    # log2FC = value - reference.
    mat - ref_mean
  } else {
    # Linear-scale input (TPM / normalized counts): classic ratio-based log2FC.
    log2((mat + 1) / (ref_mean + 1))
  }
}

# ── Main normalisation dispatcher ─────────────────────────────────────────────
# mat                : numeric matrix (genes × conditions), already averaged per condition
# method             : one of the keys in normalization_choices
# ref_condition      : character vector; required only when method = "log2fc"
# already_log_scale  : TRUE when `mat` has already been through a log-like
#                      variance-stabilizing transform upstream (e.g. DESeq2
#                      VST). VST output is already on an approximately log2
#                      scale; applying log2(mat + 1) again on top of it is a
#                      double transform that further compresses an
#                      already-compressed scale. When TRUE, the "log2" /
#                      "Log2 Z-score" options skip the extra log2() step.
#                      Default FALSE preserves existing behaviour for every
#                      other/future caller.
normalize_matrix <- function(mat, method = "zscore_row", ref_condition = NULL,
                             already_log_scale = FALSE) {
  switch(method,
         # NOTE: app.r does NOT route zscore_row/zscore_col/log2zscore here.
         # It standardises the SAMPLE matrix and averages afterwards, which is
         # not the same thing -- averaging shrinks variance, so a denominator
         # taken from group means differs from one taken from samples. These
         # branches are the group-mean variants, kept for any caller that only
         # has group means to work with. See normalization_substrate() below
         # and the "Normalise" block in app.r before changing either side.
         zscore_row = t(scale(t(mat))),
         zscore_col = scale(mat),
         log2 = {
           if (isTRUE(already_log_scale)) {
             mat  # already log-scale (e.g. VST) — do not re-log
           } else {
             # Shift to positive first if negative values are present (e.g. LFC input)
             min_val <- min(mat, na.rm = TRUE)
             if (min_val < 0) mat <- mat - min_val
             log2(mat + 1)
           }
         },
         log2zscore = {
           if (isTRUE(already_log_scale)) {
             t(scale(t(mat)))  # already log-scale — just z-score it
           } else {
             min_val <- min(mat, na.rm = TRUE)
             if (min_val < 0) mat <- mat - min_val
             t(scale(t(log2(mat + 1))))
           }
         },
         log2fc = {
           if (is.null(ref_condition) || length(ref_condition) == 0)
             stop("Log2FC vs reference requires at least one reference condition.")
           log2fc_vs_reference(mat, ref_condition, already_log_scale = already_log_scale)
         },
         none = mat,
         t(scale(t(mat)))  # default fallback
  )
}

# ── Heatmap row/column alignment helpers ──────────────────────────────────────
# These three functions each perform ONE of the synchronized drops that keep
# the clustering matrix (VST/expression-scale) and the display matrix
# (z-score/log2/etc.) row- and column-aligned after filtering. Each returns
# both matrices already filtered, plus the names of what was dropped, so the
# caller (app.r) only has to turn `dropped` into a showNotification() message
# -- no filtering logic lives in app.r anymore, so nothing there can drift out
# of sync with what's tested here.

# 1) Drop rows with NaN/Inf/NA (checked on the SCALED/display matrix, since
#    that's where e.g. a divide-by-zero in a display transform would surface).
drop_invalid_rows <- function(mat_raw, mat_scaled, gene_ids) {
  bad <- apply(mat_scaled, 1, function(x) any(is.nan(x) | is.infinite(x) | is.na(x)))
  list(
    raw      = mat_raw[!bad, , drop = FALSE],
    scaled   = mat_scaled[!bad, , drop = FALSE],
    gene_ids = gene_ids[!bad],
    dropped  = gene_ids[bad]
  )
}

# 2) Drop reference column(s) when method == "log2fc". mat_raw's VALUES are
#    left untouched (still absolute expression, not log2FC) -- only its
#    COLUMNS are dropped to match mat_scaled. This is what avoids the
#    shared-reference correlation artifact during clustering (see the
#    clustering-fix report): every log2FC display column shares the same
#    subtracted reference vector, which would otherwise inflate correlations
#    between conditions that only "look similar" because both are relative
#    to the same baseline.
drop_reference_columns <- function(mat_raw, mat_scaled, method, ref_condition) {
  out <- list(raw = mat_raw, scaled = mat_scaled,
              display_col_mask = rep(TRUE, ncol(mat_scaled)),
              dropped = character(0), error = NULL)
  if (!identical(method, "log2fc") || is.null(ref_condition) || length(ref_condition) == 0)
    return(out)
  
  ref_drop  <- intersect(ref_condition, colnames(mat_scaled))
  keep_cols <- setdiff(colnames(mat_scaled), ref_drop)
  if (length(keep_cols) == 0) {
    out$error <- "All columns are reference conditions \u2014 nothing to display."
    return(out)
  }
  out$display_col_mask <- colnames(mat_scaled) %in% keep_cols
  out$scaled  <- mat_scaled[, keep_cols, drop = FALSE]
  out$raw     <- mat_raw[, keep_cols, drop = FALSE]
  out$dropped <- ref_drop
  out
}

# 3) Drop zero-variance rows. Checked on BOTH matrices, dropping the union.
#    SCALED: a row constant after the chosen display transform is meaningless
#    to draw. RAW: a row constant on the clustering substrate has an undefined
#    Pearson correlation, which puts NA into 1 - cor(t(mat_raw)) and makes
#    hclust() error out. For every method except zscore_col the two sets are
#    identical; they differ only under zscore_col, where scale() per column
#    turns a raw-constant row into a non-constant one. Checking the raw matrix
#    as well is what stops such a row reaching the clustering call.
drop_zero_variance_rows <- function(mat_raw, mat_scaled, gene_ids, gene_labels) {
  if (ncol(mat_scaled) <= 1) {
    return(list(raw = mat_raw, scaled = mat_scaled,
                gene_ids = gene_ids, gene_labels = gene_labels,
                dropped = character(0)))
  }
  is_const <- function(m) apply(m, 1, function(x) {
    v <- var(x, na.rm = TRUE)
    is.na(v) || v < 1e-10
  })
  const_rows <- is_const(mat_scaled) | is_const(mat_raw)

  list(
    raw         = mat_raw[!const_rows, , drop = FALSE],
    scaled      = mat_scaled[!const_rows, , drop = FALSE],
    gene_ids    = gene_ids[!const_rows],
    gene_labels = gene_labels[!const_rows],
    dropped     = gene_ids[const_rows]
  )
}

# ---- The row-clustering contract, in one place ------------------------------
# Every heatmap in the app (main, cart, GO) and the k-suggestion preview must
# answer "should I cluster these rows, and how?" the same way, so this function
# is the single answer and every caller goes through it. Re-implementing the
# decision per panel invites three specific failures: accepting 2 columns where
# the main heatmap refuses them, falling through to ComplexHeatmap's default
# (euclidean, on the DISPLAY matrix) by not passing cluster_rows at all, and
# scoring a suggested k against a tree the heatmap will not cut.
#
# mat      : the CLUSTERING substrate (VST/expression scale). NEVER the display
#            matrix -- that is the whole point of the mean_mat_raw /
#            mean_mat_scaled split, and passing mat_scaled here defeats it.
# min_cols : below this, correlation distance is not informative. 3 by default:
#            over 2 points every Pearson correlation is forced to exactly +-1,
#            so the tree would reflect having 2 observations, not co-expression.
# linkage  : hclust() agglomeration method. Restricted to the four that behave
#            on a correlation distance; centroid/median are excluded because
#            they produce inversions (a merge lower than the merge beneath it),
#            which makes the dendrogram cross itself.
#
# Returns list(hc, d, reason):
#   hc     = an hclust object, or FALSE (drops straight into cluster_rows =)
#   d      = the dist object the tree was built from, or NULL. Returned so the
#            k-suggestion can score THIS tree with silhouette() instead of
#            rebuilding its own distance and drifting out of sync.
#   reason = NULL on success, else a user-facing sentence explaining the refusal
#
# Refuses rather than repairs: it never modifies `mat`, so the caller's matrix
# and the returned tree can never fall out of sync. FALSE is always safe --
# ComplexHeatmap simply draws rows in matrix order.
CLUSTER_LINKAGE_CHOICES <- c(
  "Complete \u2014 compact, similar-diameter clusters" = "complete",
  "Average (UPGMA) \u2014 conventional for expression"  = "average",
  "Ward.D2 \u2014 balanced cluster sizes"               = "ward.D2",
  "Single \u2014 chaining; diagnostic only"             = "single"
)

cluster_rows_or_reason <- function(mat, min_cols = 3L, linkage = "complete") {
  if (!linkage %in% CLUSTER_LINKAGE_CHOICES) linkage <- "complete"

  if (is.null(mat) || !is.matrix(mat) || !is.numeric(mat))
    return(list(hc = FALSE, d = NULL,
                reason = "No numeric clustering substrate available."))
  if (nrow(mat) < 2)
    return(list(hc = FALSE, d = NULL, reason = sprintf(
      "Only %d row(s) \u2014 nothing to cluster.", nrow(mat))))
  # "samples", not "conditions". Every call site that SHOWS this sentence --
  # app.r's go_btn observer and optimalKServer -- passes the sample-level
  # matrix, so the columns being counted are samples. The three sub-heatmap
  # call sites that can be handed condition means (GO heatmap, cart heatmap x2)
  # take $hc only and discard $reason, so they never display it. Saying
  # "conditions" named a unit the check was not counting: a 2-condition x
  # 3-replicate design has 6 columns here and is clustered normally.
  if (ncol(mat) < min_cols)
    return(list(hc = FALSE, d = NULL, reason = sprintf(
      paste0("Row clustering needs \u2265 %d samples to be meaningful ",
             "(correlation over %d point(s) is not informative)."),
      min_cols, ncol(mat))))
  if (anyNA(mat))
    return(list(hc = FALSE, d = NULL,
                reason = "Clustering substrate contains NA values."))

  # A constant row has sd = 0, so its Pearson correlation with anything is
  # undefined -> NA in the distance matrix -> hclust() errors. The main path
  # cannot reach this (drop_zero_variance_rows removes them upstream), but the
  # sub-heatmaps' display-matrix fallback can.
  const <- apply(mat, 1, function(x) {
    v <- var(x, na.rm = TRUE)
    is.na(v) || v < 1e-10
  })
  if (any(const))
    return(list(hc = FALSE, d = NULL, reason = sprintf(
      paste0("%d row(s) are constant across samples \u2014 Pearson ",
             "correlation is undefined for them."), sum(const))))

  d_raw <- tryCatch(1 - cor(t(mat), method = "pearson"), error = function(e) NULL)
  if (is.null(d_raw) || anyNA(d_raw))
    return(list(hc = FALSE, d = NULL,
                reason = "Correlation distance is undefined for this subset."))

  # ward.D2 squares the dissimilarities internally, which assumes they are
  # Euclidean distances. (1 - cor) is NOT a metric: it violates the triangle
  # inequality. sqrt(2 * (1 - cor)) IS the Euclidean embedding of correlation
  # distance, so Ward is applied to a proper metric. This branch cannot change
  # any pre-existing figure: ward.D2 is a new option, and complete/average/
  # single are invariant under any monotone transform of the distances.
  # To revert, delete this line and use `as.dist(d_raw)` unconditionally.
  d <- as.dist(if (identical(linkage, "ward.D2")) sqrt(2 * d_raw) else d_raw)

  hc <- tryCatch(hclust(d, method = linkage), error = function(e) NULL)
  if (is.null(hc))
    return(list(hc = FALSE, d = NULL, reason = "hclust() failed on this subset."))
  list(hc = hc, d = d, reason = NULL)
}

# ---- Which matrix a display transform is computed on -------------------------
# Two substrates are in play and they give different numbers:
#   "samples"     -- the transform is applied to the sample-level matrix and
#                    replicates are averaged AFTERWARDS. The sd therefore
#                    includes within-group (replicate) variance, which pulls in
#                    genes whose replicates disagree.
#   "group means" -- the transform is applied after replicates were averaged.
#
# Reported on the figure and in the reproducibility bundle because the two are
# indistinguishable from the plot itself, yet a gene can land on a visibly
# different colour depending on which was used.
normalization_substrate <- function(method) {
  if (identical(method, "zscore_row") || identical(method, "zscore_col") ||
      identical(method, "log2zscore")) "samples" else "group means"
}

# ---- Provenance caption printed on the figure --------------------------------
# The caption is drawn as the heatmap's column_title, so it ends up on every
# exported PDF, so every claim in it must be conditional on what actually ran.
# VST in particular is applied only when the user ticks "These are raw counts";
# asserting it unconditionally would make a TPM/CPM figure claim a transform
# that never happened. Kept as a separate function, rather than built inline in
# the go_btn observer, so the claim can be unit-tested instead of eyeballed.
#
# row_clustered / col_clustered : did each axis actually get a dendrogram?
# vst_applied                   : did DESeq2 VST run upstream?
# colour_substrate              : "samples" or "group means" -- see
#                                 normalization_substrate(). Omitted (NULL) the
#                                 colour claim stays as it was.
#
# row_substrate / col_substrate : the two axes are no longer clustered on the
#   same matrix -- rows use the row z-scored SAMPLE matrix, columns the
#   un-z-scored GROUP MEANS (see the two comment blocks in app.r's go_btn
#   observer for why the z-score deliberately stops at the row axis). Passing
#   one string per axis is what stops the caption from claiming they share a
#   substrate. Both default to the old single value, so any caller that does not
#   pass them gets exactly the previous sentence.
build_provenance_caption <- function(scale_label, row_clustered, linkage,
                                     col_clustered, n_genes, vst_applied,
                                     colour_substrate = NULL,
                                     row_substrate = NULL, col_substrate = NULL) {
  substrate <- if (isTRUE(vst_applied)) "VST" else "the uploaded values"
  row_sub <- if (is.null(row_substrate)) substrate else row_substrate
  col_sub <- if (is.null(col_substrate)) substrate else col_substrate
  paste0(
    "Colour: ", scale_label,
    if (!is.null(colour_substrate)) paste0(" (on ", colour_substrate, ")") else "",
    "  ·  Rows: ",
    if (!isTRUE(row_clustered)) "unclustered"
    else paste0("Pearson / ", linkage, " on ", row_sub),
    "  ·  Columns: ",
    if (!isTRUE(col_clustered)) "unclustered"
    else paste0("Spearman / complete on ", col_sub),
    "  ·  n = ", n_genes, " genes"
  )
}

# ---- Top-N most variable rows ------------------------------------------------
# mat : the CLUSTERING/expression substrate (VST scale), NOT the display matrix.
#       Ranking on a z-scored matrix is meaningless: every row has variance 1,
#       so "most variable" becomes an arbitrary tie and a cosmetic transform
#       silently decides which genes are shown.
# n   : how many rows to keep.
# Returns an integer vector of row indices, in ORIGINAL matrix order (not
# variance order) so the caller's row_order/clustering is unaffected by the
# filter -- it only removes rows, it never reorders them.
top_variable_rows <- function(mat, n) {
  if (is.null(mat) || !is.matrix(mat)) return(integer(0))
  n <- as.integer(n)
  if (is.na(n) || n <= 0 || n >= nrow(mat)) return(seq_len(nrow(mat)))
  v <- apply(mat, 1, var, na.rm = TRUE)
  v[is.na(v)] <- -Inf                      # constant/NA rows rank last
  sort(utils::head(order(v, decreasing = TRUE), n))
}

# ---- Cluster reproducibility (gene resampling) -------------------------------
# Hennig's clusterboot IDEA, but NOT his procedure -- and the difference matters
# when reporting the result. Draws are made with replacement, then de-duplicated
# on the next line (duplicate rows give zero correlation distance and break
# hclust). De-duplicating discards the multiplicities that make a draw a
# bootstrap, leaving a SUBSAMPLE of expected size 1 - 1/e = ~63.2% of the genes.
# So: this is subsampling without replacement at ~63%, not a bootstrap, and the
# word "bootstrap" is deliberately not used for it any more.
#
# Consequence for interpretation: Hennig's rules of thumb (> 0.75 stable,
# < 0.6 not credible) are calibrated for resampling WITH replacement and do not
# transfer as absolute verdicts. The values remain a sound RELATIVE indicator --
# a cluster carried by a handful of genes really does score lower than one
# driven by broad co-expression -- and the UI presents them that way.
#
# Genes are resampled, not conditions. Resampling conditions is what pvclust
# does, but these datasets typically have 3-6 conditions -- resampling 3
# columns produces mostly degenerate replicates. Resampling genes works at any
# number of conditions and answers the question that actually matters here: is
# this partition driven by the whole gene set, or by a handful of genes?
#
# For each original cluster C, and each replicate, we find the replicate cluster
# with the highest Jaccard overlap with C, then average that maximum over
# replicates.
#
# Seeded by default (seed = 1L). Without a seed the same partition scores
# differently on every call and a cluster sitting near a band edge changes
# verdict between runs, which reads as instability in the data rather than in
# the estimate. Pass seed = NULL for a free-running stream.
#
# mat      : clustering substrate (genes x conditions), rows aligned to `assign`
# assign   : character/numeric vector, one cluster label per row of `mat`
# k        : number of clusters to cut in each replicate
# linkage  : must match what produced `assign`
# B        : bootstrap replicates
# Returns a data.frame(cluster, n_genes, stability) or NULL if it cannot run.
cluster_stability <- function(mat, assign, k, linkage = "complete", B = 50L,
                              seed = 1L) {
  if (is.null(mat) || !is.matrix(mat) || is.null(assign)) return(NULL)
  if (nrow(mat) != length(assign)) return(NULL)
  k <- as.integer(k); B <- as.integer(B)
  if (is.na(k) || k < 2 || is.na(B) || B < 1) return(NULL)

  # ── Seeding ───────────────────────────────────────────────────────────────
  # Without a seed the same partition scores differently on every click, and a
  # cluster sitting near a band edge changes verdict between runs. seed = NULL
  # selects a free-running stream for callers that want one.
  #
  # The ambient RNG is saved and put back afterwards. A bare set.seed() would
  # re-point the whole session's random stream, so any later unrelated draw
  # would depend on whether the user happened to press this button -- trading
  # one reproducibility problem for a subtler one.
  #
  # base::assign() is qualified because this function's own `assign` argument
  # shadows it; an unqualified call here would try to invoke a character vector.
  if (!is.null(seed)) {
    if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      .prev_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
      on.exit(base::assign(".Random.seed", .prev_seed, envir = globalenv()),
              add = TRUE)
    } else {
      on.exit(suppressWarnings(rm(".Random.seed", envir = globalenv())),
              add = TRUE)
    }
    set.seed(seed)
  }

  assign   <- as.character(assign)
  clusters <- sort(unique(assign))
  jac_sum  <- setNames(numeric(length(clusters)), clusters)
  # Running sum of squares, so the spread of the per-replicate Jaccard values
  # can be reported alongside their mean. Storing the squares rather than every
  # value keeps the accumulator the same shape as jac_sum / jac_n.
  jac_sq   <- setNames(numeric(length(clusters)), clusters)
  jac_n    <- setNames(integer(length(clusters)), clusters)

  for (b in seq_len(B)) {
    idx <- sample.int(nrow(mat), nrow(mat), replace = TRUE)
    idx <- unique(idx)                    # duplicate rows -> zero distance; drop
    if (length(idx) < k + 1L) next
    cl_b <- cluster_rows_or_reason(mat[idx, , drop = FALSE], linkage = linkage)
    if (isFALSE(cl_b$hc)) next
    part_b <- cutree(cl_b$hc, k = min(k, length(idx) - 1L))

    for (cid in clusters) {
      orig <- which(assign == cid)
      orig_in_b <- intersect(orig, idx)   # the cluster as seen by this replicate
      if (length(orig_in_b) == 0) next
      best <- 0
      for (bid in unique(part_b)) {
        cand <- idx[part_b == bid]
        inter <- length(intersect(orig_in_b, cand))
        if (inter == 0) next
        uni <- length(union(orig_in_b, cand))
        if (uni > 0) best <- max(best, inter / uni)
      }
      jac_sum[cid] <- jac_sum[cid] + best
      jac_sq[cid]  <- jac_sq[cid]  + best^2
      jac_n[cid]   <- jac_n[cid] + 1L
    }
  }

  if (all(jac_n == 0)) return(NULL)

  # Mean and standard error of the per-replicate Jaccard values. The variance
  # uses the n-1 denominator and is floored at 0, because the running-sums
  # identity can return a tiny negative value through rounding when every
  # replicate scored identically. NA when a cluster survived fewer than two
  # replicates -- there is no spread to report from one number.
  mean_j <- ifelse(jac_n > 0, jac_sum / pmax(jac_n, 1L), NA_real_)
  var_j  <- ifelse(jac_n > 1,
                   pmax(0, (jac_sq - jac_n * mean_j^2) / pmax(jac_n - 1L, 1L)),
                   NA_real_)
  se_j   <- ifelse(jac_n > 1, sqrt(var_j / pmax(jac_n, 1L)), NA_real_)

  out <- data.frame(
    cluster    = clusters,
    n_genes    = as.integer(table(factor(assign, levels = clusters))),
    stability  = ifelse(jac_n > 0, round(mean_j, 3), NA_real_),
    se         = round(se_j, 4),
    # How many of the B replicates actually contained this cluster. Reported
    # because it is the denominator of the standard error: a cluster scored on
    # 12 of 50 replicates carries a far weaker estimate than one scored on 50,
    # and the mean alone does not show that.
    n_resamples = as.integer(jac_n),
    stringsAsFactors = FALSE
  )
  # Recorded on the result so the run can be reproduced and cited. Attributes
  # rather than columns: the table is rendered straight to the UI and these are
  # provenance, not per-cluster measurements.
  attr(out, "seed")    <- seed
  attr(out, "B")       <- B
  attr(out, "linkage") <- linkage
  out
}
