# =============================================================================
# demo_module.R
# -----------------------------------------------------------------------------
# Example dataset and guided tour.
#
# No real count matrix is redistributable with the app, so a small, realistic
# one is generated at runtime and labelled synthetic wherever it appears. Gene
# IDs follow the real Solyc / AT patterns so canonicalization and link-outs
# behave exactly as they do on real data; the counts are simulated.
#
# The IDs are drawn at random from those patterns, so only a small fraction
# collide with real loci. Enrichment on the demo therefore returns nothing
# rather than plausible-looking terms for a random gene set.
#
# make_demo_counts() returns a data.frame (genes x samples) ready to feed into
# the same path as an uploaded CSV.
#
# The tour uses {cicerone} if available (graceful no-op otherwise), so you can
# add it without a hard dependency.
# =============================================================================

# --- Synthetic count matrix --------------------------------------------------
# 6 samples: two conditions (control/treated) x 3 replicates. ~1500 genes.
# A subset of genes is made clearly DE so the volcano/heatmap show structure.
make_demo_counts <- function(species = NULL, n_genes = 1500L, seed = 42L) {
  species <- species %||% (if (exists("ACTIVE_SPECIES")) ACTIVE_SPECIES else "tomato")
  set.seed(seed)

  gene_ids <- if (identical(species, "tomato")) {
    # Solyc##g######.# style, real-looking, synthetic values.
    sprintf("Solyc%02dg%06d.%d",
            sample(1:12, n_genes, replace = TRUE),
            sample(5000:99999, n_genes),
            sample(1:3, n_genes, replace = TRUE))
  } else {
    chrom <- sample(1:5, n_genes, replace = TRUE)
    sprintf("AT%dG%05d", chrom, sample(1000:99999, n_genes))
  }
  gene_ids <- make.unique(gene_ids, sep = "_")

  samples <- c("control_1", "control_2", "control_3",
               "treated_1", "treated_2", "treated_3")
  is_treated <- grepl("^treated", samples)

  # Baseline expression: log-normal means -> negative binomial counts.
  base_mu <- rlnorm(n_genes, meanlog = 4, sdlog = 1.4)
  disp    <- 0.2

  mat <- matrix(0L, nrow = n_genes, ncol = length(samples),
                dimnames = list(gene_ids, samples))

  # Make ~12% of genes DE (half up, half down in treated).
  n_de <- round(n_genes * 0.12)
  de_idx <- sample(seq_len(n_genes), n_de)
  lfc <- rep(0, n_genes)
  lfc[de_idx] <- sample(c(-1, 1), n_de, replace = TRUE) * runif(n_de, 1, 3)

  for (j in seq_along(samples)) {
    fold <- if (is_treated[j]) 2^lfc else rep(1, n_genes)
    mu   <- base_mu * fold
    size <- 1 / disp
    mat[, j] <- rnbinom(n_genes, mu = mu, size = size)
  }

  as.data.frame(mat, check.names = FALSE)
}

# --- Guided tour (optional, cicerone) ----------------------------------------
# Returns a Cicerone object if the package is present, else NULL. Steps target
# element IDs that exist in app.r; adjust selectors to taste.
build_demo_tour <- function() {
  if (!requireNamespace("cicerone", quietly = TRUE)) return(NULL)
  cicerone::Cicerone$
    new()$
    step("file", "1. Upload data",
         "Upload a CSV/Excel counts matrix here \u2014 or use the demo button.")$
    step("gene_filter", "2. Filter genes",
         "Optionally restrict to a gene list. Leave empty to show everything.")$
    step("group_ui", "3. Define groups",
         "Group replicates into conditions; each becomes one heatmap column.")$
    step("normalization", "4. Normalization",
         "Pick how values are scaled (z-score, log2, log2FC vs a reference).")$
    step("go_btn", "5. Generate heatmap",
         "Build the heatmap, then BRUSH a cluster to run enrichment on it.")$
    step("enrichment_panel", "6. Enrich your selection",
         "GO, KEGG, STRING and gene details for exactly the genes you selected.")
}

# Convenience server fragment: call from app server to wire a "Load demo" button
# with id `demo_btn` and a tour trigger `tour_btn`.
#   demoServer(input, session,
#              on_load = function(df) { ... place df where mat_data reads it ... })
demoServer <- function(input, output, session, on_load) {
  observeEvent(input$demo_btn, {
    df <- make_demo_counts()
    showNotification(
      paste0("\U0001f9ea Loaded SYNTHETIC demo: ", nrow(df),
             " genes \u00d7 ", ncol(df),
             " samples (control/treated \u00d7 3). Not real data."),
      type = "message", duration = 8)
    on_load(df)
  })

  tour <- build_demo_tour()
  observeEvent(input$tour_btn, {
    if (is.null(tour)) {
      showNotification(
        "Install the 'cicerone' package to enable the guided tour.",
        type = "warning", duration = 6)
    } else {
      tour$init()$start()
    }
  })
}

