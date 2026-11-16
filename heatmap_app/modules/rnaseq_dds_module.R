# =============================================================================
# rnaseq_dds_module.R
# -----------------------------------------------------------------------------
# Shiny module responsible for the first stage of the RNA-seq pipeline:
#   1. Uploading and parsing the raw counts matrix (CSV or Excel).
#   2. Building or manually creating the sample metadata table.
#   3. Renaming metadata columns before saving.
#   4. Constructing and fitting the DESeq2 object (DESeqDataSet).
#
# This module is the entry point of the RNA-seq analysis.
# It exposes four reactive values to downstream modules:
#   $dds                  — fitted DESeqDataSet, NULL until built
#   $metadata_saved       — committed metadata table
#   $metadata_for_analysis — metadata with live column renames applied
#   $counts               — raw counts data frame
# =============================================================================


# -----------------------------------------------------------------------------
# UI
# rnaseqDdsUI(): sidebar panel with file upload, metadata controls,
#                table editor, and DESeq2 design/build section.
# -----------------------------------------------------------------------------
rnaseqDdsUI <- function(id, start_step = 2) {
ns <- NS(id)

# Step numbers for the three stages this module owns. They are passed in rather
# than hard-coded because the module knows it is a sequence, but not where that
# sequence starts -- in the app, step 1 (Choose species) lives in app.r's shared
# well panel, above both features.
s_upload <- start_step         # 2
s_meta   <- start_step + 1L    # 3
s_design <- start_step + 2L    # 4

# Everything from the input-type selector down to the DESeq2 design panel folds
# into one block, and the server closes it automatically once the DESeq2 object
# has been built (see the runjs call at the end of the build_rnaseq_dds
# observer). By that point the upload, the metadata table and the design are all
# settled, and they occupy most of the sidebar above the contrast controls the
# user is about to reach for.
#
# Closed, not removed: re-opening it is one click, and every input inside keeps
# its value because <details> only hides its content, it does not unbind it.
#
# The id has to be namespaced -- two modules could otherwise collide on it.
collapsible_panel(
"🧬 RNAseq Metadata",
id = ns("meta_panel"), heading = "h4", open = TRUE,

step_header(s_upload, "Upload counts"),

# [SALMON-UI-1] Input type selector
radioButtons(
  ns("input_data_type"),
  "Input data type",
  choices = c(
    "Raw counts (CSV / Excel)"       = "raw",
    "Salmon output (quant.sf files)" = "salmon"
  ),
  selected = "raw"
),
tags$hr(),

# --- File upload ------------------------------------------------------------
# Accepts CSV or Excel files. First column is treated as gene ID row names.
# [SALMON-UI-2] Raw counts panel — visible only when input_data_type == "raw"
conditionalPanel(
  condition = paste0("input['", ns("input_data_type"), "'] == 'raw'"),
  fileInput(
    ns("file_rnaseq"),
    "Upload counts matrix (CSV or Excel)",
    accept = c(".csv", ".xlsx", ".xls")
  ),
  helpText("Sample names must be separated with '_' and must end with the replicate identifier."),
  helpText("Supported examples: condition_1, genotype_treatment_2, genotype_treatment_batch_3")
),

# [SALMON-UI-3] Salmon upload panel — visible only when input_data_type == "salmon"
conditionalPanel(
  condition = paste0("input['", ns("input_data_type"), "'] == 'salmon'"),
  div(
    class = "alert alert-info",
    style = "font-size: 0.88em;",
    tags$b("📦 Salmon import:"),
    tags$ul(
      tags$li("Upload all quant.sf files (one per sample). Name each file so encode the sample ID, e.g. WT_1_quant.sf, mut_2_quant.sf."),
      tags$li("Upload a tx2gene table: two-column CSV/TSV with transcript IDs (col 1) and gene IDs (col 2).")
    )
  ),
  fileInput(
    ns("file_quant_sf"),
    "Upload quant.sf files (one per sample)",
    accept   = c(".sf", ".tsv", ".txt"),
    multiple = TRUE
  ),
  fileInput(
    ns("file_tx2gene"),
    "Upload tx2gene table (CSV or TSV)",
    accept = c(".csv", ".tsv", ".txt")
  ),
  selectInput(
    ns("quant_type"),
    "Quantifier type",
    choices  = c("Salmon (quant.sf)" = "salmon", "Alevin" = "alevin"),
    selected = "salmon",
    width    = "100%"
  )
),

tags$hr(),
step_header(s_meta, "Build & save metadata"),

# --- Metadata mode ----------------------------------------------------------
# "separate"  → parses sample names into two independent factor columns
# "combined"  → creates a single combined_factor column from the full name
radioButtons(
ns("rnaseq_metadata_mode"),
lbl("Metadata mode", "How sample names become factors (cut on '_', last piece = replicate). Two factors: T22_mock_R1 becomes T22 + mock, which lets you fit ~ a + b and test an interaction. Combined: one group, 'T22_mock'. The table stays editable either way."),
choices = c(
"Two separate factors" = "separate",
"Single combined factor" = "combined"
),
selected = "separate"
),

# --- Metadata generation buttons --------------------------------------------
# "Build from sample names": auto-parses sample names → calls validate + build helpers
# "Create manual table":     produces an empty editable table for manual entry
fluidRow(
column(
6,
actionButton(
ns("build_rnaseq_metadata"),
"Build from sample names",
class = "btn-primary btn-sm",
width = "100%"
)
),
column(
6,
actionButton(
ns("create_manual_metadata"),
"Create manual table",
class = "btn-default btn-sm",
width = "100%"
)
)
),

# --- metadata file upload ----------------------------------------------
# Allows importing a ready-made metadata table from CSV or Excel.
# The uploaded table is displayed in the same editor used by the existing
# workflow and can then be confirmed with "Save this table".
fileInput(
  ns("file_rnaseq_metadata"),
  "Upload metadata table (CSV or Excel)",
  accept = c(".csv", ".xlsx", ".xls")
),
helpText("The metadata file must contain a 'sample' column matching the count matrix column names."),

br(), 
br(),

# --- Dynamic UI placeholders ------------------------------------------------
uiOutput(ns("rnaseq_metadata_messages")),  # validation alert banners
uiOutput(ns("rnaseq_colname_editor")),     # column rename text inputs
uiOutput(ns("rnaseq_metadata_ui")),        # DT table editor

# --- Download + Save buttons (visible only when metadata exists) ------------
conditionalPanel(
condition = paste0("output['", ns("rnaseq_metadata_ready"), "']"),
fluidRow(
column(4, downloadButton(ns("download_rnaseq_metadata_csv"), "Download CSV")),
column(4, downloadButton(ns("download_rnaseq_metadata_xlsx"), "Download Excel")),
column(
4,
actionButton(
ns("save_rnaseq_metadata"),
"Save this table",
class = "btn-success btn-sm",
width = "100%"
)
)
)
),

# --- DESeq2 design panel (visible only after metadata is saved) -------------
uiOutput(ns("rnaseq_design_panel"))
)   # end collapsible_panel("RNAseq Metadata") -- it replaced the old tagList()
}


# -----------------------------------------------------------------------------
# SERVER
# rnaseqDdsServer(): manages the full upstream RNA-seq workflow.
# No reactive arguments — all state is internal.
# Returns a named list of four reactives for use by downstream modules.
# -----------------------------------------------------------------------------
rnaseqDdsServer <- function(id, start_step = 2) {
moduleServer(id, function(input, output, session) {
ns <- session$ns

# Step 4 of the sidebar. Lives here and not in rnaseqDdsUI() because the design
# panel is rendered server-side (it only appears once metadata is saved), so the
# same offset the UI uses has to be available on this side too.
s_design <- start_step + 2L

# ---------------------------------------------------------------------------
# Reactive state
# ---------------------------------------------------------------------------
rnaseq_metadata_rv       <- reactiveVal(NULL)  # working metadata (editable)
rnaseq_validation_rv     <- reactiveVal(NULL)  # validation results from sample name check
rnaseq_metadata_saved_rv <- reactiveVal(NULL)  # committed metadata (locked for analysis)
rnaseq_dds_rv            <- reactiveVal(NULL)  # fitted DESeqDataSet
rnaseq_colname_map_rv    <- reactiveVal(NULL)  # named vector: original → renamed columns
# [SALMON-SRV-1] Stores the structured result from import_raw_counts()
# or import_salmon_counts(). Downstream code reads $type, $counts, $txi.
import_result_rv <- reactiveVal(NULL)

# [SALMON-SRV-3a] Trigger raw import when a counts file is uploaded.
observeEvent(input$file_rnaseq, {
  req(input$file_rnaseq)
  req(input$input_data_type == "raw")
  
  result <- tryCatch(
    import_raw_counts(
      file_path = input$file_rnaseq$datapath,
      file_name = input$file_rnaseq$name
    ),
    error = function(e) {
      showNotification(
        paste0("❌ Count file import failed: ", conditionMessage(e)),
        type = "error", duration = 10
      )
      NULL
    }
  )
  import_result_rv(result)
  if (!is.null(result)) {
    # ── Surface the validation warnings ────────────────────────────────────
    # import_raw_counts() attaches these as attr(counts, "validation_warnings")
    # and, until now, nothing read them anywhere in the app -- so the RNA-seq
    # path computed them and threw them away, while the heatmap path showed the
    # same warnings via notify_and_log_validation().
    #
    # This is where they matter most. DESeq2 requires UN-normalized counts: it
    # models counting noise and estimates its own size factors. Hand it CPM or
    # already-normalized counts and it does not complain -- size factors come
    # out ~1, dispersions are wrong, and every p-value is invalid. The check
    # that detects exactly that ("library sizes differ by only 1.0x") lives in
    # validate_count_matrix() and was being discarded.
    for (w in attr(result$counts, "validation_warnings") %||% character(0)) {
      showNotification(paste0("⚠️ ", w), type = "warning", duration = 12)
      log_warning(w)
    }
    showNotification(
      paste0("✅ Counts matrix loaded: ",
             nrow(result$counts), " genes × ",
             length(result$sample_names), " samples."),
      type = "message", duration = 5
    )
  }
})

# [SALMON-SRV-3b] Trigger Salmon import when quant.sf AND tx2gene are uploaded.
observeEvent(
  list(input$file_quant_sf, input$file_tx2gene),
  {
    req(input$input_data_type == "salmon")
    req(input$file_quant_sf)
    req(input$file_tx2gene)
    
    tx2gene_ext <- tools::file_ext(input$file_tx2gene$name)
    tx2gene_sep <- if (tolower(tx2gene_ext) == "csv") "," else "\t"
    
    quant_paths <- tryCatch(
      parse_salmon_uploads(input$file_quant_sf),
      error = function(e) {
        showNotification(
          paste0("❌ Salmon file names could not be parsed: ", conditionMessage(e)),
          type = "error", duration = 10
        )
        NULL
      }
    )
    req(!is.null(quant_paths))
    
    showNotification(
      paste0("⏳ Running tximport on ", length(quant_paths), " samples…"),
      type = "message", id = "txi_msg", duration = NULL
    )
    
    result <- tryCatch(
      import_salmon_counts(
        quant_paths  = quant_paths,
        tx2gene_path = input$file_tx2gene$datapath,
        quant_type   = input$quant_type,
        tx2gene_sep  = tx2gene_sep
      ),
      error = function(e) {
        removeNotification("txi_msg")
        showNotification(
          paste0("❌ tximport failed: ", conditionMessage(e)),
          type = "error", duration = 12
        )
        NULL
      }
    )
    
    removeNotification("txi_msg")
    import_result_rv(result)
    
    if (!is.null(result)) {
      showNotification(
        paste0("✅ Salmon import complete: ",
               nrow(result$txi$counts), " genes × ",
               length(result$sample_names), " samples."),
        type = "message", duration = 6
      )
    }
  },
  ignoreNULL = TRUE,
  ignoreInit = TRUE
)

# [SALMON-SRV-2] mat_data_rnaseq: unified reactive exposing a data.frame of
# counts regardless of input type. Used only for metadata building and sample
# name extraction — NOT for constructing the DESeqDataSet.
mat_data_rnaseq <- reactive({
  res <- import_result_rv()
  req(!is.null(res))
  if (res$type == "raw") {
    res$counts
  } else {
    as.data.frame(res$txi$counts)
  }
})

# ---------------------------------------------------------------------------
# rnaseq_sample_names: extracts column names from the counts matrix.
# These are the sample identifiers used to auto-build the metadata table.
# ---------------------------------------------------------------------------
rnaseq_sample_names <- reactive({
req(mat_data_rnaseq())
colnames(mat_data_rnaseq())
})

# ---------------------------------------------------------------------------
# imported_rnaseq_metadata: reads a user-supplied metadata table.
# CSV and Excel are supported. The table must contain a 'sample' column.
# ---------------------------------------------------------------------------
imported_rnaseq_metadata <- reactive({
  req(input$file_rnaseq_metadata)
  
  ext <- tolower(tools::file_ext(input$file_rnaseq_metadata$name))
  
  meta <- if (ext == "csv") {
    read.csv(input$file_rnaseq_metadata$datapath, check.names = FALSE)
  } else {
    as.data.frame(readxl::read_excel(input$file_rnaseq_metadata$datapath))
  }
  
  meta
})

# ---------------------------------------------------------------------------
# observeEvent: file_rnaseq_metadata
# Loads a metadata table provided by the user and reuses the current editor/
# save workflow. The table is validated before being exposed to the UI.
# ---------------------------------------------------------------------------
observeEvent(input$file_rnaseq_metadata, {
  req(imported_rnaseq_metadata())
  req(rnaseq_sample_names())
  
  meta <- imported_rnaseq_metadata()
  
  if (!"sample" %in% colnames(meta)) {
    showNotification(
      "The metadata file must contain a 'sample' column.",
      type = "error",
      duration = 8
    )
    return()
  }
  
  meta$sample <- as.character(meta$sample)
  
  if (!setequal(meta$sample, rnaseq_sample_names())) {
    showNotification(
      "The metadata file sample names do not match the count matrix column names exactly.",
      type = "error",
      duration = 10
    )
    return()
  }
  
  # Reorder metadata rows to match the count matrix column order.
  meta <- meta[match(rnaseq_sample_names(), meta$sample), , drop = FALSE]
  
  # Convert all non-sample columns to factors.
  # This ensures categorical design variables imported from CSV/Excel
  # are treated as factors by DESeq2.
  non_sample_cols <- setdiff(colnames(meta), "sample")
  if (length(non_sample_cols) > 0) {
    meta[non_sample_cols] <- lapply(meta[non_sample_cols], as.character)
  }
  
  rnaseq_validation_rv(NULL)
  rnaseq_metadata_rv(meta)
  rnaseq_metadata_saved_rv(NULL) # reset downstream state
  rnaseq_dds_rv(NULL)
  rnaseq_colname_map_rv(stats::setNames(colnames(meta), colnames(meta)))
  
  showNotification(
    "Metadata table loaded from file. Review it and click 'Save this table' to use it for the analysis.",
    type = "message",
    duration = 6
  )
}, ignoreNULL = TRUE)

# ---------------------------------------------------------------------------
# observeEvent: build_rnaseq_metadata
# Triggered by the "Build from sample names" button.
# Validates sample name format and auto-parses them into a metadata table.
# Resets saved metadata and dds to NULL to prevent stale downstream state.
# ---------------------------------------------------------------------------
observeEvent(input$build_rnaseq_metadata, {
req(rnaseq_sample_names())

validation <- validate_rnaseq_sample_names(rnaseq_sample_names())
meta <- build_rnaseq_metadata_table(
sample_names = rnaseq_sample_names(),
mode = input$rnaseq_metadata_mode
)

rnaseq_validation_rv(validation)
rnaseq_metadata_rv(meta)
rnaseq_metadata_saved_rv(NULL)  # reset downstream state
rnaseq_dds_rv(NULL)
rnaseq_colname_map_rv(stats::setNames(colnames(meta), colnames(meta)))

showNotification(
"RNA-seq metadata table generated from sample names.",
type = "message",
duration = 4
)
})

# ---------------------------------------------------------------------------
# observeEvent: create_manual_metadata
# Triggered by the "Create manual table" button.
# Produces an empty table with sample names already filled in.
# Useful when sample names do not follow a parseable convention.
# ---------------------------------------------------------------------------
observeEvent(input$create_manual_metadata, {
req(rnaseq_sample_names())

meta <- build_empty_rnaseq_metadata_table(
sample_names = rnaseq_sample_names(),
mode = input$rnaseq_metadata_mode
)

rnaseq_validation_rv(NULL)
rnaseq_metadata_rv(meta)
rnaseq_metadata_saved_rv(NULL)  # reset downstream state
rnaseq_dds_rv(NULL)
rnaseq_colname_map_rv(stats::setNames(colnames(meta), colnames(meta)))

showNotification(
"Empty RNA-seq metadata table created for manual editing.",
type = "message",
duration = 4
)
})

# ---------------------------------------------------------------------------
# rnaseq_current_colname_map: builds a named vector of original → current names.
# Reads text inputs from the column rename editor in real time.
# Falls back to the original column name if the input is empty or whitespace.
# The "sample" column is always protected and returned unchanged.
# ---------------------------------------------------------------------------
rnaseq_current_colname_map <- reactive({
req(rnaseq_metadata_rv())

cols <- colnames(rnaseq_metadata_rv())

stats::setNames(
vapply(cols, function(col) {
if (col == "sample") return("sample")
val <- input[[paste0("rnaseq_colname_", col)]]
if (is.null(val) || !nzchar(trimws(val))) col else trimws(val)
}, character(1)),
cols
)
})

# ---------------------------------------------------------------------------
# rnaseq_metadata_for_analysis: applies the current column name map to the
# working metadata. Used everywhere the renamed column names are needed
# (table display, downloads, design formula, DESeq2 object).
# ---------------------------------------------------------------------------
rnaseq_metadata_for_analysis <- reactive({
req(rnaseq_metadata_rv())
req(rnaseq_current_colname_map())

meta <- rnaseq_metadata_rv()
new_names <- unname(rnaseq_current_colname_map()[colnames(meta)])
colnames(meta) <- new_names
meta
})

# ---------------------------------------------------------------------------
# default_rnaseq_design: proposes a DESeq2 design formula automatically.
# Priority: combined_factor > factor1 + factor2 > first available column.
# At most two terms are included to avoid overfitting on small datasets.
# ---------------------------------------------------------------------------
# rnaseq_design_terms: the de-aliased factor(s) that go into the formula.
# Shared by default_rnaseq_design(), rnaseq_interaction_ok() and the
# interaction checkbox observer, so all three always agree on the terms.
rnaseq_design_terms <- reactive({
  req(rnaseq_metadata_for_analysis())
  
  cols <- setdiff(colnames(rnaseq_metadata_for_analysis()), "sample")
  preferred <- c("combined_factor", "factor1", "factor2", "observations")
  
  original_to_current <- rnaseq_current_colname_map()
  current_preferred <- unname(
    original_to_current[preferred[preferred %in% names(original_to_current)]]
  )
  current_preferred <- current_preferred[current_preferred %in% cols]
  
  selected <- current_preferred
  
  # combined_factor is the concatenation of factor1/factor2, so it is aliased
  # with either of them: `~ combined_factor + factor1` has linearly dependent
  # model-matrix columns and DESeq() stops with "model matrix is not full rank".
  # Drop it whenever a component factor is also available; keep it when it is
  # the only preferred column present (then it is not aliased with anything).
  if (length(selected) > 1) {
    cf <- unname(original_to_current["combined_factor"])
    if (length(cf) == 1 && !is.na(cf) && cf %in% selected) {
      selected <- setdiff(selected, cf)
    }
  }
  
  if (length(selected) == 0) {
    selected <- cols[1]
  } else if (length(selected) > 2) {
    selected <- selected[1:2]
  }
  
  selected
})

# ---------------------------------------------------------------------------
# default_rnaseq_design: ADDITIVE by default, deliberately.
# Under `~ f1 + f2`, results(contrast = c("f1","B","A")) is the f1 effect
# AVERAGED over f2. Under `~ f1 * f2` the SAME UI selection returns the f1
# effect at the REFERENCE LEVEL of f2 only. Making `*` the default would
# silently change the meaning of every saved contrast, so the interaction is
# opt-in via the checkbox below instead.
# ---------------------------------------------------------------------------
default_rnaseq_design <- reactive({
  selected <- rnaseq_design_terms()
  if (length(selected) == 2) {
    paste("~", paste(selected, collapse = " + "))
  } else {
    paste("~", selected[1])
  }
})

# ---------------------------------------------------------------------------
# rnaseq_interaction_ok: TRUE only when an interaction model is actually
# fittable. Every cell of the f1 x f2 table needs >= 2 samples: an empty cell
# makes the model rank-deficient, a singleton cell leaves no residual df.
# ---------------------------------------------------------------------------
rnaseq_interaction_ok <- reactive({
  selected <- rnaseq_design_terms()
  if (length(selected) != 2) return(FALSE)
  meta_now <- rnaseq_metadata_for_analysis()
  cells <- tryCatch(
    table(as.character(meta_now[[selected[1]]]),
          as.character(meta_now[[selected[2]]])),
    error = function(e) NULL
  )
  !is.null(cells) && length(cells) > 0 && all(cells >= 2)
})

# ---------------------------------------------------------------------------
# output$rnaseq_reference_ui: one selector per metadata column (except sample)
# to choose the reference level used for DESeq2 factor coding.
# Only columns with at least two distinct values are shown.
# ---------------------------------------------------------------------------
output$rnaseq_reference_ui <- renderUI({
  req(rnaseq_metadata_saved_rv())
  
  meta <- rnaseq_metadata_saved_rv()
  cols <- setdiff(colnames(meta), c("sample", "replicate"))
  
  valid_cols <- cols[sapply(cols, function(x) {
    vals <- unique(as.character(meta[[x]]))
    vals <- vals[!is.na(vals) & nzchar(vals)]
    length(vals) >= 2
  })]
  
  if (length(valid_cols) == 0) return(NULL)
  
  tagList(
    tags$b("Control levels"),
    tags$p(
      style = "font-size: 0.85em; color: #666; margin-bottom: 6px;",
      "Choose the control level for each factor before building the DESeq2 object."
    ),
    lapply(valid_cols, function(col) {
      lvls <- unique(as.character(meta[[col]]))
      lvls <- lvls[!is.na(lvls) & nzchar(lvls)]

      selectInput(
        # inputId stays ref_level_* -- it is read back by the build observer and
        # by contrast_params()$reference_levels, which the bundle records.
        inputId = ns(paste0("ref_level_", col)),
        label = paste("Control for", col),
        choices = lvls,
        selected = lvls[1],
        width = "100%"
      )
    })
  )
})

# ---------------------------------------------------------------------------
# observeEvent: save_rnaseq_metadata
# Commits the current metadata state to rnaseq_metadata_saved_rv.
# Validates that column names are unique before saving.
# Once saved, the DESeq2 design panel becomes visible.
# ---------------------------------------------------------------------------
observeEvent(input$save_rnaseq_metadata, {
req(rnaseq_metadata_rv())
req(rnaseq_metadata_for_analysis())
req(rnaseq_current_colname_map())

final_meta <- rnaseq_metadata_for_analysis()

if (anyDuplicated(colnames(final_meta)) > 0) {
showNotification(
"Column names must be unique before saving the table.",
type = "error",
duration = 8
)
return()
}

rnaseq_colname_map_rv(rnaseq_current_colname_map())
rnaseq_metadata_saved_rv(final_meta)

showNotification(
"RNA-seq metadata table saved for future analysis.",
type = "message",
duration = 4
)
})

# ---------------------------------------------------------------------------
# output$rnaseq_metadata_messages: validation alert banners.
# Three conditions are checked by validate_rnaseq_sample_names():
#   - missing_separator: sample names that do not contain '_'
#   - too_short: names that cannot be split into factor + replicate
#   - empty_replicate: names where the replicate token is empty
# All warnings are shown simultaneously. Green success if all pass.
# ---------------------------------------------------------------------------
output$rnaseq_metadata_messages <- renderUI({
validation <- rnaseq_validation_rv()
if (is.null(validation)) return(NULL)

msgs <- list()

if (length(validation$missing_separator) > 0) {
msgs[[length(msgs) + 1]] <- div(
class = "alert alert-warning",
tags$b("Warning: "),
paste0(
"The following sample names do not contain '_': ",
paste(validation$missing_separator, collapse = ", ")
)
)
}

if (length(validation$too_short) > 0) {
msgs[[length(msgs) + 1]] <- div(
class = "alert alert-warning",
tags$b("Warning: "),
paste0(
"The following sample names do not contain enough tokens to split factor and replicate: ",
paste(validation$too_short, collapse = ", ")
)
)
}

if (length(validation$empty_replicate) > 0) {
msgs[[length(msgs) + 1]] <- div(
class = "alert alert-warning",
tags$b("Warning: "),
paste0(
"The following sample names have an empty replicate field: ",
paste(validation$empty_replicate, collapse = ", ")
)
)
}

if (length(msgs) == 0) {
return(div(class = "alert alert-success", "Sample names passed the basic validation rules."))
}

tagList(msgs)
})

# ---------------------------------------------------------------------------
# output$rnaseq_colname_editor: renders one textInput per non-sample column.
# Hidden once metadata has been saved (to prevent editing after committing).
# Pre-fills each input with the current value from rnaseq_colname_map_rv.
# ---------------------------------------------------------------------------
output$rnaseq_colname_editor <- renderUI({
req(rnaseq_metadata_rv())

if (!is.null(rnaseq_metadata_saved_rv())) return(NULL)

cols <- setdiff(colnames(rnaseq_metadata_rv()), "sample")
if (length(cols) == 0) return(NULL)

current_map <- rnaseq_colname_map_rv()
if (is.null(current_map)) {
current_map <- stats::setNames(colnames(rnaseq_metadata_rv()), colnames(rnaseq_metadata_rv()))
}

tagList(
tags$hr(),
tags$b("Rename columns"),
tags$p(
style = "font-size: 0.85em; color: #666; margin-bottom: 6px;",
"Edit the column names displayed in the table, exported files, and DESeq2 design."
),
lapply(cols, function(col) {
textInput(
inputId = ns(paste0("rnaseq_colname_", col)),
label = col,
value = current_map[[col]] %||% col,
width = "100%"
)
})
)
})

# ---------------------------------------------------------------------------
# output$rnaseq_design_panel: DESeq2 configuration section.
# Visible only after metadata has been saved.
# Provides: design formula input, count filter threshold, build button.
# ---------------------------------------------------------------------------
output$rnaseq_design_panel <- renderUI({
req(rnaseq_metadata_saved_rv())

tagList(
tags$hr(),
step_header(s_design, "Design & build the model"),
textInput(
ns("rnaseq_design_formula"),
lbl("Design formula", "The model DESeq2 fits - it decides what every log2FC and p-value mean. ~ A + B is additive (use it to control for a batch). ~ A * B adds the interaction, and the main effects then apply only at the other factor's control level. Term ORDER does not matter here: you pick the comparison yourself in Contrast below, and the count filter is sized on your smallest group whichever way round you type the formula."),
value = default_rnaseq_design(),
width = "100%"
),
uiOutput(ns("rnaseq_interaction_ui")),
uiOutput(ns("rnaseq_reference_ui")),
numericInput(
ns("rnaseq_min_count_filter"),
lbl("Minimum count threshold for pre-filtering", "Drops genes too lowly expressed to test, before the model is fitted. Not cosmetic: it changes which genes exist, which changes DESeq2's independent filtering, which changes every padj. 10 is typical for raw counts."),
value = 10,
min = 1,
step = 1
),
helpText("Genes are kept if they reach the selected count threshold in at least as many samples as your SMALLEST experimental group — the smallest across every variable in the design, so a gene switched on in only one condition is never filtered out."),
actionButton(
ns("build_rnaseq_dds"),
"Create DESeq2 object",
class = "btn-primary btn-sm",
width = "100%"
),
br(), br(),
uiOutput(ns("rnaseq_dds_status"))
)
})

# ---------------------------------------------------------------------------
# output$rnaseq_interaction_ui: opt-in switch for the interaction model.
# Shown only when an interaction is actually fittable. Unchecked by default so
# the additive model — and therefore the meaning of every existing contrast —
# is preserved unless the user explicitly asks for the interaction.
# ---------------------------------------------------------------------------
output$rnaseq_interaction_ui <- renderUI({
  if (!isTRUE(rnaseq_interaction_ok())) return(NULL)
  terms <- rnaseq_design_terms()
  tagList(
    checkboxInput(
      ns("rnaseq_test_interaction"),
      sprintf("Test interaction between %s and %s", terms[1], terms[2]),
      value = FALSE
    ),
    helpText(
      HTML(sprintf(
        paste0("Rewrites the formula to <code>~ %1$s * %2$s</code> and unlocks the ",
               "interaction coefficients in the DE tab.<br><b>Warning:</b> under an ",
               "interaction model a simple contrast on <code>%1$s</code> is the effect ",
               "at the <i>control level</i> of <code>%2$s</code> only, not the effect ",
               "averaged over %2$s."),
        terms[1], terms[2]))
    )
  )
})

# ---------------------------------------------------------------------------
# Rewrite the formula when the checkbox is toggled. updateTextInput (rather
# than making default_rnaseq_design() depend on the checkbox) keeps the design
# box editable by hand: a manual formula survives until the box is toggled.
# ---------------------------------------------------------------------------
observeEvent(input$rnaseq_test_interaction, {
  terms <- rnaseq_design_terms()
  if (length(terms) != 2) return(NULL)
  op <- if (isTRUE(input$rnaseq_test_interaction)) " * " else " + "
  updateTextInput(
    session,
    "rnaseq_design_formula",
    value = paste("~", paste(terms, collapse = op))
  )
}, ignoreInit = TRUE)

# ---------------------------------------------------------------------------
# output$rnaseq_metadata_ready: logical flag used by conditionalPanel.
# TRUE once a metadata table exists (whether saved or not).
# ---------------------------------------------------------------------------
output$rnaseq_metadata_ready <- reactive({
!is.null(rnaseq_metadata_rv())
})
outputOptions(output, "rnaseq_metadata_ready", suspendWhenHidden = FALSE)

# ---------------------------------------------------------------------------
# output$rnaseq_metadata_ui: switches between the DT editor and a
# locked info banner once the table has been saved.
# ---------------------------------------------------------------------------
output$rnaseq_metadata_ui <- renderUI({
req(rnaseq_metadata_rv())

if (!is.null(rnaseq_metadata_saved_rv())) {
return(div(class = "alert alert-info", "Metadata table saved for analysis."))
}

tagList(
h4("📋 Metadata editor"),
DT::DTOutput(ns("rnaseq_metadata_table"))
)
})

# ---------------------------------------------------------------------------
# output$rnaseq_metadata_table: interactive DT table.
# All cells are editable except column 0 (sample names, index 0 in JS).
# server = FALSE so that cell edits are returned to R immediately.
# ---------------------------------------------------------------------------
output$rnaseq_metadata_table <- DT::renderDT({
req(rnaseq_metadata_for_analysis())

DT::datatable(
rnaseq_metadata_for_analysis(),
rownames = FALSE,
editable = list(
target = "cell",
disable = list(columns = c(0))  # protect the sample column
),
selection = "none",
options = list(
pageLength = 10,
scrollX = TRUE,
dom = "tip"
)
)
}, server = FALSE)

# ---------------------------------------------------------------------------
# observeEvent: rnaseq_metadata_table_cell_edit
# Updates the working metadata reactive when a cell is edited in the DT.
# Column index 0 (sample) is protected and triggers a warning if edited.
# DT uses 0-based column indices; R data frame uses 1-based, hence col + 1.
# ---------------------------------------------------------------------------
observeEvent(input$rnaseq_metadata_table_cell_edit, {
req(rnaseq_metadata_rv())

info <- input$rnaseq_metadata_table_cell_edit
meta <- rnaseq_metadata_rv()

row   <- info$row
col   <- info$col
value <- info$value

if (col == 0) {
showNotification(
"The sample column is protected and cannot be edited.",
type = "warning",
duration = 4
)
return()
}

meta[row, col + 1] <- value
rnaseq_metadata_rv(meta)
}, ignoreNULL = TRUE)

# ---------------------------------------------------------------------------
# Download handlers: CSV and Excel export of the current metadata table.
# File names include a timestamp to avoid overwriting previous exports.
# ---------------------------------------------------------------------------
output$download_rnaseq_metadata_csv <- downloadHandler(
filename = function() {
paste0("rnaseq_metadata_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
},
content = function(file) {
req(rnaseq_metadata_for_analysis())
write.csv(rnaseq_metadata_for_analysis(), file, row.names = FALSE)
}
)

output$download_rnaseq_metadata_xlsx <- downloadHandler(
filename = function() {
paste0("rnaseq_metadata_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".xlsx")
},
content = function(file) {
req(rnaseq_metadata_for_analysis())
openxlsx::write.xlsx(rnaseq_metadata_for_analysis(), file = file, overwrite = TRUE)
}
)

# ---------------------------------------------------------------------------
# observeEvent: build_rnaseq_dds
# Triggered by the "Create DESeq2 object" button.
# Full pipeline:
#   1. Safety checks (sample column present, all samples found in counts).
#   2. Align metadata row order to match counts column order.
#   3. Parse and validate the design formula.
#   4. Build DESeqDataSetFromMatrix (counts rounded to integers).
#   5. Pre-filter low-expression genes:
#        keep if count >= min_count in at least min_group_size samples,
#        where min_group_size = size of smallest group in first design term.
#   6. Fit the full DESeq2 model: estimate size factors, dispersions, GLM.
# On success, rnaseq_dds_rv is updated with the fitted object.
# ---------------------------------------------------------------------------
observeEvent(input$build_rnaseq_dds, {
req(import_result_rv())
req(rnaseq_metadata_saved_rv())
req(input$rnaseq_design_formula)

count_df <- if (!is.null(import_result_rv()) && import_result_rv()$type == "raw") {
  import_result_rv()$counts
} else {
  NULL
}
meta     <- rnaseq_metadata_saved_rv()

# Remove technical DESeq2 columns carried over from previous builds.
meta <- meta[, !colnames(meta) %in% c("sizeFactor", "replaceable", "basepairs"), drop = FALSE]
# --- Safety checks ----------------------------------------------------------
if (!"sample" %in% colnames(meta)) {
showNotification("The saved metadata does not contain the required sample column.", type = "error", duration = 8)
return()
}

# Raw counts only: validate + align here. For Salmon, the equivalent
# validation and row alignment happen in the tximport branch below, where the
# correct count matrix (import_res$txi$counts) is available. count_df is NULL
# for Salmon, so running this block unconditionally aborts every Salmon build.
if (import_result_rv()$type == "raw") {
  if (!setequal(meta$sample, colnames(count_df))) {
    showNotification(
      "Counts matrix columns and metadata sample names do not match exactly.",
      type = "error",
      duration = 8
    )
    return()
  }
  meta <- meta[match(colnames(count_df), meta$sample), , drop = FALSE]
  rownames(meta) <- meta$sample
}

# --- Parse design formula ---------------------------------------------------
design_formula <- tryCatch(as.formula(input$rnaseq_design_formula), error = function(e) NULL)
if (is.null(design_formula)) {
showNotification("Invalid design formula.", type = "error", duration = 8)
return()
}

# --- Apply selected reference levels ---------------------------------------
# All non-sample metadata columns are converted to factors before DESeq2
# object construction. Selected reference levels are applied via relevel().
meta_cols <- setdiff(colnames(meta), c("sample", "replicate"))

for (col in meta_cols) {
  vals <- unique(as.character(meta[[col]]))
  vals <- vals[!is.na(vals) & nzchar(vals)]
  
  if (length(vals) >= 2) {
    meta[[col]] <- as.factor(as.character(meta[[col]]))
    
    ref_val <- input[[paste0("ref_level_", col)]]
    if (!is.null(ref_val) && ref_val %in% levels(meta[[col]])) {
      meta[[col]] <- stats::relevel(meta[[col]], ref = ref_val)
    }
  }
}

design_terms <- all.vars(design_formula)
if (length(design_terms) == 0) {
showNotification("The design formula must contain at least one metadata variable.", type = "error", duration = 8)
return()
}

missing_terms <- setdiff(design_terms, colnames(meta))
if (length(missing_terms) > 0) {
showNotification(
paste0("The following design variables are not present in the saved metadata: ", paste(missing_terms, collapse = ", ")),
type = "error",
duration = 10
)
return()
}

# --- Build DESeqDataSet ------------------------------------------------------
# [SALMON-SRV-4] Branch on import type to build the DESeqDataSet.
# Raw:    DESeqDataSetFromMatrix (unchanged behaviour)
# Salmon: DESeqDataSetFromTximport (uses the txi list from import_result_rv)
import_res <- import_result_rv()

if (is.null(import_res)) {
  showNotification("No counts data loaded. Please upload a file first.",
                   type = "error", duration = 8)
  return()
}

if (import_res$type == "raw") {
  dds <- tryCatch({
    DESeq2::DESeqDataSetFromMatrix(
      countData = { m <- round(as.matrix(count_df)); m[m < 0] <- 0L; m },
      colData   = meta,
      design    = design_formula
    )
  }, error = function(e) {
    showNotification(
      paste0("Failed to create DESeq2 object: ", conditionMessage(e)),
      type = "error", duration = 10
    )
    NULL
  })
  
} else {
  # Salmon branch: DESeqDataSetFromTximport does not need round() —
  # tximport already handles count scaling internally.
  # Reference: https://bioconductor.org/packages/release/bioc/vignettes/
  #            tximport/inst/doc/tximport.html#DESeq2
  if (!setequal(meta$sample, colnames(import_res$txi$counts))) {
    showNotification(
      paste0("Sample names in Salmon output do not match the metadata. ",
             "Salmon: ", paste(sort(colnames(import_res$txi$counts)), collapse = ", "),
             ". Metadata: ", paste(sort(meta$sample), collapse = ", ")),
      type = "error", duration = 12
    )
    return()
  }
  meta <- meta[match(colnames(import_res$txi$counts), meta$sample), , drop = FALSE]
  rownames(meta) <- meta$sample
  
  dds <- tryCatch({
    DESeq2::DESeqDataSetFromTximport(
      txi     = import_res$txi,
      colData = meta,
      design  = design_formula
    )
  }, error = function(e) {
    showNotification(
      paste0("Failed to create DESeq2 object from Salmon data: ", conditionMessage(e)),
      type = "error", duration = 10
    )
    NULL
  })
}

if (is.null(dds)) return()

if (!design_terms[1] %in% colnames(SummarizedExperiment::colData(dds))) {
showNotification("The first design variable could not be identified for count filtering.", type = "error", duration = 8)
return()
}

# --- Pre-filter genes -------------------------------------------------------
# Strategy: a gene is kept if it reaches min_count reads in at least
# as many samples as the size of the smallest experimental group.
# --- Fit DESeq2 model -------------------------------------------------------
# DESeq() runs: estimateSizeFactors → estimateDispersions → nbinomWaldTest
# Sized on the smallest group across EVERY design variable, not just the first.
# DESeq2's rule is "at least as many samples as the smallest group of the
# variable being TESTED", because a gene switched on in one condition and silent
# in the rest is real biology and must survive the filter. The contrast is not
# known here -- the dds is built before the Contrast selectors are populated from
# it -- so the minimum over all design terms is the conservative stand-in: it can
# never demand a gene appear in more samples than the smallest condition
# contains. Taking design_terms[1] instead made the threshold depend on the
# ORDER the formula was typed in, and a first factor with larger groups deleted
# exactly those condition-specific genes before the model ever saw them.
# max(1L, ...) guards a factor carrying an unused level, whose table() entry is 0.
# grouping_var still names the variable the size came from, so the record written
# into metadata(dds)$cheski_prefilter below stays accurate.
group_sizes    <- vapply(design_terms, function(v)
                    min(table(SummarizedExperiment::colData(dds)[[v]])), integer(1))
grouping_var   <- names(group_sizes)[which.min(group_sizes)]
min_group_size <- max(1L, min(group_sizes))
min_count      <- as.numeric(input$rnaseq_min_count_filter)

keep <- rowSums(DESeq2::counts(dds) >= min_count) >= min_group_size
if (sum(keep) < 2) {
  showNotification(
    paste0("Only ", sum(keep), " gene(s) pass the count filter (need \u2265 2). ",
           "Lower the minimum count threshold."),
    type = "error", duration = 10)
  return()
}
dds <- dds[keep, ]

# ── Record the pre-filter so it travels with the object ─────────────────────
# min_count decides which genes exist, therefore which baseMean quantile
# DESeq2's independent filtering picks, therefore EVERY padj in the results
# table. It was read from the input and discarded, so rerun_DE.R could only
# emit "apply the same pre-filtering you used in the app".
# S4Vectors::metadata() is used because it survives both the dds[keep, ]
# subset above and DESeq() below, so downstream code can ask the fitted object
# instead of re-reading a live input that may have changed since.
S4Vectors::metadata(dds)$cheski_prefilter <- list(
  min_count      = min_count,
  min_group_size = as.integer(min_group_size),
  grouping_var   = grouping_var,
  n_genes_before = length(keep),
  n_genes_after  = sum(keep)
)

dds <- tryCatch(
  # Wrapped so DESeq()'s own narration (estimating size factors -> gene-wise
  # dispersion estimates -> fitting model and testing) is shown live instead of
  # only reaching the console. See with_progress_messages() in log_utils.R.
  with_progress_messages(DESeq2::DESeq(dds), "Fitting the DESeq2 model"),
  error = function(e) {
    showNotification(paste0("DESeq() failed: ", conditionMessage(e)),
                     type = "error", duration = 12)
    NULL
  }
)
if (is.null(dds)) return()
rnaseq_dds_rv(dds)

# Fold the metadata block away now that the object exists. Only on SUCCESS --
# every early return above leaves it open, which is exactly when the user still
# needs the inputs that produced the error.
#
# Setting .open on the element rather than re-rendering the panel: a renderUI
# would rebuild every input inside and lose what the user typed. This just
# closes the disclosure; the values are untouched and one click brings them
# back. The optional-chaining-free guard keeps it silent if the id ever moves.
shinyjs::runjs(sprintf(
  "var d = document.getElementById('%s'); if (d) { d.open = false; }",
  session$ns("meta_panel")))

showNotification(
paste0("DESeq2 object created successfully. ", sum(keep),
       " genes retained after filtering. The metadata section has been collapsed — ",
       "click its title to reopen it."),
type = "message",
duration = 6
)
})

# ---------------------------------------------------------------------------
# output$rnaseq_dds_status: green badge shown once the dds is ready.
# ---------------------------------------------------------------------------
output$rnaseq_dds_status <- renderUI({
if (is.null(rnaseq_dds_rv())) return(NULL)
div(class = "alert alert-success", "DESeq2 object ready.")
})

# ---------------------------------------------------------------------------
# Return value: named list of four reactives consumed by downstream modules.
# ---------------------------------------------------------------------------
list(
dds                  = reactive(rnaseq_dds_rv()),
metadata_saved       = reactive(rnaseq_metadata_saved_rv()),
metadata_for_analysis = reactive(rnaseq_metadata_for_analysis()),
counts               = reactive(mat_data_rnaseq()),
reset                = function() {
  import_result_rv(NULL); rnaseq_metadata_rv(NULL)
  rnaseq_metadata_saved_rv(NULL); rnaseq_dds_rv(NULL)
  rnaseq_validation_rv(NULL); rnaseq_colname_map_rv(NULL)
}
)
})
}
