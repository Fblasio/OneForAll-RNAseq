library(shiny)
library(colourpicker)
library(ComplexHeatmap)
library(InteractiveComplexHeatmap)
library(circlize)
library(readxl)
#library(org.At.tair.db)
library(AnnotationDbi)
library(clusterProfiler)
library(DOSE)
library(dplyr)
library(DT)
library(ggplot2)
library(httr)
library(jsonlite)
library(digest)
library(shinyjs)
library(cluster)
library(shinycssloaders) 
library(pathview) 
library(base64enc)
library(httr)
library(jsonlite)
library(DESeq2)
library(ashr)
library(openxlsx)
library(ggrepel)
library(later)
library(tibble)
library(scales)
library(grDevices)
library(ggVennDiagram)
library(sf)
library(sp)
library(ggvenn)
library(stringr)

MAX_GROUPS        <- 12
MAX_ANNOTS        <- 4
MAX_GENES_MEMORY  <- 5000

`%||%` <- function(a, b) {
  if (is.null(a) || length(a) == 0) return(b)
  if (is.character(a) && length(a) == 1L && !nzchar(a)) return(b)
  a
}

# ── Inline "?" help ─────────────────────────────────────────────────────────
# Defined HERE, above the source() calls at the bottom of this file, so both
# app.r and every module UI function can see it (module files are sourced into
# .GlobalEnv, so anything defined only inside app.r would be invisible to them).
#
# Scope on purpose: these are attached ONLY to controls where a wrong setting
# changes the RESULT silently -- no error, no empty panel, just different
# numbers -- and where the label alone cannot tell you which option to pick.
# Cosmetic controls (colours, font sizes, plot dimensions) and controls that
# already carry an adequate helpText() are deliberately left alone, so a "?"
# in this app always means "this one can change your conclusions".
#
# Pure CSS, no new package: the text rides on a data-tip attribute and is
# revealed by :hover/:focus (see the .cheski-help rules in app.r's tags$head).
# tabindex makes it keyboard-reachable, which is why :focus is styled too.
help_tip <- function(text) {
  shiny::tags$span(class = "cheski-help", `data-tip` = text,
                   tabindex = "0", role = "note", "?")
}

# Wraps a control label so the "?" sits next to it. Shiny renders a tagList
# passed as `label` inside the control's own <label> element, so this works for
# selectInput / numericInput / checkboxInput / radioButtons / sliderInput alike.
lbl <- function(label, help) {
  shiny::tagList(label, help_tip(help))
}

# ── Numbered workflow step ──────────────────────────────────────────────────
# The sidebar of the Heatmap and RNAseq features is a SEQUENCE -- species,
# upload, describe the experiment, run -- but nothing said so, and the headings
# read as an unordered pile of settings. This puts a number next to each one.
#
# Defined here, above the source() calls at the bottom of this file, for the
# same reason help_tip() and collapsible_panel() are: module files are sourced
# into .GlobalEnv, so anything defined only inside app.r would be invisible to
# them, and rnaseq_dds_module.R / rnaseq_de_module.R both need this.
#
# n        : the badge text. A string, not an integer, so sub-steps like "4a"
#            work without a second function.
# optional : pale badge plus an "optional" tag, so numbering a step cannot
#            imply it is required. Steps that genuinely gate the next one are
#            left solid.
# help     : optional help_tip() text, same convention as lbl().
step_header <- function(n, title, help = NULL, optional = FALSE) {
  shiny::tags$div(
    style = "display:flex; align-items:center; gap:8px; margin:14px 0 6px;",
    shiny::tags$span(
      style = paste0(
        "flex:0 0 auto; width:22px; height:22px; border-radius:50%; ",
        "background:", if (optional) "#dfe4ee" else "#1F3864", "; ",
        "color:", if (optional) "#1F3864" else "#ffffff", "; ",
        "font-size:12px; font-weight:700; line-height:22px; text-align:center;"),
      n),
    shiny::tags$span(style = "font-size:15px; font-weight:bold; color:#1F3864;",
                     title),
    if (isTRUE(optional))
      shiny::tags$small(style = "color:#888; font-weight:normal;", "optional"),
    if (!is.null(help)) help_tip(help)
  )
}

# ── Foldable panel ──────────────────────────────────────────────────────────
# A native <details>/<summary> disclosure block: click the title to collapse or
# expand. No JS accordion, no new package, and it keeps working mid-stream.
#
# This is the same construction enrichPlotStyleUI() already uses for the GO and
# KEGG "Plot style" blocks, hoisted here so any module can reach it -- module
# files are sourced into .GlobalEnv, so anything defined only inside app.r would
# be invisible to them (same reason help_tip() lives here).
#
# htmltools::tag() rather than tags$details(): tag() is the generic constructor
# present in every htmltools version, unlike the tags$xxx() shorthands for the
# newer HTML5 elements.
#
# title   : the always-visible summary line
# open    : TRUE renders it expanded, which is the default so that adding this
#           wrapper around an existing block changes nothing on first load.
#           `open` is an HTML boolean attribute, hence the NA rather than TRUE.
# heading : NULL (default) uses the body size, for a control block inside a tab.
#           "h3" / "h4" size the summary like that heading, for folding a whole
#           SECTION whose title sits alongside real headings -- otherwise the
#           folded section's title renders smaller than its unfolded siblings
#           and the page hierarchy reads wrong. Bootstrap 3 on a 14px body:
#           h3 = 24px (1.7em), h4 = 18px (1.29em).
# id      : optional DOM id, so the server can open or close the panel later
#           with shinyjs::runjs(). Namespace it with ns() from inside a module.
collapsible_panel <- function(title, ..., open = TRUE, heading = NULL, id = NULL) {
  size <- switch(heading %||% "", h3 = " font-size:1.7em;", h4 = " font-size:1.29em;", "")
  args <- list(
    style = paste0("background:#f8f9fa; border:1px solid #dee2e6; border-radius:4px; ",
                   "padding:6px 10px 10px 10px; margin:8px 0;"),
    htmltools::tag("summary", list(
      style = paste0("cursor:pointer; font-weight:bold; color:#1F3864; padding:4px 0;", size),
      title
    )),
    shiny::div(style = "margin-top:8px;", ...)
  )
  if (!is.null(id))   args$id   <- id
  if (isTRUE(open))   args$open <- NA
  htmltools::tag("details", args)
}

# ── Local-launch auto-close support ───────────────────────────────────────
# Tracks how many browser sessions are currently connected to this single
# local R process. Only active when CHESKI_AUTOCLOSE=1 (set by run_app.R),
# so running app.r directly in RStudio for development is unaffected.
.cheski_launch_state <- new.env(parent = emptyenv())
.cheski_launch_state$active_sessions <- 0L

# Call once, at the very top of server(). Increments the active-session
# counter on connect; on disconnect, waits `grace_seconds` to see whether
# a new session reconnects (e.g. a page refresh) before treating it as a
# real close and calling stopApp().
register_autoclose_on_browser_close <- function(session, grace_seconds = 3) {
  if (!identical(Sys.getenv("CHESKI_AUTOCLOSE"), "1")) return(invisible(NULL))
  
  .cheski_launch_state$active_sessions <- .cheski_launch_state$active_sessions + 1L
  
  session$onSessionEnded(function() {
    .cheski_launch_state$active_sessions <- .cheski_launch_state$active_sessions - 1L
    later::later(function() {
      if (.cheski_launch_state$active_sessions <= 0L) {
        message("One for All: no browser session reconnected within ",
                grace_seconds, "s \u2014 shutting down.")
        shiny::stopApp()
      }
    }, delay = grace_seconds)
  })
  invisible(NULL)
}

source("utils/enrichment_functions.R")
source("utils/gsea_functions.R") 
source("utils/color_utils.R")
source("utils/normalization_utils.R")
source("utils/memory_utils.R")
source("utils/import_utils.R")
source("utils/rnaseq_metadata.R")
source("utils/species_config.R")
source("utils/consistency_utils.R")  
source("utils/validation_utils.R")   
source("utils/trust_utils.R")         
source("utils/log_utils.R")           
source("modules/go_module.R")
source("modules/gsea_module.R")
source("modules/kegg_module.R")
source("modules/selection_history_module.R")
source("modules/string_module.R")
source("modules/tair_module.R")
source("modules/gene_detail_module.R")
source("modules/optimal_k_module.R")
source("modules/rnaseq_viz_module.R")
source("modules/rnaseq_eda_module.R")
source("modules/rnaseq_dds_module.R")
source("modules/rnaseq_de_module.R")
source("modules/rnaseq_volcano_module.R")
source("modules/venn_module_v2.r")
source("modules/venn_manual_module.r")
source("modules/repro_module.R")    
source("modules/de_comparison_log.R")
#source("modules/de_cross_comparison_module.R")
source("modules/demo_module.R")       
source("modules/log_module.R")
source("modules/go_section_module.R")
Sys.setenv(LANGUAGE = "en")
Sys.setenv(LANGUAGE = "en") 