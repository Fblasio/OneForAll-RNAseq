# =============================================================================
# log_module.R
# -----------------------------------------------------------------------------
# Tiny module exposing ONE download button for the session log built by
# utils/log_utils.R. Put the UI anywhere in app.R's sidebar (e.g. next to
# "Clear enrichment cache") and call the server once from the main server.
#
# Usage in app.R:
#   UI:     appLogUI("applog")
#   Server: appLogServer("applog")
# =============================================================================

appLogUI <- function(id) {
  ns <- NS(id)
  tagList(
    downloadButton(ns("dl_log"), "\u2b07 Download session log (.txt)",
                   class = "btn-sm btn-default", width = "100%"),
    tags$p(style = "font-size:0.8em;color:#888;margin-top:4px;",
          "Records uploads, validation warnings/errors, and empty-result reasons from this session.")
  )
}

appLogServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    output$dl_log <- downloadHandler(
      filename = function() {
        paste0("OneForAll_session_log_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt")
      },
      content = function(file) {
        lines <- get_session_log(session)
        if (length(lines) == 0) {
          lines <- c(
            "One for All session log",
            "No events were recorded in this session.",
            "(Events are logged by validate_count_matrix(), enrichment empty-result",
            " checks, and other checkpoints — see utils/log_utils.R log_event().)"
          )
        } else {
          lines <- c(
            "One for All session log",
            paste0("Downloaded: ", format(Sys.time())),
            strrep("-", 60),
            lines
          )
        }
        writeLines(lines, file)
      }
    )
  })
}
