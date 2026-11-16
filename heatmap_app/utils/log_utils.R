# =============================================================================
# log_utils.R
# -----------------------------------------------------------------------------
# Session-scoped event log, independent from the reproducibility bundle.
#
# 
#   - manifest.json / session_info.txt (repro bundle) = the PARAMETERS of one
#     analysis (heatmap config, DE contrast, versions) — reproducibility.
#   - this file = a running TIMELINE of what happened this session (uploads,
#     validation warnings, empty-result reasons, errors) — troubleshooting.
#   A user asking "why did this fail" wants the log; a reviewer asking
#     "how do I reproduce this" wants the bundle. Both are one click.
#
# HOW STORAGE WORKS:
#   session$userData is a plain environment that Shiny creates once per
#   session and shares, UNCHANGED, across every module's proxy session. That
#   means log_event() called from inside any module (go_module, kegg_module,
#   rnaseq_de_module, ...) and app_log_download() called from the top-level
#   server both read/write the SAME buffer, with no wiring needed beyond
#   sourcing this file. It is a plain environment, not a reactiveVal, so
#   appending to it does not trigger any invalidation/re-render.
# =============================================================================

# Append one line to the session log. Safe to call with no session (no-ops
# storage, still respects the console debug switch from consistency_utils.R).
#   type: "info" | "warning" | "error" | "debug"
log_event <- function(msg, type = "info",
                      session = shiny::getDefaultReactiveDomain()) {
  line <- sprintf("[%s] %-7s %s",
                  format(Sys.time(), "%Y-%m-%d %H:%M:%S"), toupper(type), msg)

  if (!is.null(session) && !is.null(session$userData)) {
    if (is.null(session$userData$cheski_log)) {
      session$userData$cheski_log <- character(0)
    }
    session$userData$cheski_log <- c(session$userData$cheski_log, line)
  }

  # Reuse the existing debug switch (consistency_utils.R) for console echo.
  if (exists(".cheski_debug_on", mode = "function") && .cheski_debug_on()) {
    message(line)
  }
  invisible(line)
}

# Retrieve the accumulated log for the current (or given) session.
get_session_log <- function(session = shiny::getDefaultReactiveDomain()) {
  if (is.null(session) || is.null(session$userData$cheski_log)) return(character(0))
  session$userData$cheski_log
}

# Convenience wrappers so call sites read naturally.
log_info    <- function(msg, session = shiny::getDefaultReactiveDomain()) log_event(msg, "info",    session)
log_warning <- function(msg, session = shiny::getDefaultReactiveDomain()) log_event(msg, "warning", session)
log_error   <- function(msg, session = shiny::getDefaultReactiveDomain()) log_event(msg, "error",   session)

# --- Show a slow computation's own progress messages on screen ---------------
# DESeq2, clusterProfiler and fgsea all narrate what they are doing with
# message() -- "estimating size factors", "gene-wise dispersion estimates",
# "fitting model and testing". Those lines go to the console, where the user
# never sees them, so a long step looks like a frozen app.
#
# This runs `expr` with a calling handler that intercepts each message and
# writes it into a Shiny progress panel as it arrives. Two things make that
# work, and neither is obvious:
#
#   * Shiny is single-threaded, so anything written to a reactiveVal during a
#     long observer only reaches the browser when the observer RETURNS -- all
#     the messages would appear at once, at the end, which is useless.
#     Progress$set() is the exception: it writes to the websocket immediately.
#   * withCallingHandlers() (not tryCatch) lets the computation carry on after
#     the handler runs. invokeRestart("muffleMessage") then stops the line also
#     going to the console, so the message is moved rather than duplicated.
#
# The computation itself is untouched: this only redirects its narration. Errors
# propagate unchanged, so existing tryCatch() wrappers keep working exactly as
# before.
#
# expr    : the slow call. Lazily evaluated -- forced inside the handler context.
# label   : the progress panel's headline ("Fitting the DESeq2 model").
# session : NULL (no Shiny session, e.g. the test suite) just runs expr.
#
# The number of messages is not fixed -- DESeq() emits 6 on the happy path but
# adds more for a dispersion-fit fallback or an outlier refit -- so the bar is
# capped at 0.95 and never claims to be finished early.
with_progress_messages <- function(expr, label,
                                   session = shiny::getDefaultReactiveDomain(),
                                   log = TRUE) {
  if (is.null(session)) return(expr)
  # A session that cannot make a Progress (a mock session in the test suite)
  # must not take the analysis down with it: fall back to running expr plainly.
  p <- tryCatch(shiny::Progress$new(session, min = 0, max = 1),
                error = function(e) NULL)
  if (is.null(p)) return(expr)
  on.exit(try(p$close(), silent = TRUE), add = TRUE)
  try(p$set(message = label, detail = "starting…", value = 0.05), silent = TRUE)

  n <- 0L
  withCallingHandlers(
    expr,
    message = function(m) {
      txt <- trimws(sub("\n$", "", conditionMessage(m)))
      if (nzchar(txt)) {
        n <<- n + 1L
        try(p$set(detail = txt, value = min(0.95, 0.05 + n * 0.15)), silent = TRUE)
        # Also written to the downloadable session log, so the sequence of
        # steps survives the run instead of vanishing with the progress panel.
        if (isTRUE(log)) log_info(paste0(label, ": ", txt), session)
      }
      invokeRestart("muffleMessage")
    })
}

# --- Optional integration with validation_utils.R ---------------------------
# Drop-in replacement for notify_validation() that ALSO writes to the log, so
# every warning/error shown as a notification is also captured for download.
# Use instead of notify_validation() where you want both.
notify_and_log_validation <- function(v, session = shiny::getDefaultReactiveDomain()) {
  for (w in v$warnings) {
    shiny::showNotification(paste0("\u26a0\ufe0f ", w), type = "warning", duration = 9)
    log_warning(w, session)
  }
  for (e in v$errors) {
    shiny::showNotification(paste0("\u274c ", e), type = "error", duration = 12)
    log_error(e, session)
  }
  invisible(isTRUE(v$ok))
}
