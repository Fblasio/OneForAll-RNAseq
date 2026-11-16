selectionHistoryUI <- function(id, enable_load = TRUE) {
  ns <- NS(id)
  tagList(
    br(),
    fluidRow(
      column(8, textInput(ns("sel_name"), label = NULL,
                          placeholder = "Name this selection...")),
      column(4, actionButton(ns("save_sel"), "💾 Save",
                             class = "btn-primary btn-sm", width = "100%"))
    ),
    br(),
    helpText("Click a row then use Load/Delete."),
    DTOutput(ns("history_table")),
    br(),
    if (isTRUE(enable_load)) {
      fluidRow(
        column(6, actionButton(ns("load_sel"), "📂 Load selected",
                               class = "btn-info btn-sm", width = "100%")),
        column(6, actionButton(ns("delete_sel"), "🗑️ Delete selected",
                               class = "btn-danger btn-sm", width = "100%"))
      )
    } else {
      fluidRow(
        column(12, actionButton(ns("delete_sel"), "🗑️ Delete selected",
                                class = "btn-danger btn-sm", width = "100%"))
      )
    },
    br(),
    downloadButton(ns("dl_history"), "⬇ Export all as CSV",
                   class = "btn-sm btn-success")
  )
}

selectionHistoryServer <- function(id, current_genes) {
  moduleServer(id, function(input, output, session) {
    
    history      <- reactiveVal(list())
    loaded_genes <- reactiveVal(NULL)
    
    observeEvent(input$save_sel, {
      genes <- current_genes()
      req(genes)
      name <- if (nzchar(trimws(input$sel_name %||% ""))) {
        trimws(input$sel_name)
      } else {
        paste0("Selection_", length(history()) + 1)
      }
      history(c(history(), list(list(
        name  = name,
        genes = genes,
        n     = length(genes),
        time  = format(Sys.time(), "%H:%M:%S")
      ))))
      updateTextInput(session, "sel_name", value = "")
      showNotification(paste0("✅ '", name, "' saved (", length(genes), " genes)"),
                       type = "message", duration = 3)
    })
    
    output$history_table <- renderDT({
      h <- history()
      if (length(h) == 0)
        return(datatable(data.frame(Message = "No selections saved yet."),
                         rownames = FALSE, options = list(dom = "t")))
      df <- data.frame(
        `#`    = seq_along(h),
        Name   = sapply(h, `[[`, "name"),
        Genes  = sapply(h, `[[`, "n"),
        Time   = sapply(h, `[[`, "time"),
        check.names = FALSE
      )
      datatable(df, rownames = FALSE, selection = "single",
                options = list(dom = "tp", pageLength = 10))
    })
    
    observeEvent(input$load_sel, {
      sel <- input$history_table_rows_selected
      req(sel)
      h <- history()
      loaded_genes(h[[sel]]$genes)
      showNotification(
        paste0("📂 Loaded '", h[[sel]]$name, "' — ", h[[sel]]$n, " genes"),
        type = "message", duration = 3)
    })
    
    observeEvent(input$delete_sel, {
      sel <- input$history_table_rows_selected
      req(sel)
      h    <- history()
      name <- h[[sel]]$name
      history(h[-sel])
      showNotification(paste0("🗑️ Deleted '", name, "'"),
                       type = "warning", duration = 3)
    })
    
    output$dl_history <- downloadHandler(
      filename = function() paste0("selection_history_",
                                   format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
      content  = function(file) {
        h <- history()
        if (length(h) == 0) {
          write.csv(data.frame(Message = "No selections."), file, row.names = FALSE)
          return()
        }
        df <- do.call(rbind, lapply(h, function(e)
          data.frame(Name = e$name, Time = e$time,
                     GeneID = e$genes, stringsAsFactors = FALSE)
        ))
        write.csv(df, file, row.names = FALSE)
      }
    )
    
    return(list(
      loaded_genes = loaded_genes,
      history      = history
    ))
  })
}
