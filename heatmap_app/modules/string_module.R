stringUI <- function(id) {
  ns <- NS(id)
  tagList(
    br(),
    fluidRow(
      column(4,
             numericInput(ns("score"), "Min interaction score = confidence.",
                          value = 400, min = 150, max = 999, step = 50),
             helpText("150=low | 400=medium | 700=high | 900=highest"),
             numericInput(ns("add_nodes"), "Extra proteins from STRING",
                          value = 0, min = 0, max = 50, step = 5),
             helpText("0 = only your input genes"),
             actionButton(ns("run"), "▶ Query STRING",
                          class = "btn-primary btn-sm", width = "100%"),
             br(), br(),
             downloadButton(ns("dl_table"), "⬇ Interactions CSV",
                            class = "btn-sm btn-success")
      ),
      column(8,
             uiOutput(ns("string_link")),
             br(),
             DTOutput(ns("interaction_table"))
      )
    )
  )
}

stringServer <- function(id, selected_genes) {
  moduleServer(id, function(input, output, session) {
    
    result <- eventReactive(input$run, {
      genes <- selected_genes()
      req(genes)
      
      withProgress(message = "Querying STRING...", value = 0.2, {
        
        taxon <- UNIPROT_TAXID
        base_url <- "https://string-db.org/api/json/"
        
        genes_clean <- unique(canon_ids(genes))
        genes_clean <- genes_clean[!is.na(genes_clean) & nzchar(genes_clean)]
        
        if (length(genes_clean) == 0) {
          return(list(error = "❌ No valid gene IDs available for STRING query."))
        }
        
        # ── Translate gene IDs -> protein names for STRING ──────────────────
        # STRING resolves tomato on protein name (e.g. A0A3Q7E704), not on the
        # Solyc gene ID, so the raw IDs map poorly or not at all. When a side
        # map is available (STRING_MAP, see species_config.R) query with the
        # protein names instead and keep a reverse lookup so the network labels
        # still show the gene the user actually selected.
        # No map -> query_ids stays as the gene IDs = previous behaviour.
        query_ids <- genes_clean
        prot2gene <- NULL
        n_unmapped <- 0L
        if (exists("STRING_MAP") && !is.null(STRING_MAP)) {
          sel <- STRING_MAP[STRING_MAP[[1]] %in% genes_clean, , drop = FALSE]
          if (nrow(sel) > 0) {
            query_ids  <- unique(sel[[2]])
            prot2gene  <- stats::setNames(as.character(sel[[1]]),
                                          as.character(sel[[2]]))
            n_unmapped <- length(setdiff(genes_clean, unique(sel[[1]])))
          }
        }
        
        setProgress(0.3, detail = "Mapping gene IDs...")
        
        map_resp <- tryCatch(
          httr::POST(
            paste0(base_url, "get_string_ids"),
            body = list(
              identifiers = paste(query_ids, collapse = "\r"),
              species = taxon,
              limit = 1,
              caller_identity = "shiny_heatmap_app"
            ),
            encode = "form"
          ),
          error = function(e) NULL
        )
        
        if (is.null(map_resp) || httr::status_code(map_resp) != 200) {
          return(list(error = "❌ Cannot connect to STRING. Check internet connection."))
        }
        
        mapped <- tryCatch(
          jsonlite::fromJSON(httr::content(map_resp, "text", encoding = "UTF-8")),
          error = function(e) NULL
        )
        
        if (is.null(mapped) || !is.data.frame(mapped) || nrow(mapped) < 1) {
          return(list(error = "❌ No STRING identifiers found for these genes."))
        }
        
        if (!all(c("stringId", "queryItem") %in% colnames(mapped))) {
          return(list(error = "❌ STRING mapping response is missing expected columns."))
        }
        
        setProgress(0.6, detail = "Fetching interactions...")
        
        net_resp <- tryCatch(
          httr::POST(
            paste0(base_url, "network"),
            body = list(
              identifiers = paste(mapped$stringId, collapse = "\r"),
              species = taxon,
              required_score = input$score,
              add_nodes = input$add_nodes,
              caller_identity = "shiny_heatmap_app"
            ),
            encode = "form"
          ),
          error = function(e) NULL
        )
        
        if (is.null(net_resp) || httr::status_code(net_resp) != 200) {
          return(list(error = "❌ Failed to fetch STRING network."))
        }
        
        interactions <- tryCatch(
          jsonlite::fromJSON(httr::content(net_resp, "text", encoding = "UTF-8")),
          error = function(e) NULL
        )
        
        if (is.null(interactions) || !is.data.frame(interactions) || nrow(interactions) < 1) {
          return(list(error = paste0("No interactions found above score ", input$score, ".")))
        }
        
        if (!all(c("stringId_A", "stringId_B", "preferredName_A", "preferredName_B", "score") %in% colnames(interactions))) {
          return(list(error = "❌ STRING network response is missing expected columns."))
        }
        
        setProgress(1.0)
        
        str_url <- paste0(
          "https://string-db.org/cgi/network?identifiers=",
          paste(mapped$stringId, collapse = "%0d"),
          "&species=", taxon,
          "&required_score=", input$score
        )
        
        # queryItem holds whatever we sent — protein names when STRING_MAP was
        # used. Translate back so labels read "preferredName (Solyc...)" rather
        # than exposing the protein accession the user never asked about.
        query_item <- as.character(mapped$queryItem)
        if (!is.null(prot2gene)) {
          hit <- query_item %in% names(prot2gene)
          query_item[hit] <- unname(prot2gene[query_item[hit]])
        }
        id_lookup <- stats::setNames(query_item, mapped$stringId)
        
        make_label <- function(str_id, pref_name) {
          gene_id <- id_lookup[[str_id]]
          if (!is.null(gene_id) && !is.na(gene_id) && nzchar(gene_id)) {
            paste0(pref_name, " (", gene_id, ")")
          } else {
            pref_name
          }
        }
        
        idf <- data.frame(
          `Gene A` = mapply(make_label, interactions$stringId_A, interactions$preferredName_A),
          `Gene B` = mapply(make_label, interactions$stringId_B, interactions$preferredName_B),
          `Score` = interactions$score,
          stringsAsFactors = FALSE,
          check.names = FALSE
        )
        
        list(
          interactions = idf,
          url = str_url,
          n_mapped = nrow(mapped),
          n_unmapped = n_unmapped,
          n_input = length(genes_clean),
          error = NULL
        )
      })
    })
    
    output$string_link <- renderUI({
      res <- string_result <- result()
      req(res)
      if (!is.null(res$error))
        return(tags$p(style = "color:red; font-weight:bold;", res$error))
      tagList(
        tags$p(paste0("✅ ", res$n_mapped, " genes mapped to STRING. ",
                      nrow(res$interactions), " interactions found.")),
        if (!is.null(res$n_unmapped) && res$n_unmapped > 0)
          tags$p(style = "color:#8a6d3b;",
                 sprintf("Note: %d of %d selected gene(s) have no protein name in the ID map and could not be sent to STRING.",
                         res$n_unmapped, res$n_input)),
        tags$a(href = res$url, target = "_blank",
               class = "btn btn-info btn-sm",
               "🔗 Open full network on STRING website")
      )
    })
    
    output$interaction_table <- renderDT({
      res <- result()
      req(res, is.null(res$error))
      datatable(res$interactions,
                rownames = FALSE,
                options  = list(pageLength = 15, scrollX = TRUE,
                                order = list(list(2, "desc"))))
    })
    
    output$dl_table <- downloadHandler(
      filename = function() paste0("STRING_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
      content  = function(file) {
        res <- result()
        req(res, is.null(res$error))
        write.csv(res$interactions, file, row.names = FALSE)
      }
    )
  })
}
