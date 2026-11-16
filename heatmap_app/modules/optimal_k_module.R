# optimal_k_module.R
# Suggests optimal k (Elbow + Silhouette) for dendrogram cutting.
# Renders compact plots inside the sidebar, below the k_clusters slider.

optimalKUI <- function(id) {
  ns <- NS(id)
  uiOutput(ns("panel"))
}

optimalKServer <- function(id, cluster_matrix_rv, linkage_rv = reactive("complete")) {
  moduleServer(id, function(input, output, session) {
    
    result <- reactive({
      mat <- cluster_matrix_rv()
      req(mat)
      req(nrow(mat) >= 4)
      
      k_max <- min(10, nrow(mat) - 1)
      k_min <- 2
      if (k_max < k_min) {
        log_warning(sprintf(
          "k-suggestion preview: only %d gene(s) available, not enough to test multiple k.",
          nrow(mat)))
        return(NULL)
      }
      
      # Score the EXACT tree the main heatmap cuts via cutree() when "Cut
      # dendrogram into k clusters" is checked, by going through the same
      # contract it does: cluster_rows_or_reason() (normalization_utils.R).
      # Re-implementing the tree here with its own `hclust(1 - cor, ...)`
      # would make the suggestion agree with the cut only by coincidence, and
      # silently diverge the moment either side changed its distance or
      # linkage. Routing through the shared contract means both can only ever
      # come from one place.
      # `d` is reused from the helper for the same reason: rebuilding it here
      # would be a second chance to disagree.
      # Named cl_res, not cl: `cl` is already the cutree() loop variable below.
      cl_res <- cluster_rows_or_reason(mat, linkage = linkage_rv())
      if (isFALSE(cl_res$hc)) {
        log_warning(paste("k-suggestion preview:", cl_res$reason))
        return(NULL)
      }
      d  <- cl_res$d
      hc <- cl_res$hc
      
      # Within-cluster dispersion generalizes k-means' WCSS to an arbitrary
      # distance matrix: for each cluster, sum of squared distances (in the
      # SAME correlation-distance space as `d`) from every member to the
      # cluster's medoid (the member minimizing total in-cluster distance).
      .dispersion <- function(cl) {
        d_mat <- as.matrix(d)
        total <- 0
        for (cid in unique(cl)) {
          members <- which(cl == cid)
          if (length(members) <= 1) next
          sub_d      <- d_mat[members, members, drop = FALSE]
          medoid_idx <- which.min(rowSums(sub_d))
          total      <- total + sum(sub_d[medoid_idx, ]^2)
        }
        total
      }
      
      wcss       <- numeric(k_max)
      sil_scores <- numeric(k_max)
      
      withProgress(message = "Computing optimal k...", value = 0, session = session, {
        for (k in 1:k_max) {
          incProgress(1 / k_max, detail = paste("Testing k =", k))
          cl      <- cutree(hc, k = k)
          wcss[k] <- .dispersion(cl)
          if (k >= k_min) {
            sil <- cluster::silhouette(cl, d)
            sil_scores[k] <- mean(sil[, "sil_width"])
          }
        }
      })
      
      # Elbow: knee point via second derivative of within-cluster dispersion
      wcss_diff2   <- diff(diff(wcss))
      elbow_k      <- which.max(wcss_diff2) + 1
      
      # Silhouette: k with highest average silhouette width
      sil_valid    <- sil_scores[k_min:k_max]
      silhouette_k <- which.max(sil_valid) + k_min - 1
      
      list(
        elbow_k      = elbow_k,
        silhouette_k = silhouette_k,
        df_elbow     = data.frame(k = 1:k_max,     wcss       = wcss),
        df_sil       = data.frame(k = k_min:k_max, silhouette = sil_valid),
        k_max        = k_max,
        k_min        = k_min
      )
    })
    
    # Dynamic UI: suggestion text + two small plots
    output$panel <- renderUI({
      res <- result()
      if (is.null(res)) {
        return(tags$p(
          style = "color:gray; font-size:0.82em; font-style:italic; margin-top:8px;",
          "No k-suggestion available for the current selection (too few genes with ",
          "variation, or no groups assigned yet). Check the session log for details."
        ))
      }
      ns <- session$ns
      tagList(
        tags$div(
          style = "margin-top: 8px; padding: 6px; background: #f8f9fa;
                   border-left: 3px solid #888; border-radius: 3px;",
          tags$small(tags$b("\U0001f50d Suggested k:")),
          tags$p(
            style = "margin: 2px 0; font-size: 0.82em; color: #c0392b;",
            paste0("Elbow method \u2192 k = ", res$elbow_k)
          ),
          tags$p(
            style = "margin: 2px 0; font-size: 0.82em; color: #27ae60;",
            paste0("Silhouette score \u2192 k = ", res$silhouette_k)
          ),
          fluidRow(
            column(6,
                   tags$div(
                     title = "Click to enlarge",
                     style = "cursor: pointer;",
                     shinycssloaders::withSpinner(
                       plotOutput(ns("plot_elbow"), height = "160px",
                                  click = ns("click_elbow")),
                       type = 8, color = "#5B9BD5", size = 0.6
                     )
                   )
            ),
            column(6,
                   tags$div(
                     title = "Click to enlarge",
                     style = "cursor: pointer;",
                     shinycssloaders::withSpinner(
                       plotOutput(ns("plot_sil"), height = "160px",
                                  click = ns("click_sil")),
                       type = 8, color = "#E07B54", size = 0.6
                     )
                   )
            )
          )
        )
      )
    })
    
    # Elbow plot
    output$plot_elbow <- renderPlot({
      res <- result()
      req(res)
      df <- res$df_elbow
      ggplot2::ggplot(df, ggplot2::aes(x = k, y = wcss)) +
        ggplot2::geom_line(color = "steelblue", linewidth = 0.8) +
        ggplot2::geom_point(size = 2, color = "steelblue") +
        ggplot2::geom_vline(xintercept = res$elbow_k,
                            linetype = "dashed", color = "firebrick", linewidth = 0.7) +
        ggplot2::annotate("text", x = res$elbow_k + 0.3, y = max(df$wcss) * 0.95,
                          label = paste("k =", res$elbow_k),
                          color = "firebrick", hjust = 0, size = 2.8) +
        ggplot2::scale_x_continuous(breaks = 1:res$k_max) +
        ggplot2::labs(title = "Elbow", x = "k", y = "WCSS") +
        ggplot2::theme_minimal(base_size = 8)
    }, bg = "transparent")
    
    # Silhouette plot
    output$plot_sil <- renderPlot({
      res <- result()
      req(res)
      df <- res$df_sil
      ggplot2::ggplot(df, ggplot2::aes(x = k, y = silhouette)) +
        ggplot2::geom_line(color = "darkorange", linewidth = 0.8) +
        ggplot2::geom_point(size = 2, color = "darkorange") +
        ggplot2::geom_vline(xintercept = res$silhouette_k,
                            linetype = "dashed", color = "darkgreen", linewidth = 0.7) +
        ggplot2::annotate("text", x = res$silhouette_k + 0.3, y = max(df$silhouette) * 0.99,
                          label = paste("k =", res$silhouette_k),
                          color = "darkgreen", hjust = 0, size = 2.8) +
        ggplot2::scale_x_continuous(breaks = res$k_min:res$k_max) +
        ggplot2::labs(title = "Silhouette", x = "k", y = "Avg width") +
        ggplot2::theme_minimal(base_size = 8)
    }, bg = "transparent")
    
    # Open elbow plot in full-size modal on click
    observeEvent(input$click_elbow, {
      res <- result()
      req(res)
      df <- res$df_elbow
      p  <- ggplot2::ggplot(df, ggplot2::aes(x = k, y = wcss)) +
        ggplot2::geom_line(color = "steelblue", linewidth = 1) +
        ggplot2::geom_point(size = 3, color = "steelblue") +
        ggplot2::geom_vline(xintercept = res$elbow_k,
                            linetype = "dashed", color = "firebrick", linewidth = 0.8) +
        ggplot2::annotate("text", x = res$elbow_k + 0.3, y = max(df$wcss) * 0.95,
                          label = paste("Elbow k =", res$elbow_k),
                          color = "firebrick", hjust = 0, size = 4) +
        ggplot2::scale_x_continuous(breaks = 1:res$k_max) +
        ggplot2::labs(title = "Elbow Method",
                      x = "Number of Clusters (k)", y = "Within-Cluster Sum of Squares") +
        ggplot2::theme_minimal(base_size = 13)
      
      showModal(modalDialog(
        title     = "Elbow Method — full view",
        size      = "l",
        easyClose = TRUE,
        footer    = modalButton("Close"),
        renderPlot(p, height = 400)
      ))
    })
    
    # Open silhouette plot in full-size modal on click
    observeEvent(input$click_sil, {
      res <- result()
      req(res)
      df <- res$df_sil
      p  <- ggplot2::ggplot(df, ggplot2::aes(x = k, y = silhouette)) +
        ggplot2::geom_line(color = "darkorange", linewidth = 1) +
        ggplot2::geom_point(size = 3, color = "darkorange") +
        ggplot2::geom_vline(xintercept = res$silhouette_k,
                            linetype = "dashed", color = "darkgreen", linewidth = 0.8) +
        ggplot2::annotate("text", x = res$silhouette_k + 0.3, y = max(df$silhouette) * 0.99,
                          label = paste("Best k =", res$silhouette_k),
                          color = "darkgreen", hjust = 0, size = 4) +
        ggplot2::scale_x_continuous(breaks = res$k_min:res$k_max) +
        ggplot2::labs(title = "Silhouette Score",
                      x = "Number of Clusters (k)", y = "Average Silhouette Width") +
        ggplot2::theme_minimal(base_size = 13)
      
      showModal(modalDialog(
        title     = "Silhouette Score — full view",
        size      = "l",
        easyClose = TRUE,
        footer    = modalButton("Close"),
        renderPlot(p, height = 400)
      ))
    })
  })
}
    
