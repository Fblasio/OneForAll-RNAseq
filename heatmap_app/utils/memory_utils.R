check_memory <- function(n_genes, n_samples) {
  msgs <- c()
  if (n_genes > MAX_GENES_MEMORY)
    msgs <- c(msgs, paste0("⚠️ ", n_genes, " genes — consider filtering to < ",
                           MAX_GENES_MEMORY, " for better performance."))
  est_mb <- round((n_genes * n_samples * 8) / 1024^2)
  if (est_mb > 200)
    msgs <- c(msgs, paste0("⚠️ Estimated memory: ~", est_mb,
                           " MB. Clustering may be slow."))
  msgs
}