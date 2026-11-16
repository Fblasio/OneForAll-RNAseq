validate_rnaseq_sample_names <- function(sample_names) {
  sample_names <- as.character(sample_names)
  
  missing_separator <- sample_names[!grepl("_", sample_names, fixed = TRUE)]
  
  token_lengths <- vapply(
    strsplit(sample_names, "_", fixed = TRUE),
    length,
    integer(1)
  )
  
  too_short <- sample_names[token_lengths < 2]
  
  replicate_tokens <- vapply(
    strsplit(sample_names, "_", fixed = TRUE),
    function(x) {
      x <- trimws(x)
      if (length(x) == 0) return(NA_character_)
      x[length(x)]
    },
    character(1)
  )
  
  empty_replicate <- sample_names[is.na(replicate_tokens) | !nzchar(replicate_tokens)]
  
  list(
    has_warning = length(missing_separator) > 0 ||
      length(too_short) > 0 ||
      length(empty_replicate) > 0,
    missing_separator = unique(missing_separator),
    too_short = unique(too_short),
    empty_replicate = unique(empty_replicate)
  )
}

build_rnaseq_metadata_table <- function(sample_names, mode = "separate") {
  sample_names <- as.character(sample_names)
  
  parsed_list <- lapply(sample_names, function(x) {
    tokens <- strsplit(x, "_", fixed = TRUE)[[1]]
    tokens <- trimws(tokens)
    
    if (length(tokens) < 2) {
      return(data.frame(
        sample = x,
        factor1 = tokens[1],
        factor2 = NA_character_,
        combined_factor = tokens[1],
        replicate = NA_character_,
        observations = NA_character_,
        stringsAsFactors = FALSE
      ))
    }
    
    replicate <- tokens[length(tokens)]
    prefix_tokens <- tokens[-length(tokens)]
    
    factor1 <- if (length(prefix_tokens) >= 1) prefix_tokens[1] else NA_character_
    factor2 <- if (length(prefix_tokens) >= 2) prefix_tokens[2] else NA_character_
    
    observations <- if (length(prefix_tokens) > 2) {
      paste(prefix_tokens[3:length(prefix_tokens)], collapse = "_")
    } else {
      NA_character_
    }
    
    combined_factor <- if (!is.na(factor2) && nzchar(factor2)) {
      paste(factor1, factor2, sep = "_")
    } else {
      factor1
    }
    
    data.frame(
      sample = x,
      factor1 = factor1,
      factor2 = factor2,
      combined_factor = combined_factor,
      replicate = replicate,
      observations = observations,
      stringsAsFactors = FALSE
    )
  })
  
  meta <- do.call(rbind, parsed_list)
  
  if (all(is.na(meta$factor2) | meta$factor2 == "")) {
    meta$factor2 <- NULL
  }
  
  if (all(is.na(meta$observations) | meta$observations == "")) {
    meta$observations <- NULL
  }
  
  if (identical(mode, "combined")) {
    meta$factor1 <- NULL
    if ("factor2" %in% colnames(meta)) meta$factor2 <- NULL
    meta$combined_factor <- as.character(meta$combined_factor)
  } else {
    meta$combined_factor <- NULL
    meta$factor1 <- as.character(meta$factor1)
    if ("factor2" %in% colnames(meta)) {
      meta$factor2 <- as.character(meta$factor2)
    }
  }
  
  meta$replicate <- as.character(meta$replicate)
  
  meta
}

build_empty_rnaseq_metadata_table <- function(sample_names, mode = "separate") {
  sample_names <- as.character(sample_names)
  
  if (identical(mode, "combined")) {
    data.frame(
      sample = sample_names,
      combined_factor = NA_character_,
      replicate = NA_character_,
      observations = NA_character_,
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      sample = sample_names,
      factor1 = NA_character_,
      factor2 = NA_character_,
      replicate = NA_character_,
      observations = NA_character_,
      stringsAsFactors = FALSE
    )
  }
}