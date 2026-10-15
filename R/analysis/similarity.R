# =============================================================================
# similarity.R — similarity of two vectors and its permutation test
# =============================================================================
# Used by the analyses that relate plasticity scores to summary statistics.

# Original similarity function using various methods.
calc_similarity <- function(x, y, method = "pearson") {
  if (method == "pearson") {
    return(cor(x, y, method = "pearson", use = "complete.obs"))
  } else if (method == "spearman") {
    return(cor(x, y, method = "spearman", use = "complete.obs"))
  } else if (method == "kendall") {
    return(cor(x, y, method = "kendall", use = "complete.obs"))
  } else if (method == "dcor") {
    return(energy::dcor(x, y))
  } else if (method == "cosine") {
    return(sum(x * y) / (sqrt(sum(x^2)) * sqrt(sum(y^2))))
  } else {
    stop("Unknown method: ", method)
  }
}

calc_similarity_perm <- function(x, y, method = "pearson", n_perm = 999) {
  # Convert inputs to numeric (if not already)
  x <- as.numeric(x)
  y <- as.numeric(y)
  
  # Determine indices with complete cases
  complete_idx <- complete.cases(x, y)
  if (sum(complete_idx) < 2) {
    return(list(estimate = NA, p.value = NA))
  }
  
  # Use only complete pairs
  x_complete <- x[complete_idx]
  y_complete <- y[complete_idx]
  
  # Compute the observed similarity measure
  obs <- calc_similarity(x_complete, y_complete, method)
  
  # Perform permutation test by shuffling y_complete
  perm_stats <- replicate(n_perm, {
    permuted_y <- sample(y_complete, length(y_complete), replace = FALSE)
    calc_similarity(x_complete, permuted_y, method)
  })
  
  # Two-sided test: count how many permuted similarities are as extreme as observed.
  p_val <- (sum(abs(perm_stats) >= abs(obs)) + 1) / (n_perm + 1)
  
  return(list(estimate = obs, p.value = p_val))
}
