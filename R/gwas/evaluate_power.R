# =============================================================================
# evaluate_power.R — how many causal SNPs each PPI recovers, per shape
# =============================================================================

# sig_hits: significant hits of ONE shape (Score, SNP); truth_shape: causal SNPs of
# the SAME shape (Parameter, SNP). Hits on SNPs that are not causal for this shape
# count as false positives.
compute_power <- function(sig_hits, truth_shape, scores) {
  rows <- list()
  for (s in scores) {
    hits <- unique(sig_hits$SNP[sig_hits$Score == s])
    for (p in unique(truth_shape$Parameter)) {
      true_snps <- truth_shape$SNP[truth_shape$Parameter == p]
      found <- length(intersect(hits, true_snps))
      rows[[length(rows) + 1]] <- data.frame(Score = s, Parameter = p, True_SNPs = length(true_snps),
                                             Found = found, Power = found / length(true_snps),
                                             stringsAsFactors = FALSE)
    }
    rows[[length(rows) + 1]] <- data.frame(Score = s, Parameter = "False_Positives", True_SNPs = NA_integer_,
                                           Found = length(setdiff(hits, truth_shape$SNP)), Power = NA_real_,
                                           stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}

# Bonferroni threshold alpha / number of SNPs; one block of rows per shape
evaluate_power <- function(gwas_results, truth, n_snps, alpha = 0.05) {
  sig <- gwas_results[gwas_results$P_Value < alpha / n_snps, ]
  do.call(rbind, lapply(unique(gwas_results$Type), function(shape) {
    # only the PPIs whose GWAS was run for this shape
    scores <- unique(gwas_results$Score[gwas_results$Type == shape])
    res <- compute_power(sig[sig$Type == shape, ], truth[truth$Shape == shape, ], scores)
    res$Shape <- shape
    res
  }))
}
