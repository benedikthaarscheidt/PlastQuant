# =============================================================================
# snp_genotypes.R — SNP genotypes, causal SNPs and the SNP score of a genotype
# =============================================================================

# Dosages x_is ~ Bin(2, p_s) with allele frequencies p_s ~ U(0.2, 0.8).
# Filled in column blocks to keep memory low for large SNP numbers.
simulate_snp_genotypes <- function(n_genotypes, n_snps, block = 10000) {
  p <- runif(n_snps, 0.2, 0.8)
  X <- matrix(0L, n_genotypes, n_snps,
    dimnames = list(paste0("Genotype_", seq_len(n_genotypes)), paste0("SNP_", seq_len(n_snps)))
  )
  for (b in split(seq_len(n_snps), ceiling(seq_len(n_snps) / block))) {
    X[, b] <- rbinom(n_genotypes * length(b), 2, rep(p[b], each = n_genotypes))
  }
  list(X = X, p = p)
}

# Mc causal SNPs for every (shape, parameter) pair incl. the repeatability trait V,
# drawn without replacement so that the sets are disjoint; effects are +1 or -1.
draw_causal_snps <- function(cfg) {
  pairs <- do.call(rbind, lapply(cfg$shapes, function(s) {
    data.frame(Shape = s, Parameter = c(reaction_norm_parameters(cfg, s), "V"), stringsAsFactors = FALSE)
  }))
  need <- nrow(pairs) * cfg$Mc
  if (need > cfg$n_snps) stop("need ", need, " causal SNPs but only ", cfg$n_snps, " SNPs simulated")
  idx <- sample.int(cfg$n_snps, need)
  data.frame(
    Shape     = rep(pairs$Shape, each = cfg$Mc),
    Parameter = rep(pairs$Parameter, each = cfg$Mc),
    SNP_index = idx,
    SNP       = paste0("SNP_", idx),
    Effect    = sample(c(-1L, 1L), need, replace = TRUE),
    stringsAsFactors = FALSE
  )
}

# SNP score S_i = (1/sqrt(Mc)) * sum_s beta_s (x_is - 2 p_s) / sqrt(2 p_s (1 - p_s))
snp_score <- function(X, p, rows, snp_idx, beta) {
  ps <- p[snp_idx]
  Z <- sweep(X[rows, snp_idx, drop = FALSE], 2, 2 * ps, "-")
  Z <- sweep(Z, 2, sqrt(2 * ps * (1 - ps)), "/")
  as.vector(Z %*% beta) / sqrt(length(snp_idx))
}
