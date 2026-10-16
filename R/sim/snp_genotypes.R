# =============================================================================
# snp_genotypes.R — SNP genotypes, causal SNPs and the SNP score of a genotype
# =============================================================================

# Dosages of every genotype at every SNP, allele frequencies p_s ~ U(0.2, 0.8).
# ld_strength = 0: independent SNPs, x_is ~ Bin(2, p_s), filled in column blocks to keep
#   memory low for large SNP numbers.
# ld_strength = q > 0: linkage disequilibrium along the SNPs. Every genotype has two
#   haplotypes (allele 0 or 1 at every SNP). Along each haplotype, SNP s copies the allele
#   of SNP s - 1 with probability q, otherwise its allele is drawn anew (1 with
#   probability p_s); the dosage is the sum of the two haplotypes. LD between SNPs d apart
#   decays roughly like r^2 = q^(2d). Copying shifts the allele frequency to
#   f_s = q f_(s-1) + (1 - q) p_s, which is returned as p (the dosages are Bin(2, f_s)).
simulate_snp_genotypes <- function(n_genotypes, n_snps, ld_strength = 0, block = 10000) {
  p <- runif(n_snps, 0.2, 0.8)
  X <- matrix(0L, n_genotypes, n_snps,
    dimnames = list(paste0("Genotype_", seq_len(n_genotypes)), paste0("SNP_", seq_len(n_snps)))
  )
  if (ld_strength == 0) {
    for (b in split(seq_len(n_snps), ceiling(seq_len(n_snps) / block))) {
      X[, b] <- rbinom(n_genotypes * length(b), 2, rep(p[b], each = n_genotypes))
    }
    return(list(X = X, p = p))
  }
  n_hap <- 2 * n_genotypes
  first <- seq_len(n_genotypes)            # haplotype 1 of every genotype
  second <- n_genotypes + first            # haplotype 2
  f <- numeric(n_snps)
  hap <- as.integer(runif(n_hap) < p[1])   # alleles of all haplotypes at the current SNP
  f[1] <- p[1]
  X[, 1] <- hap[first] + hap[second]
  for (s in seq_len(n_snps)[-1]) {
    redraw <- which(runif(n_hap) >= ld_strength)
    hap[redraw] <- as.integer(runif(length(redraw)) < p[s])
    f[s] <- ld_strength * f[s - 1] + (1 - ld_strength) * p[s]
    X[, s] <- hap[first] + hap[second]
  }
  list(X = X, p = f)
}

# Mean r^2 between the dosages of SNPs `distances` apart, over n_pairs SNP pairs spread
# evenly along the SNPs (no random numbers, so the simulation streams stay unchanged)
ld_decay <- function(X, distances, n_pairs = 2000) {
  data.frame(Distance = distances, R2 = vapply(distances, function(d) {
    j <- unique(round(seq(1, ncol(X) - d, length.out = min(n_pairs, ncol(X) - d))))
    A <- scale(X[, j, drop = FALSE]); B <- scale(X[, j + d, drop = FALSE])
    mean((colSums(A * B) / (nrow(X) - 1))^2, na.rm = TRUE)
  }, numeric(1)))
}

# LD window: the largest distance (in SNPs) at which the mean r^2 is still >= r2, i.e.
# the first distance with mean r^2 < r2, minus 1; 0 for independent SNPs
ld_window <- function(X, r2 = 0.1, max_distance = 1000) {
  for (d in seq_len(min(max_distance, ncol(X) - 2))) {
    if (ld_decay(X, d)$R2 < r2) return(d - 1L)
  }
  stop("the mean r^2 stays >= ", r2, " up to ", d, " SNPs apart; use weaker LD (ld_strength)")
}

# Mc causal SNPs for every (shape, parameter) pair incl. the repeatability trait V,
# drawn without replacement so that the sets are disjoint; effects are +1 or -1.
# spacing > 0 (with LD): the causal SNPs of one shape are at least `spacing` SNPs apart,
# so they are not in LD with each other and their LD windows do not overlap.
draw_causal_snps <- function(cfg, spacing = 0) {
  pairs <- do.call(rbind, lapply(cfg$shapes, function(s) {
    data.frame(Shape = s, Parameter = c(reaction_norm_parameters(cfg, s), "V"), stringsAsFactors = FALSE)
  }))
  need <- nrow(pairs) * cfg$Mc
  if (need > cfg$n_snps) stop("need ", need, " causal SNPs but only ", cfg$n_snps, " SNPs simulated")
  if (spacing == 0) {
    idx <- sample.int(cfg$n_snps, need)
  } else {
    idx <- integer(0)
    for (s in cfg$shapes) {   # shapes in the order of `pairs`
      n_shape <- sum(pairs$Shape == s) * cfg$Mc
      idx <- c(idx, draw_spaced_snps(n_shape, setdiff(seq_len(cfg$n_snps), idx), spacing))
    }
  }
  data.frame(
    Shape     = rep(pairs$Shape, each = cfg$Mc),
    Parameter = rep(pairs$Parameter, each = cfg$Mc),
    SNP_index = idx,
    SNP       = paste0("SNP_", idx),
    Effect    = sample(c(-1L, 1L), need, replace = TRUE),
    stringsAsFactors = FALSE
  )
}

# n SNPs out of `available`, in random order, each at least `spacing` SNPs from the others
draw_spaced_snps <- function(n, available, spacing) {
  taken <- integer(0)
  for (s in available[sample.int(length(available))]) {
    if (all(abs(s - taken) >= spacing)) taken <- c(taken, s)
    if (length(taken) == n) return(taken)
  }
  stop("cannot place ", n, " causal SNPs at least ", spacing, " SNPs apart; ",
       "use more SNPs, fewer causal SNPs or weaker LD (ld_strength)")
}

# SNP score S_i = (1/sqrt(Mc)) * sum_s beta_s (x_is - 2 p_s) / sqrt(2 p_s (1 - p_s));
# p = NULL (a real SNP panel, see snp_panel.R): the dosages are standardised with the mean
# and SD of each SNP over the genotypes in `rows` instead
snp_score <- function(X, p, rows, snp_idx, beta) {
  if (is.null(p)) {
    Z <- scale(X[rows, snp_idx, drop = FALSE])
  } else {
    ps <- p[snp_idx]
    Z <- sweep(X[rows, snp_idx, drop = FALSE], 2, 2 * ps, "-")
    Z <- sweep(Z, 2, sqrt(2 * ps * (1 - ps)), "/")
  }
  as.vector(Z %*% beta) / sqrt(length(snp_idx))
}
