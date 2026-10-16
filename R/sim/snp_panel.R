# =============================================================================
# snp_panel.R — simulate the reaction norms on the SNPs of a real panel
# =============================================================================
# simulate_population_from_snps(cfg, snps, dir) is simulate_population() with the lines of
# a real SNP panel as genotypes instead of simulated SNPs, so that LD, allele frequencies
# and relatedness are those of a real population. It differs from simulate_norms() only in:
#   - every shape is simulated on all lines of the panel (genotype IDs = line IDs); the
#     shapes stay separate (own reaction norms, PPIs and GWAS)
#   - causal SNPs: allele frequency MAF >= min_maf, and within a shape no two causal SNPs
#     with r^2 >= max_r2. L, a and V have the same causal SNPs (and effects) in every
#     shape, so the shapes are comparable; the shape parameters (c, w, k, psi) get their own
#   - the dosages of a causal SNP are standardised with their mean and SD over the lines
#     (no Hardy-Weinberg assumption: inbred lines carry almost only 0 and 2)
#   - tag SNPs: for every causal SNP the SNPs with r^2 >= max_r2 to it, written to
#     data/causal_tags.csv; a GWAS hit on a tag SNP finds the causal SNP (evaluate_power.R)
# From data/ on, the population folder is the same as for a simulated population, and
# sample_resolution(), compute_population_ppis() and run_population_gwas() run unchanged.
# =============================================================================

# SNP panel from a matrix or a file (.rds with the matrix, or .csv with the column
# Genotype and one column per SNP): lines x SNPs, dosages 0/1/2, without missing values
read_snp_panel <- function(snps) {
  X <- if (is.character(snps)) {
    if (grepl("\\.rds$", snps, ignore.case = TRUE)) {
      readRDS(snps)
    } else {
      tab <- read_genotype_table(snps)
      m <- as.matrix(tab[, -1, drop = FALSE])
      rownames(m) <- tab$Genotype
      m
    }
  } else snps
  if (!is.matrix(X) || is.null(rownames(X)) || is.null(colnames(X))) {
    stop("the SNP panel must be a matrix with line IDs as rownames and SNP IDs as colnames")
  }
  if (!is_dosage_matrix(X)) stop("the SNP panel must contain dosages 0, 1, 2 without missing values")
  X
}

# r^2 between the dosages of SNP j and the SNPs `cols`, in column blocks
snp_r2 <- function(X, j, cols = seq_len(ncol(X)), block = 50000) {
  x <- X[, j]
  unlist(lapply(split(cols, ceiling(seq_along(cols) / block)), function(b) {
    r <- suppressWarnings(cor(x, X[, b, drop = FALSE]))
    r[is.na(r)] <- 0   # monomorphic SNPs
    as.vector(r)^2
  }))
}

# n SNPs out of `candidates`, in random order, none with r^2 >= max_r2 to another one or
# to the SNPs already `taken`
draw_unlinked_snps <- function(X, n, candidates, taken, max_r2) {
  chosen <- integer(0)
  for (s in candidates[sample.int(length(candidates))]) {
    others <- c(taken, chosen)
    if (length(others) == 0 || max(snp_r2(X, s, others)) < max_r2) chosen <- c(chosen, s)
    if (length(chosen) == n) return(chosen)
  }
  stop("cannot draw ", n, " causal SNPs with r^2 < ", max_r2, " to each other; use fewer causal SNPs (Mc)")
}

# Causal SNPs of every (shape, parameter) pair: L, a and V shared by all shapes, the shape
# parameters drawn per shape; effects +1 or -1. Columns as draw_causal_snps().
draw_panel_causal_snps <- function(X, cfg, min_maf, max_r2) {
  f <- colMeans(X) / 2
  eligible <- which(pmin(f, 1 - f) >= min_maf)
  shared <- c("L", "a", "V")
  shared_idx <- draw_unlinked_snps(X, length(shared) * cfg$Mc, eligible, integer(0), max_r2)
  shared_effect <- sample(c(-1L, 1L), length(shared_idx), replace = TRUE)
  used <- shared_idx
  do.call(rbind, lapply(cfg$shapes, function(s) {
    own <- setdiff(reaction_norm_parameters(cfg, s), shared)
    own_idx <- if (length(own)) draw_unlinked_snps(X, length(own) * cfg$Mc, setdiff(eligible, used), shared_idx, max_r2) else integer(0)
    used <<- c(used, own_idx)
    idx <- c(shared_idx, own_idx)
    data.frame(
      Shape     = s,
      Parameter = rep(c(shared, own), each = cfg$Mc),
      SNP_index = idx,
      SNP       = colnames(X)[idx],
      Effect    = c(shared_effect, sample(c(-1L, 1L), length(own_idx), replace = TRUE)),
      stringsAsFactors = FALSE
    )
  }))
}

# Tag SNPs: every SNP with r^2 >= max_r2 to a causal SNP (the causal SNP itself included)
causal_tag_snps <- function(X, causal_idx, max_r2) {
  do.call(rbind, lapply(unique(causal_idx), function(j) {
    r2 <- snp_r2(X, j)
    tags <- which(r2 >= max_r2)
    data.frame(SNP = colnames(X)[j], Tag = colnames(X)[tags], R2 = r2[tags], stringsAsFactors = FALSE)
  }))
}

# The norms of a population on a SNP panel; same result as simulate_norms()
simulate_norms_from_snps <- function(cfg, X, min_maf = 0.05, max_r2 = 0.1) {
  cfg$n_per_shape <- nrow(X)
  cfg$n_snps <- ncol(X)
  cfg$ld_strength <- 0
  cfg$ld_window <- 0L
  cfg$snp_panel <- list(min_maf = min_maf, max_r2 = max_r2)

  # Genetics from their own seed, as in simulate_norms()
  set.seed(cfg$seed)
  causal <- draw_panel_causal_snps(X, cfg, min_maf, max_r2)
  # p = NULL: snp_score() standardises with the mean and SD of each SNP over the lines
  geno <- list(X = X, p = NULL)

  set.seed(cfg$seed + 1)
  rows <- seq_len(nrow(X))
  per_shape <- lapply(cfg$shapes, function(shape) simulate_shape_norms(cfg, shape, geno, causal, rows))
  names(per_shape) <- cfg$shapes

  list(
    snp_matrix        = X,
    truth             = causal[, c("Shape", "Parameter", "SNP", "Effect")],
    causal_tags       = causal_tag_snps(X, causal$SNP_index, max_r2),
    mother_parameters = do.call(rbind, lapply(per_shape, `[[`, "mother_parameters")),
    plant_parameters  = lapply(per_shape, `[[`, "plant_parameters"),
    latent            = lapply(per_shape, `[[`, "latent"),
    config            = cfg
  )
}

# Simulate the true reaction norms on the lines of a SNP panel (matrix or file, see
# read_snp_panel()) and save them as a population folder, as simulate_population() does,
# plus data/causal_tags.csv. cfg: sim_config(); its n_per_shape, n_snps and ld_strength
# are replaced by the panel. Returns dir.
simulate_population_from_snps <- function(cfg, snps, dir, min_maf = 0.05, max_r2 = 0.1) {
  X <- read_snp_panel(snps)
  norms <- simulate_norms_from_snps(cfg, X, min_maf, max_r2)
  rm(X)
  write_simulated_population(norms, dir)
}
