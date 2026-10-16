# =============================================================================
# evaluate_power.R — how many causal SNPs each PPI recovers, per shape
# =============================================================================

# sig_hits: significant hits of ONE shape (Score, SNP); truth_shape: causal SNPs of
# the SAME shape (Parameter, SNP). Hits on SNPs that are not causal for this shape
# count as false positives.
# With LD (window w > 0, snp_pos: position of every SNP, named by SNP) a causal SNP is
# found if a hit lies within w SNPs of it, and the false positives are counted as loci:
# hits further than w from every causal SNP of the shape, where hits at most w apart
# form one locus. w = 0 is the exact-SNP rule above.
# On a real SNP panel (tags: table SNP, Tag of causal_tags.csv, see snp_panel.R) a causal
# SNP is found if a hit is one of its tag SNPs; the false positives are the hits that tag
# no causal SNP of the shape, counted as loci (hits with r^2 >= tag_r2 to each other form
# one locus; needs snp_matrix, otherwise every hit counts).
compute_power <- function(sig_hits, truth_shape, scores, window = 0, snp_pos = NULL,
                          tags = NULL, snp_matrix = NULL, tag_r2 = 0.1) {
  if (window > 0 && is.null(snp_pos)) stop("an LD window needs the SNP positions (snp_pos)")
  tag_sets <- if (!is.null(tags)) split(tags$Tag, tags$SNP)
  rows <- list()
  for (s in scores) {
    hits <- unique(sig_hits$SNP[sig_hits$Score == s])
    for (p in unique(truth_shape$Parameter)) {
      true_snps <- truth_shape$SNP[truth_shape$Parameter == p]
      found <- if (!is.null(tags)) {
        sum(vapply(true_snps, function(cs) any(tag_sets[[cs]] %in% hits), logical(1)))
      } else if (window == 0) {
        length(intersect(hits, true_snps))
      } else {
        sum(vapply(snp_pos[true_snps], function(x) any(abs(snp_pos[hits] - x) <= window), logical(1)))
      }
      rows[[length(rows) + 1]] <- data.frame(Score = s, Parameter = p, True_SNPs = length(true_snps),
                                             Found = found, Power = found / length(true_snps),
                                             stringsAsFactors = FALSE)
    }
    false_hits <- if (!is.null(tags)) {
      count_loci(setdiff(hits, unlist(tag_sets[truth_shape$SNP])), snp_matrix, tag_r2)
    } else if (window == 0) {
      length(setdiff(hits, truth_shape$SNP))
    } else {
      x <- sort(snp_pos[hits])
      causal <- snp_pos[truth_shape$SNP]
      x <- x[vapply(x, function(h) all(abs(causal - h) > window), logical(1))]
      if (length(x) == 0) 0L else 1L + sum(diff(x) > window)   # loci of the remaining hits
    }
    rows[[length(rows) + 1]] <- data.frame(Score = s, Parameter = "False_Positives", True_SNPs = NA_integer_,
                                           Found = false_hits, Power = NA_real_,
                                           stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}

# Number of loci among the SNPs `snps`: greedily, the first remaining SNP and all SNPs with
# r^2 >= r2 to it form one locus. Without snp_matrix every SNP is its own locus.
count_loci <- function(snps, snp_matrix, r2) {
  if (is.null(snp_matrix)) return(length(snps))
  rest <- match(snps, colnames(snp_matrix))
  rest <- sort(rest[!is.na(rest)])
  n <- 0L
  while (length(rest)) {
    rest <- rest[snp_r2(snp_matrix, rest[1], rest) < r2]   # removes rest[1] itself (r^2 = 1)
    n <- n + 1L
  }
  n
}

# Bonferroni threshold alpha / number of SNPs; one block of rows per shape.
# window, snp_pos, tags, snp_matrix, tag_r2: hits in LD with a causal SNP, see compute_power()
evaluate_power <- function(gwas_results, truth, n_snps, alpha = 0.05, window = 0, snp_pos = NULL,
                           tags = NULL, snp_matrix = NULL, tag_r2 = 0.1) {
  sig <- gwas_results[gwas_results$P_Value < alpha / n_snps, ]
  do.call(rbind, lapply(unique(gwas_results$Type), function(shape) {
    # only the PPIs whose GWAS was run for this shape
    scores <- unique(gwas_results$Score[gwas_results$Type == shape])
    res <- compute_power(sig[sig$Type == shape, ], truth[truth$Shape == shape, ], scores, window, snp_pos,
                         tags, snp_matrix, tag_r2)
    res$Shape <- shape
    res
  }))
}
