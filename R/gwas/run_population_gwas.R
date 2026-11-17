# =============================================================================
# run_population_gwas.R — Step 4 (optional): the GWAS of the PPIs of a resolution
# =============================================================================
# Steps on a population folder (see population_folder.R and the README):
#   simulate_population() -> sample_resolution() -> compute_population_ppis() -> run_population_gwas()

# GWAS of every PPI of a resolution and replicate handling (dir/<replicate_handling>/ppis)
# against the population's SNPs -> dir/<replicate_handling>/gwas/gwas_results.csv
# (kinship from the population's data/kinship.csv if present, otherwise from the
# population's SNPs); with a truth table also the power per shape (Bonferroni alpha /
# number of tested SNPs; with LD per LD window, see evaluate_power.R) ->
# gwas_power.csv and its heatmaps -> dir/<replicate_handling>/plots/gwas_power_heatmap.pdf.
# Always: the minor allele frequency of the significant hits (Bonferroni) and, with a
# truth table, of the causal SNPs (hit_maf.R) -> gwas/gwas_hit_maf.csv,
# gwas/gwas_causal_maf.csv and plots/gwas_hit_maf.pdf.
#   gwas_snps        NULL: test the population's SNPs (data/snp_matrix.rds). Otherwise
#                    a larger SNP set of the same lines (matrix or file, see
#                    read_snp_panel(): lines x SNPs, dosages 0/1/2, line IDs as rownames),
#                    e.g. the whole panel before pruning; it must hold every line of the
#                    population. On a SNP panel the tag SNPs of the causal SNPs are then
#                    recomputed on this set (r^2 >= max_r2 of the panel). The kinship
#                    stays that of the population's SNPs.
#   min_maf          SNPs with a minor allele frequency below it (over the lines of the
#                    population) are not tested, before any draw (n_gwas_snps); a causal
#                    SNP below it is an error, as every causal SNP must be tested
#   true_parameters  TRUE: the parameters of the mother reaction norms are GWAS traits
#                    next to the PPIs, named true_<parameter> (true_parameter_traits())
#   n_gwas_snps      number of SNPs tested: all causal SNPs plus a random draw of the
#                    other SNPs (MAF >= min_maf) up to this number, in genomic order; the
#                    draw keeps the LD structure of the SNP set at that density and comes
#                    from its own stream (seed + 3 of the simulation). Inf, or a number
#                    above the SNPs available: all of them.
# n_cores > 1 runs the SNP tests on that many cores (forked processes: from the terminal
# with Rscript, not in RStudio); the results are the same.
run_population_gwas <- function(dir, replicate_handling = c("mean_before", "replicates"), alpha = 0.05, n_cores = 1,
                                gwas_snps = NULL, min_maf = 0.05, true_parameters = FALSE, n_gwas_snps = Inf) {
  replicate_handling <- match.arg(replicate_handling)
  set_dir <- ppi_set_path(dir, replicate_handling)
  ppi_file <- file.path(set_dir, "ppis", "ppi_scores.rds")
  if (!file.exists(ppi_file)) {
    stop("no PPIs in ", set_dir, " (run compute_population_ppis(dir, replicate_handling = \"", replicate_handling, "\") first)")
  }
  pop <- read_population(dir)
  if (is.null(pop$snp_matrix) && is.null(gwas_snps)) {
    stop("the GWAS needs a SNP matrix in ", file.path(pop$population_dir, "data"), " (snp_matrix.csv) or gwas_snps")
  }
  traits <- readRDS(ppi_file)
  if (true_parameters) {
    if (is.null(pop$mother_parameters)) stop("true_parameters = TRUE needs a simulated population (data/mother_parameters.csv)")
    traits <- c(traits, true_parameter_traits(pop$mother_parameters))
  }

  # the SNPs tested: the population's own or a larger set of the same lines
  if (is.null(gwas_snps)) {
    X <- pop$snp_matrix
  } else {
    if (!is.null(pop$truth) && is.null(pop$causal_tags)) {
      stop("gwas_snps is for populations on a SNP panel (simulate_population_from_snps()) or own data, ",
           "not for populations with simulated SNPs")
    }
    X <- read_snp_panel(gwas_snps)
    lines <- if (is.null(pop$snp_matrix)) unique(unlist(lapply(traits, `[[`, "Genotype_ID"))) else rownames(pop$snp_matrix)
    missing <- setdiff(lines, rownames(X))
    if (length(missing)) stop("gwas_snps lacks ", length(missing), " lines of the population, e.g. ", paste(head(missing, 5), collapse = ", "))
    X <- X[lines, , drop = FALSE]
  }
  maf <- snp_maf(X)
  X <- X[, maf >= min_maf, drop = FALSE]
  maf <- maf[maf >= min_maf]
  # the causal SNPs are always tested: none of them may have MAF < min_maf or be missing
  causal <- if (is.null(pop$truth)) character(0) else unique(pop$truth$SNP)
  lost <- setdiff(causal, colnames(X))
  if (length(lost)) {
    stop(length(lost), " causal SNPs have MAF < ", min_maf, " over the lines of the population or are not in ",
         "the SNP set, e.g. ", paste(head(lost, 5), collapse = ", "))
  }
  # a given number of SNPs: the causal SNPs plus a random draw of the others
  if (n_gwas_snps < ncol(X)) {
    if (n_gwas_snps < length(causal)) {
      warning("n_gwas_snps (", n_gwas_snps, ") is below the number of causal SNPs (", length(causal), "); only these are tested")
    }
    set.seed(if (is.null(pop$config$seed)) 3 else pop$config$seed + 3)
    others <- setdiff(colnames(X), causal)
    keep <- colnames(X) %in% c(causal, others[sample.int(length(others), max(0, n_gwas_snps - length(causal)))])
    X <- X[, keep, drop = FALSE]
    maf <- maf[keep]
  }
  n_snps <- ncol(X)
  message("GWAS of ", length(traits), " traits on ", n_snps, " SNPs (MAF >= ", min_maf, ")")

  gwas_dir <- file.path(set_dir, "gwas")
  unlink(gwas_dir, recursive = TRUE)
  dir.create(gwas_dir, recursive = TRUE, showWarnings = FALSE)
  kin <- if (!is.null(pop$kinship)) pop$kinship else kinship_matrix(if (is.null(pop$snp_matrix)) X else pop$snp_matrix)
  res <- run_gwas(traits, X, kin, n_cores, min_maf)
  write.csv(res, file.path(gwas_dir, "gwas_results.csv"), row.names = FALSE)

  cfg <- pop$config
  threshold <- alpha / n_snps
  subtitle <- paste0(
    if (!is.null(cfg)) paste0("h2 = ", cfg$h2, if (!is.null(cfg$q)) paste0(", q = ", cfg$q), ", ",
                              cfg$Mc, " causal SNPs per parameter, ",
                              cfg$n_per_shape, " genotypes per shape, ", if (!is.null(cfg$seed)) paste0("seed ", cfg$seed, ", ")),
    n_snps, " SNPs (MAF >= ", min_maf, if (!is.null(gwas_snps)) ", external SNP set", "), ",
    replicate_handling, "  (Bonferroni P < ", format(threshold, digits = 2, scientific = TRUE), ")"
  )
  window <- 0
  snp_pos <- NULL
  tags <- NULL
  if (!is.null(pop$truth)) {
    # LD window of a simulation with LD (0 otherwise): hits within it of a causal SNP find it;
    # positions in the population's SNP matrix, as simulated
    window <- if (is.null(cfg$ld_window)) 0 else cfg$ld_window
    if (!is.null(pop$snp_matrix)) snp_pos <- setNames(seq_len(ncol(pop$snp_matrix)), colnames(pop$snp_matrix))
    # on a real SNP panel: hits on the tag SNPs of a causal SNP find it (data/causal_tags.csv,
    # or recomputed on the tested SNPs when they are not the population's)
    tag_r2 <- if (is.null(cfg$snp_panel)) 0.1 else cfg$snp_panel$max_r2
    tags <- pop$causal_tags
    if (!is.null(tags) && !is.null(gwas_snps)) {
      tags <- causal_tag_snps(X, match(causal, colnames(X)), tag_r2)
      write.csv(tags, file.path(gwas_dir, "causal_tags_tested_snps.csv"), row.names = FALSE)
    }
    power <- evaluate_power(res, pop$truth, n_snps, alpha = alpha, window = window, snp_pos = snp_pos,
                            tags = tags, snp_matrix = X, tag_r2 = tag_r2)
    write.csv(power, file.path(gwas_dir, "gwas_power.csv"), row.names = FALSE)
    power_subtitle <- paste0(subtitle, if (window > 0) paste0(", LD window ", window, " SNPs"),
                             if (!is.null(tags)) ", SNP panel (hits on tag SNPs)")
    unlink(file.path(set_dir, "plots", "gwas_power_heatmap.pdf"))
    plot_gwas_power(power, file.path(set_dir, "plots", "gwas_power_heatmap.pdf"), power_subtitle)
  }

  # minor allele frequency of the hits and of the causal SNPs
  maf_tabs <- hit_maf_tables(res, maf, threshold, pop$truth, tags, window, snp_pos)
  write.csv(maf_tabs$hits, file.path(gwas_dir, "gwas_hit_maf.csv"), row.names = FALSE)
  if (!is.null(maf_tabs$causal)) write.csv(maf_tabs$causal, file.path(gwas_dir, "gwas_causal_maf.csv"), row.names = FALSE)
  unlink(file.path(set_dir, "plots", "gwas_hit_maf.pdf"))
  plot_hit_maf(maf_tabs, maf, file.path(set_dir, "plots", "gwas_hit_maf.pdf"), subtitle)

  message("GWAS written to ", gwas_dir)
  dir
}
