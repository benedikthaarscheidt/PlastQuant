# =============================================================================
# parameter_ppi_correlation.R — true reaction norm parameters against the PPIs
# =============================================================================
# WHAT IT DOES: For every shape separately, correlates every true reaction norm parameter
#   with every PPI across the reaction norms of one resolution:
#     REPLICATE_HANDLING "replicates"   one point per replicate plant: the plant's own
#                                       parameters (data/plant_parameters.csv) against the
#                                       PPI of that plant's reaction norm
#     REPLICATE_HANDLING "mean_before"  one point per genotype: the mother parameters
#                                       (data/mother_parameters.csv) against the PPI of the
#                                       replicate-averaged reaction norm
#   A population-level PPI (FW) has one value per genotype also with "replicates"; it is
#   correlated with the mother parameters (Level = "genotype").
#   Columns: r (Pearson), rho_s (Spearman), r_within (Pearson after centring parameter and
#   PPI per genotype: the replicate-to-replicate part only; "replicates" only), n.
# REQUIRES:     the PPIs of the resolution and replicate handling (compute_population_ppis())
#               of a simulated population (it needs the true parameters).
# PRODUCES:     RESOLUTION_DIR/<REPLICATE_HANDLING>/analyses/parameter_ppi_correlation.csv
#               RESOLUTION_DIR/<REPLICATE_HANDLING>/plots/parameter_ppi_correlation.pdf
# HOW TO RUN:   setwd("~/PlastQuant")
#               RESOLUTION_DIR <- here::here("output", "populations", "Eve", "n_6")
#               REPLICATE_HANDLING <- "replicates"
#               source(here::here("R", "analysis", "parameter_ppi_correlation.R"))
# -----------------------------------------------------------------------------
# PARAMETERS — set these before sourcing (guarded defaults)
#   RESOLUTION_DIR      path, REQUIRED                             resolution folder n_<k>  [COMMON]
#   REPLICATE_HANDLING  "replicates" (default) or "mean_before"   which PPIs of the resolution
# =============================================================================
source(here::here("R", "pipeline", "load_pipeline.R"))

if (!exists("RESOLUTION_DIR")) stop("set RESOLUTION_DIR (a resolution folder n_<k>) first")
if (!exists("REPLICATE_HANDLING")) stop("set REPLICATE_HANDLING (a replicate handeling) first")
PPI_SET_DIR <- file.path(RESOLUTION_DIR, REPLICATE_HANDLING)
population_dir <- dirname(normalizePath(RESOLUTION_DIR))

ppi_file <- file.path(PPI_SET_DIR, "ppis", "ppi_scores.rds")
if (!file.exists(ppi_file)) {
  stop("no PPIs in ", PPI_SET_DIR, " (run compute_population_ppis(RESOLUTION_DIR, \"", REPLICATE_HANDLING, "\") first)")
}
traits <- readRDS(ppi_file)
plants <- read.csv(file.path(population_dir, "data", "plant_parameters.csv"), stringsAsFactors = FALSE)
mothers <- read.csv(file.path(population_dir, "data", "mother_parameters.csv"), stringsAsFactors = FALSE)
parameter_names <- setdiff(names(plants), c("Shape", "Genotype", "Replicate"))

# Correlations of one PPI with every parameter of one shape; d has the columns Score,
# Genotype and the parameters
correlate <- function(d, shape, ppi, level) {
  do.call(rbind, lapply(parameter_names, function(p) {
    ok <- is.finite(d$Score) & is.finite(d[[p]])
    x <- d[[p]][ok]; y <- d$Score[ok]; g <- d$Genotype[ok]
    if (sum(ok) < 3 || sd(x) == 0 || sd(y) == 0) return(NULL)   # parameter not used by the shape or constant
    r_within <- NA_real_
    if (level == "replicate") {
      xc <- x - ave(x, g); yc <- y - ave(y, g)
      if (sd(xc) > 0 && sd(yc) > 0) r_within <- cor(xc, yc)
    }
    data.frame(Shape = shape, PPI = ppi, Parameter = p, Level = level, n = sum(ok),
               r = cor(x, y), rho_s = cor(x, y, method = "spearman"), r_within = r_within,
               stringsAsFactors = FALSE)
  }))
}

res <- do.call(rbind, lapply(names(traits), function(ppi) {
  tab <- traits[[ppi]]
  
  do.call(rbind, lapply(unique(tab$Type), function(shape) {
    t <- tab[tab$Type == shape, ]
    
    if (all(is.na(t$Replicate))) {   # one value per genotype: mother parameters
      d <- merge(data.frame(Genotype = t$Genotype_ID, Score = t$Score, stringsAsFactors = FALSE),
                 mothers[mothers$Shape == shape, ], by = "Genotype")
      correlate(d, shape, ppi, "genotype")
    } else {                         # one value per replicate plant: its own parameters
      d <- merge(data.frame(Genotype = t$Genotype_ID, Replicate = t$Replicate, Score = t$Score, stringsAsFactors = FALSE),
                 plants[plants$Shape == shape, ], by = c("Genotype", "Replicate"))
      correlate(d, shape, ppi, "replicate")
    }
  }))
}))
rownames(res) <- NULL

out_dir <- population_subdir(PPI_SET_DIR, "analyses")
write.csv(res, file.path(out_dir, "parameter_ppi_correlation.csv"), row.names = FALSE)

# Heatmaps: PPI x parameter per shape, one page for r (and one for r_within with replicates)
plots_dir <- population_subdir(PPI_SET_DIR, "plots")
heat <- function(col, title) {
  d <- res[is.finite(res[[col]]), ]
  d$value <- d[[col]]
  ggplot2::ggplot(d, ggplot2::aes(x = Parameter, y = PPI, fill = value)) +
    ggplot2::geom_tile(colour = "white") +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", value)), size = 2.3) +
    ggplot2::scale_fill_gradient2(low = "#4575B4", mid = "white", high = "#D73027", midpoint = 0, limits = c(-1, 1)) +
    ggplot2::facet_grid(. ~ Shape, scales = "free_x", space = "free_x") +
    ggplot2::labs(title = title, subtitle = paste0(basename(population_dir), ", ", basename(RESOLUTION_DIR), ", ", REPLICATE_HANDLING),
                  x = "True parameter", y = "PPI", fill = col) +
    ggplot2::theme_minimal(base_size = 9)
}
setup <- paste(resolution_setup(RESOLUTION_DIR, REPLICATE_HANDLING),
               if (REPLICATE_HANDLING == "replicates")
                 "One point per replicate plant: its own true parameters against the PPI of its reaction norm (population-level PPIs such as FW: one point per genotype, against the genotype's parameters)."
               else "One point per genotype: its true (mother) parameters against the PPI of its replicate-averaged reaction norm.",
               "Per shape; empty cells: parameter not part of the shape.")
param_reading <- "Shows which true parameter each PPI actually measures: a high correlation with a (relative amplitude) means the PPI tracks the magnitude of plasticity, with L the level, with the shape parameters (c, w, k, psi) the pattern of the response. A PPI that correlates with no parameter cannot find its SNPs in a GWAS."
grDevices::pdf(file.path(plots_dir, "parameter_ppi_correlation.pdf"), width = 14, height = 10)
print(with_caption(heat("r", "Pearson correlation of the true parameters with the PPIs"),
  paste(setup, "Pearson r between each true parameter and each PPI."), param_reading, width = 180))
print(with_caption(heat("rho_s", "Spearman correlation of the true parameters with the PPIs"),
  paste(setup, "Spearman rank correlation between each true parameter and each PPI."),
  paste(param_reading, "Differences from Pearson r point to monotone but non-linear relations."), width = 180))
if (any(is.finite(res$r_within))) print(with_caption(heat("r_within", "Within-genotype correlation (replicate deviations only)"),
  paste(setup, "Pearson r after centring parameter and PPI within each genotype: only the replicate-to-replicate deviations."),
  "Shows whether a PPI also follows the non-genetic differences between replicates of the same genotype; together with r it separates genetic from within-genotype signal.", width = 180))
invisible(grDevices::dev.off())

# --- Groups of PPIs by their correlation profiles, compared with the Kendall clustering ---
# Two profiles of every PPI (ppi_profiles.R), each a vector of signed Spearman correlations:
#   parameters          with the true parameters of the reaction norm (L, a and the shape
#                       parameters, one value per shape and parameter; not the replicate
#                       stability rho)
#   summary statistics  with five summary statistics of the reaction norm (mean, min, max,
#                       range, slope), over all genotypes pooled (as Supplementary figure S2)
# Both are clustered and compared with each other and with the Kendall clustering of the
# PPIs (Mantel test between the distance matrices, adjusted Rand index between the groups).
# Only for the PPIs of the replicate-averaged reaction norms.
if (REPLICATE_HANDLING == "mean_before") {
  wide <- ppi_wide_table(traits, shapes = unique(res$Shape))
  df_k <- wide[, setdiff(names(wide), c("Shape", "Genotype")), drop = FALSE]
  df_k <- df_k[, colSums(is.na(df_k)) == 0, drop = FALSE]
  kend <- kendall_ppi_clustering(df_k, split(seq_len(nrow(wide)), wide$Shape))
  kend_col <- group_colours(kend$hc, kend$groups)

  pr <- res[res$Parameter != "rho" & is.finite(res$rho_s), ]
  pr$Property <- paste(pr$Shape, pr$Parameter, sep = ": ")
  P_par <- with(pr, tapply(rho_s, list(PPI, Property), identity))
  ppis <- intersect(rownames(P_par), names(df_k))
  prof_par <- profile_clustering(P_par[ppis, , drop = FALSE])
  prof_stat <- profile_clustering(summary_stat_profile(population_score_table(RESOLUTION_DIR), ppis))

  cmp <- list(parameters_vs_kendall = compare_ppi_clusterings(prof_par, kend),
              statistics_vs_kendall = compare_ppi_clusterings(prof_stat, kend),
              parameters_vs_statistics = compare_ppi_clusterings(prof_par, prof_stat))
  cmp_tab <- data.frame(comparison = names(cmp),
                        mantel_r = vapply(cmp, `[[`, numeric(1), "mantel_r"),
                        mantel_p = vapply(cmp, `[[`, numeric(1), "mantel_p"),
                        ARI = vapply(cmp, `[[`, numeric(1), "ari"))
  write.csv(cmp_tab, file.path(out_dir, "ppi_profile_comparison.csv"), row.names = FALSE)
  write.csv(data.frame(PPI = ppis, parameter_group = prof_par$groups[ppis],
                       statistic_group = prof_stat$groups[ppis], kendall_group = kend$groups[ppis]),
            file.path(out_dir, "ppi_profile_groups.csv"), row.names = FALSE)

  p_prof <- patchwork::wrap_plots(
    plot_profile_clustering(prof_par, kend_col, cmp$parameters_vs_kendall, "(a) Correlation with the true parameters"),
    plot_profile_clustering(prof_stat, kend_col, cmp$statistics_vs_kendall, "(b) Correlation with summary statistics"),
    nrow = 1) +
    patchwork::plot_annotation(
      title = "PPIs grouped by their correlation profiles",
      subtitle = sprintf("(a) vs (b): Mantel r = %.2f, p = %.3f; adjusted Rand index = %.2f. Mantel tests: rank correlation of the distances, 999 permutations.",
                         cmp$parameters_vs_statistics$mantel_r, cmp$parameters_vs_statistics$mantel_p,
                         cmp$parameters_vs_statistics$ari),
      theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 11),
                             plot.subtitle = ggplot2::element_text(size = 8)))
  p_prof <- with_caption(p_prof,
    paste(setup, "Clustering of the PPIs on their signed Spearman correlations (Euclidean distance, average linkage, number of groups by silhouette width) with (a) the true parameters of every shape and (b) five summary statistics of the reaction norms. Name colours: profile groups; squares: groups of the Kendall clustering of the PPIs (all shapes pooled)."),
    "PPIs in one profile group measure the same properties of the reaction norm; agreement with the Kendall groups shows that PPIs rank genotypes alike because of that.")
  grDevices::pdf(file.path(plots_dir, "ppi_profile_clustering.pdf"), width = 9, height = 6)
  print(p_prof)
  invisible(grDevices::dev.off())
}
message("Parameter-PPI correlations written to ", out_dir, " and ", plots_dir)
