# =============================================================================
# resolution_analysis.R — PPIs of the same reaction norms at different sampling
#                         resolutions, and how their relation to the summary
#                         statistics depends on the number of sampling points
# =============================================================================
# Every resolution is a folder n_<k> of the population: the same true reaction norms
# measured at k points (sample_resolution(); same genotypes and plants), so differences
# between resolutions come from the sampling only.

summary_stat_names <- c("min", "max", "mean", "median", "slope", "range", "variance", "mean_lower", "mean_upper")

# Nine summary statistics of every reaction norm (rows of M, genotypes x environments):
# slope = slope of the linear regression of the trait on the environment values;
# mean_lower / mean_upper = mean of the lower / upper half of the sampling points
summary_statistics <- function(M, environments) {
  n <- ncol(M)
  if (n < 2) stop("summary statistics need at least 2 sampling points")
  half <- floor(n / 2)
  data.frame(
    Genotype   = rownames(M),
    min        = apply(M, 1, min, na.rm = TRUE),
    max        = apply(M, 1, max, na.rm = TRUE),
    mean       = apply(M, 1, mean, na.rm = TRUE),
    median     = apply(M, 1, median, na.rm = TRUE),
    slope      = apply(M, 1, function(y) unname(coef(lm(y ~ environments))[2])),
    range      = apply(M, 1, function(y) max(y, na.rm = TRUE) - min(y, na.rm = TRUE)),
    variance   = apply(M, 1, var, na.rm = TRUE),
    mean_lower = apply(M, 1, function(y) mean(y[seq_len(half)], na.rm = TRUE)),
    mean_upper = apply(M, 1, function(y) mean(y[(half + 1):n], na.rm = TRUE)),
    row.names = NULL, stringsAsFactors = FALSE
  )
}

# Sampling designs as in the earlier analyses: the points start, start + interval, ...
# and end of a grid of grid_length points (end is always included)
grid_design <- function(start, end, interval) unique(c(seq(start, end, by = interval), end))

# Grid indices -> points t in [0, 1] of the gradient (environment 1 + 9 t)
design_points <- function(indices, grid_length) (indices - 1) / (grid_length - 1)

# All resolutions of the simulated population in dir at once: for every n in n_points
# the resolution folder n_<n> is sampled (sample_resolution()), its PPIs computed
# (compute_population_ppis()) and its summary statistics saved
# (n_<n>/analyses/summary_statistics.csv) - exactly as when done one by one. A
# resolution that already has its PPIs is kept. Returns the resolution folders.
compute_resolutions <- function(dir, n_points, registry = default_ppi_registry()) {
  vapply(n_points, function(n) {
    res <- resolution_path(dir, n)
    if (!file.exists(file.path(ppi_set_path(res, "mean_before"), "ppis", "ppi_scores.rds"))) {
      message("Resolution: ", n, " sampling points")
      sample_resolution(dir, n)
      compute_population_ppis(res, "mean_before", registry = registry)
    }
    pop <- read_population(res)
    stats <- do.call(rbind, lapply(names(pop$phenotypes), function(s) {
      cbind(Shape = s, summary_statistics(replicate_means(pop$phenotypes[[s]], trait_environments(pop$environments, s)),
                                           trait_environments(pop$environments, s)),
            stringsAsFactors = FALSE)
    }))
    write.csv(stats, file.path(population_subdir(res, "analyses"), "summary_statistics.csv"), row.names = FALSE)
    res
  }, character(1))
}

# The score tables (population_score_table()) of the resolution folders of a population
# that have PPIs, optionally only those with n_points; one table per resolution
# ("n<points>") with the column N_points
resolution_score_tables <- function(dir, n_points = NULL) {
  res <- resolution_dirs(dir)
  res <- res[file.exists(file.path(res, "mean_before", "ppis", "ppi_scores.rds"))]
  n <- as.integer(sub("^n_", "", basename(res)))
  if (!is.null(n_points)) { keep <- n %in% n_points; res <- res[keep]; n <- n[keep] }
  if (!length(res)) stop("no resolution with PPIs in ", dir, " (run regression_sumstats_dataprep.R first)")
  setNames(lapply(seq_along(res), function(k) {
    tab <- population_score_table(res[k])
    cbind(tab[, c("Shape", "Genotype")], N_points = n[k], tab[, setdiff(names(tab), c("Shape", "Genotype"))])
  }), paste0("n", n))
}

# A score table in the layout of the earlier regression-data tables read by the
# analyses and figures: Genotype (<shape>_<genotype>), the scores, the summary statistics
as_regression_data <- function(tab) {
  cbind(data.frame(Genotype = paste(tab$Shape, tab$Genotype, sep = "_"), stringsAsFactors = FALSE),
        tab[, c(ppi_columns(tab), summary_stat_names)])
}

ppi_columns <- function(tab) setdiff(names(tab), c("Shape", "Genotype", "N_points", summary_stat_names))

# Pearson correlation across genotypes between every PPI and every summary statistic,
# for every resolution
correlations_by_resolution <- function(tables) {
  do.call(rbind, lapply(tables, function(tab) {
    grid <- expand.grid(PPI = ppi_columns(tab), Stat = summary_stat_names, stringsAsFactors = FALSE)
    grid$r <- mapply(function(p, s) suppressWarnings(cor(tab[[p]], tab[[s]], use = "pairwise.complete.obs")),
                     grid$PPI, grid$Stat)
    cbind(data.frame(N_points = tab$N_points[1]), grid)
  }))
}

# For every PPI x summary statistic: linear regression of the correlation on the number
# of sampling points (needs at least 3 resolutions with a correlation)
correlation_regression <- function(cors) {
  keys <- unique(cors[, c("PPI", "Stat")])
  rows <- lapply(seq_len(nrow(keys)), function(i) {
    d <- cors[cors$PPI == keys$PPI[i] & cors$Stat == keys$Stat[i] & !is.na(cors$r), ]
    out <- data.frame(PPI = keys$PPI[i], Stat = keys$Stat[i], intercept = NA_real_, slope = NA_real_,
                      se_slope = NA_real_, t_stat = NA_real_, p_value = NA_real_, nobs = nrow(d),
                      stringsAsFactors = FALSE)
    if (nrow(d) >= 3) {
      cf <- summary(lm(r ~ N_points, data = d))$coefficients
      out$intercept <- cf["(Intercept)", "Estimate"]
      out$slope <- cf["N_points", "Estimate"]
      out$se_slope <- cf["N_points", "Std. Error"]
      out$t_stat <- out$slope / out$se_slope
      out$p_value <- 2 * pt(-abs(out$t_stat), df = nrow(d) - 2)
    }
    out
  })
  do.call(rbind, rows)
}

# For every PPI: linear mixed model PPI ~ N_points + summary statistics + (1 | genotype)
# over all resolutions. The N_points estimate is the change of the PPI per additional
# sampling point after the summary statistics; p values BH-adjusted over the PPIs.
# A PPI is "robust" if |estimate| < 0.001 or p.adj >= 0.05, else "sample-size-sensitive".
sample_size_mixed_model <- function(tables) {
  all <- do.call(rbind, tables)
  ppis <- ppi_columns(all)
  all$GenotypeKey <- paste(all$Shape, all$Genotype)
  fm_rhs <- paste("N_points +", paste(summary_stat_names, collapse = " + "), "+ (1 | GenotypeKey)")
  rows <- lapply(ppis, function(p) {
    d <- all[is.finite(all[[p]]), c(p, "N_points", summary_stat_names, "GenotypeKey")]
    names(d)[1] <- "value"
    est <- c(estimate = NA_real_, std.error = NA_real_, t.value = NA_real_, p.value = NA_real_)
    if (length(unique(d$N_points)) >= 2) {
      fit <- tryCatch(suppressMessages(suppressWarnings(lmerTest::lmer(as.formula(paste("value ~", fm_rhs)), data = d))),
                      error = function(e) NULL)
      if (!is.null(fit)) {
        cf <- summary(fit)$coefficients
        est <- c(estimate = cf["N_points", "Estimate"], std.error = cf["N_points", "Std. Error"],
                 t.value = cf["N_points", "t value"], p.value = cf["N_points", "Pr(>|t|)"])
      }
    }
    data.frame(PPI = p, t(est), stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  out$p.adj <- p.adjust(out$p.value, method = "BH")
  out$ci_lower <- out$estimate - qnorm(0.975) * out$std.error
  out$ci_upper <- out$estimate + qnorm(0.975) * out$std.error
  out$class <- ifelse(is.na(out$estimate), NA_character_,
                      ifelse(abs(out$estimate) < 0.001 | out$p.adj >= 0.05, "robust", "sample-size-sensitive"))
  out
}

# The PPIs of a resolution folder (<replicate handling>/ppis/ppi_scores.rds of
# compute_population_ppis(); one value per genotype, so "mean_before") and the nine
# summary statistics of its replicate-mean reaction norms: one row per genotype
# (Shape, Genotype, PPIs..., statistics...)
population_score_table <- function(dir, replicate_handling = "mean_before") {
  ppi_file <- file.path(ppi_set_path(dir, replicate_handling), "ppis", "ppi_scores.rds")
  if (!file.exists(ppi_file)) stop("no PPIs in ", ppi_set_path(dir, replicate_handling), " (run compute_population_ppis(dir) first)")
  pop <- read_population(dir, with_snps = FALSE)
  wide <- ppi_wide_table(readRDS(ppi_file), shapes = names(pop$phenotypes))
  stats_file <- file.path(dir, "analyses", "summary_statistics.csv")
  stats <- if (file.exists(stats_file)) {
    read.csv(stats_file, stringsAsFactors = FALSE, colClasses = c(Genotype = "character"))
  } else {
    do.call(rbind, lapply(names(pop$phenotypes), function(s) {
      cbind(Shape = s, summary_statistics(replicate_means(pop$phenotypes[[s]], trait_environments(pop$environments, s)),
                                           trait_environments(pop$environments, s)),
            stringsAsFactors = FALSE)
    }))
  }
  cbind(wide, stats[match(paste(wide$Shape, wide$Genotype), paste(stats$Shape, stats$Genotype)), summary_stat_names])
}

# Inputs of the analyses that relate the PPIs of one resolution to its reaction norms
# (summary_stats_vs_scores.R, clustering_rn_summary_stats.R): the PPIs
# (compute_population_ppis()), the replicate-mean reaction norm of every genotype per
# shape (genotypes x environments) and the shapes
score_analysis_inputs <- function(res_dir, replicate_handling = "mean_before") {
  ppi_file <- file.path(ppi_set_path(res_dir, replicate_handling), "ppis", "ppi_scores.rds")
  if (!file.exists(ppi_file)) stop("no PPIs in ", res_dir, " (run compute_population_ppis() first)")
  pop <- read_population(res_dir, with_snps = FALSE)
  list(
    scores_list    = readRDS(ppi_file),
    rn_matrix_list = lapply(setNames(names(pop$phenotypes), names(pop$phenotypes)), function(s) {
      replicate_means(pop$phenotypes[[s]], trait_environments(pop$environments, s))
    }),
    genotype_forms = names(pop$phenotypes)
  )
}
