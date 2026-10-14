# =============================================================================
# resolution_analysis.R — PPIs of the same reaction norms at different sampling
#                         resolutions, and how their relation to the summary
#                         statistics depends on the number of sampling points
# =============================================================================
# The reaction norms of a simulated population are sampled again at every resolution
# (sample_population(): same genotypes and plants, new measurement), so differences
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

# One table for a measured population: one row per genotype (Shape, Genotype,
# N_points), one column per PPI (replicates averaged before the PPI) and per summary
# statistic (of the replicate-mean reaction norm)
resolution_table <- function(measured, registry = default_ppi_registry()) {
  traits <- compute_ppis(measured$phenotypes, measured$environments, registry)
  do.call(rbind, lapply(names(measured$phenotypes), function(shape) {
    M <- replicate_means(measured$phenotypes[[shape]], measured$environments)
    tab <- data.frame(Shape = shape, Genotype = rownames(M), N_points = length(measured$environments),
                      stringsAsFactors = FALSE)
    for (p in names(traits)) {
      sc <- traits[[p]][traits[[p]]$Type == shape, ]
      tab[[p]] <- sc$Score[match(tab$Genotype, sc$Genotype_ID)]
    }
    cbind(tab, summary_statistics(M, measured$environments)[, summary_stat_names])
  }))
}

# The simulated population in dir sampled at every number of points in n_points;
# a named list of resolution_table()s ("n<points>")
resolution_tables <- function(dir, n_points, registry = default_ppi_registry()) {
  setNames(lapply(n_points, function(n) {
    message("Resolution: ", n, " sampling points")
    resolution_table(sample_population(dir, n), registry)
  }), paste0("n", n_points))
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

# The population's PPIs (ppis/ppi_scores.rds) and the nine summary statistics of its
# replicate-mean reaction norms, at the population's own sampling resolution: one row
# per genotype (Shape, Genotype, PPIs..., statistics...)
population_score_table <- function(dir) {
  ppi_file <- file.path(dir, "ppis", "ppi_scores.rds")
  if (!file.exists(ppi_file)) stop("no PPIs in ", dir, " (run run_pipeline(dir, gwas = FALSE) first)")
  pop <- read_population(dir)
  wide <- ppi_wide_table(readRDS(ppi_file), shapes = names(pop$phenotypes))
  stats <- do.call(rbind, lapply(names(pop$phenotypes), function(s) {
    cbind(Shape = s, summary_statistics(replicate_means(pop$phenotypes[[s]], pop$environments), pop$environments),
          stringsAsFactors = FALSE)
  }))
  cbind(wide, stats[match(paste(wide$Shape, wide$Genotype), paste(stats$Shape, stats$Genotype)), summary_stat_names])
}
