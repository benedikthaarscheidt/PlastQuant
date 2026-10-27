# =============================================================================
# correlation_resolution_experiment.R — resolution robustness with confounders
# =============================================================================
# WHAT IT DOES: A resolution-robustness analysis like correlation_resolution.R, but with
#   the summary statistics included in the regression as confounders.
# REQUIRES:     The resolution folders POPULATION_DIR/n_2 ... n_50 with their PPIs
#               (regression_sumstats_dataprep.R).
# PRODUCES:     confoundermodel_sumstats.pdf in POPULATION_DIR/plots/.
# HOW TO RUN:   setwd("~/PlastQuant")
#               POPULATION_DIR <- here::here("output", "populations", "ppi_comparison")
#               source(here::here("R", "analysis", "correlation_resolution_experiment.R"))
# -----------------------------------------------------------------------------
# PARAMETERS (edit the named assignment near the top of the script)
#   POPULATION_DIR  path, REQUIRED (set before sourcing)   population folder                [COMMON]
#   files       character vector  the input CSVs
#   n_samples   integer vector, default c(2,3,4,5,6,11,26,50)   sample counts
#   stat_cols   character vector  summary-stat columns to use
# =============================================================================
options(warn = -1)  # silence warnings even if the project .Rprofile was not loaded

# here we do a similar robustness analysis as in correlation_resolution.R, however in the regression we additionally consider confounders (the summary statistics) 
#in order to distinguish effects which occur due to the sampling effect on the summary statistics and might unknowingly effect the slope of the sampling resolution

library(tidyverse)
library(lme4)
library(lmerTest)
library(dplyr)   # must come last — re-asserts dplyr::select/count over MASS


if (!exists("POPULATION_DIR")) stop("set POPULATION_DIR (a population folder with the resolutions of regression_sumstats_dataprep.R) first")
source(here::here("R", "pipeline", "load_pipeline.R"))
# PPIs and summary statistics of the resolution folders n_2 ... n_50 (one table per
# resolution in the layout Genotype, scores, statistics; same order as n_samples)
regression_tables <- lapply(resolution_score_tables(POPULATION_DIR, c(2, 3, 4, 5, 6, 11, 26, 50)), as_regression_data)
files <- names(regression_tables)
n_samples <- c(2, 3, 4, 5, 6, 11, 26, 50)


stat_cols    <- c("min","max","mean","median","slope","range","variance","mean_lower","mean_upper")
first_dat    <- tibble::as_tibble(regression_tables[[files[1]]])
score_cols   <- setdiff(colnames(first_dat), c(colnames(first_dat)[1], stat_cols))
summary_cols <- stat_cols


long <- tibble(
  variable    = character(),
  value       = numeric(),
  sample_size = numeric(),
  replicate   = integer()
)

for (i in seq_along(files)) {
  dat <- tibble::as_tibble(regression_tables[[files[i]]])
  scores_df <- dat %>% dplyr::select(all_of(score_cols))
  n   <- nrow(scores_df)
  nv  <- ncol(scores_df)
  
  long <- bind_rows(long,
                    tibble(
                      variable    = rep(score_cols, each = n),
                      value       = as.vector(as.matrix(scores_df)),
                      sample_size = n_samples[i],
                      replicate   = rep(seq_len(n), times = nv)
                    )
  )
}

long <- long %>%
  mutate(
    replicate   = factor(replicate),
    sample_size = as.numeric(sample_size)
  )


wide_summ <- tibble()

for (i in seq_along(files)) {
  dat <- tibble::as_tibble(regression_tables[[files[i]]])
  
  sums_df <- dat %>%
    dplyr::select(all_of(summary_cols)) %>%
    mutate(
      replicate   = seq_len(n()),
      sample_size = n_samples[i]
    )
  
  wide_summ <- bind_rows(wide_summ, sums_df)
}

wide_summ <- wide_summ %>%
  mutate(
    replicate   = factor(replicate),
    sample_size = as.numeric(sample_size)
  )


long2 <- left_join(long, wide_summ, by = c("replicate", "sample_size"))


vars <- unique(long2$variable)
out  <- tibble(
  variable   = vars,
  estimate   = NA_real_,
  std.error  = NA_real_,
  t.value    = NA_real_,
  p.value    = NA_real_
)


conf_terms <- paste(summary_cols, collapse = " + ")

for (j in seq_along(vars)) {
  vv  <- vars[j]
  sub <- filter(long2, variable == vv, !is.na(value))
  if (length(unique(sub$sample_size)) < 2) next
  
  # build formula: value ~ sample_size + <all summary_cols> + (1|replicate)
  fm <- as.formula(
    paste0("value ~ sample_size + ", conf_terms, " + (1 | replicate)")
  )
  
  mod <- lmer(fm, data = sub)
  cf  <- summary(mod)$coefficients
  
  out$estimate[j]  <- cf["sample_size", "Estimate"]
  out$std.error[j] <- cf["sample_size", "Std. Error"]
  out$t.value[j]   <- cf["sample_size", "t value"]
  out$p.value[j]   <- cf["sample_size", "Pr(>|t|)"]
}


out <- out %>%
  mutate(
    p.adj    = p.adjust(p.value, method = "BH"),
    ci_lower = estimate - qnorm(0.975) * std.error,
    ci_upper = estimate + qnorm(0.975) * std.error,
    class    = case_when(
      abs(estimate) < 0.001        ~ "robust",
      p.adj         < 0.05         ~ "sample-size-sensitive",
      TRUE                          ~ "robust"
    )
  )

print(out)

plot_tbl <- arrange(out, estimate) %>%
  mutate(variable = factor(variable, levels = variable))

p_forest <- ggplot(plot_tbl, aes(x = estimate, y = variable, colour = class)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_pointrange(aes(xmin = ci_lower, xmax = ci_upper), size = 0.5) +
  scale_colour_manual(
    values = c("robust" = "steelblue", "sample-size-sensitive" = "tomato"),
    name   = NULL
  ) +
  labs(
    x = "Adjusted slope for sample size",
    y = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(
    legend.position  = "bottom",
    panel.grid.minor = element_blank()
  )

p_forest <- with_caption(p_forest,
  paste(resolutions_setup(POPULATION_DIR, n_samples), "All shapes pooled.",
        "For every PPI: effect of n in the mixed model PPI ~ n + nine summary statistics + (1 | genotype) over all resolutions, with 95 % CI; sample-size-sensitive (red): BH-adjusted p < 0.05 and |estimate| >= 0.001, otherwise robust (blue)."),
  "Sensitive PPIs change their values with the number of points even after accounting for the reaction norm statistics. The estimate is in the units of each PPI, so the 0.001 threshold is not comparable between PPIs; PPIs computed on the index of the measurement point change with n by construction.")
print(p_forest)

path <- file.path(POPULATION_DIR, "plots")
if (!dir.exists(path)) dir.create(path, recursive = TRUE)

ggsave(create.dir = TRUE, 
  filename = file.path(path, "confoundermodel_sumstats.pdf"),
  plot     = p_forest,
  width    = 6.3,
  height   = 6,
  units    = "in",
  dpi      = 900,
  device   = "pdf"
)
