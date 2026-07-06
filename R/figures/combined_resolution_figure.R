# =============================================================================
# combined_resolution_figure.R — combined resolution figure
# =============================================================================
# WHAT IT DOES: Combines the resolution heatmap (left) and the confounder/summary-stat
#   panel (right) into one publication figure.
# REQUIRES:     Regression-summary CSVs under output/regression_summary_stats/
#               (produced by regression_sumstats_dataprep.R / regression_summarystats.R).
# PRODUCES:     A combined resolution figure PDF.
# HOW TO RUN:   setwd("~/PlastQuant"); source(here::here("R","figures","combined_resolution_figure.R"))
# -----------------------------------------------------------------------------
# PARAMETERS (edit the named assignment near the top of the script)
#   directory   path, default output/regression_summary_stats  input CSV folder           [COMMON]
#   files       character vector  the specific input CSVs to combine
#   n_samples   integer vector, default c(2,3,4,5,6,11,26,50)   sample counts per resolution
#   stat_cols   character vector  summary-stat columns to use
# =============================================================================

# Combines slope_heatmap_resolution (left) and confoundermodel_sumstats (right)
# into a single side-by-side figure.

library(tidyverse)
library(patchwork)
library(viridis)
library(lme4)
library(lmerTest)
library(ggplot2)
library(forcats)
library(dplyr)

# ── Shared inputs ────────────────────────────────────────────────────────────

directory <- here::here("data", "regression_summary_stats")
files <- c(
  "regression_data_full_interval_49_indices_1_50.csv",
  "regression_data_full_interval_25_indices_1_26_50.csv",
  "regression_data_full_interval_20_indices_1_21_41_50.csv",
  "regression_data_full_interval_15_indices_1_16_31_46_50.csv",
  "regression_data_full_interval_10_indices_1_11_21_31_41_50.csv",
  "regression_data_full_interval_5_indices_1_6_11_16_21_26_31_36_41_46_50.csv",
  "regression_data_full_interval_2_indices_1_3_5_7_9_11_13_15_17_19_21_23_25_27_29_31_33_35_37_39_41_43_45_47_49_50.csv",
  "regression_data_full_interval_1_indices_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27_28_29_30_31_32_33_34_35_36_37_38_39_40_41_42_43_44_45_46_47_48_49_50.csv"
)
n_samples <- c(2, 3, 4, 5, 6, 11, 26, 50)
stat_cols <- c("min", "max", "mean", "median", "slope", "range", "variance",
               "mean_lower", "mean_upper")

# ── LEFT panel: slope heatmap (from correlation_resolution.R) ────────────────

wide_list <- vector("list", length(files))
for (i in seq_along(files)) {
  df      <- read.csv(file.path(directory, files[i]), stringsAsFactors = FALSE)
  PS      <- df[, !colnames(df) %in% c(colnames(df)[1], stat_cols), drop = FALSE]
  summary <- df[, stat_cols, drop = FALSE]

  cors <- matrix(NA, nrow = ncol(PS) * ncol(summary), ncol = 3,
                 dimnames = list(NULL, c("PS", "Stat", "r")))
  idx <- 1
  for (p in colnames(PS)) {
    for (s in colnames(summary)) {
      cors[idx, "PS"]   <- p
      cors[idx, "Stat"] <- s
      cors[idx, "r"]    <- cor(PS[[p]], summary[[s]],
                                use = "pairwise.complete.obs", method = "pearson")
      idx <- idx + 1
    }
  }
  wide_list[[i]]       <- data.frame(cors, stringsAsFactors = FALSE)
  wide_list[[i]]$r     <- as.numeric(wide_list[[i]]$r)
  names(wide_list)[i]  <- paste0("n", n_samples[i])
}

wide_df <- setNames(wide_list[[1]], c("PS", "Stat", "n2"))
for (i in 2:length(wide_list)) {
  wide_df <- merge(wide_df, wide_list[[i]],
                   by = c("PS", "Stat"), all = TRUE,
                   suffixes = c("", paste0(".n", i)))
  colnames(wide_df)[ncol(wide_df)] <- paste0("n", n_samples[i])
}

nr       <- nrow(wide_df)
model_df <- data.frame(
  PS        = wide_df$PS,
  Stat      = wide_df$Stat,
  intercept = NA_real_, slope   = NA_real_, se_slope = NA_real_,
  nobs      = NA_integer_, df   = NA_integer_,
  t_stat    = NA_real_,    p_value = NA_real_,
  stringsAsFactors = FALSE
)

for (j in seq_len(nr)) {
  y   <- as.numeric(wide_df[j, paste0("n", n_samples)])
  x   <- n_samples
  ok  <- !is.na(y)
  n_j <- sum(ok)
  if (n_j >= 3) {
    fit                    <- lm(y[ok] ~ x[ok])
    cf                     <- summary(fit)$coefficients
    model_df$intercept[j]  <- cf["(Intercept)", "Estimate"]
    model_df$slope[j]      <- cf["x[ok]", "Estimate"]
    model_df$se_slope[j]   <- cf["x[ok]", "Std. Error"]
    model_df$nobs[j]       <- n_j
    model_df$df[j]         <- n_j - 2
    tval                   <- model_df$slope[j] / model_df$se_slope[j]
    model_df$t_stat[j]     <- tval
    model_df$p_value[j]    <- 2 * pt(-abs(tval), df = model_df$df[j])
  }
}

slope_heat <- model_df %>%
  filter(!is.na(slope)) %>%
  mutate(
    p_label = ifelse(is.na(p_value), "",
                ifelse(p_value < 0.001, "***",
                  ifelse(p_value < 0.01, "**",
                    ifelse(p_value < 0.05, "*", "")))),
    PS   = factor(PS,   levels = rev(sort(unique(PS)))),
    Stat = factor(Stat, levels = sort(unique(Stat)))
  )

p_left <- ggplot(slope_heat, aes(x = Stat, y = PS, fill = slope)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = p_label), size = 2.5, colour = "black") +
  scale_fill_gradient2(
    low      = "#4575B4",
    mid      = "white",
    high     = "#D73027",
    midpoint = 0,
    name     = "Slope"
  ) +
  labs(x = "Summary Statistic", y = "Plasticity Score") +
  theme_minimal(base_size = 8) +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1),
    axis.text.y     = element_text(size = 7),
    legend.position = "right",
    panel.grid      = element_blank()
  )

# ── RIGHT panel: confounder model forest plot (from correlation_resolution_experiment.R) ──

first_dat  <- read_csv(file.path(directory, files[1]), show_col_types = FALSE)
score_cols <- setdiff(colnames(first_dat), c(colnames(first_dat)[1], stat_cols))

long <- tibble(variable = character(), value = numeric(),
               sample_size = numeric(), replicate = integer())

for (i in seq_along(files)) {
  dat       <- read_csv(file.path(directory, files[i]), show_col_types = FALSE)
  scores_df <- dat %>% dplyr::select(all_of(score_cols))
  n         <- nrow(scores_df)
  nv        <- ncol(scores_df)
  long <- bind_rows(long,
    tibble(
      variable    = rep(score_cols, each = n),
      value       = as.vector(as.matrix(scores_df)),
      sample_size = n_samples[i],
      replicate   = rep(seq_len(n), times = nv)
    ))
}
long <- long %>% mutate(replicate = factor(replicate), sample_size = as.numeric(sample_size))

wide_summ <- tibble()
for (i in seq_along(files)) {
  dat <- read_csv(file.path(directory, files[i]), show_col_types = FALSE)
  wide_summ <- bind_rows(wide_summ,
    dat %>%
      dplyr::select(all_of(stat_cols)) %>%
      mutate(replicate = seq_len(n()), sample_size = n_samples[i])
  )
}
wide_summ <- wide_summ %>%
  mutate(replicate = factor(replicate), sample_size = as.numeric(sample_size))

long2 <- left_join(long, wide_summ, by = c("replicate", "sample_size"))

vars <- unique(long2$variable)
out  <- tibble(
  variable  = vars,
  estimate  = NA_real_, std.error = NA_real_,
  t.value   = NA_real_, p.value   = NA_real_
)

conf_terms <- paste(stat_cols, collapse = " + ")

for (j in seq_along(vars)) {
  vv  <- vars[j]
  sub <- filter(long2, variable == vv, !is.na(value))
  if (length(unique(sub$sample_size)) < 2) next
  fm  <- as.formula(paste0("value ~ sample_size + ", conf_terms, " + (1 | replicate)"))
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
      abs(estimate) < 0.001 ~ "robust",
      p.adj < 0.05          ~ "sample-size-sensitive",
      TRUE                  ~ "robust"
    )
  )

plot_tbl <- arrange(out, estimate) %>%
  mutate(variable = factor(variable, levels = variable))

p_right <- ggplot(plot_tbl, aes(x = estimate, y = variable, colour = class)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_pointrange(aes(xmin = ci_lower, xmax = ci_upper), size = 0.3) +
  scale_colour_manual(
    values = c("robust" = "steelblue", "sample-size-sensitive" = "tomato"),
    name   = NULL
  ) +
  labs(x = "Adjusted slope for sample size", y = NULL) +
  theme_minimal(base_size = 10) +
  theme(
    legend.position  = "bottom",
    panel.grid.minor = element_blank()
  )

# ── Combine and save ──────────────────────────────────────────────────────────

combined <- p_left + p_right +
  plot_layout(ncol = 2, widths = c(1, 1)) +
  plot_annotation(tag_levels = "A")

out_path <- here::here("output", "plots", "figures")
if (!dir.exists(out_path)) dir.create(out_path, recursive = TRUE)

ggsave(
  filename = file.path(out_path, "combined_resolution_figure.pdf"),
  plot     = combined,
  width    = 16,
  height   = 20,
  units    = "cm",
  dpi      = 900,
  device   = "pdf"
)

print(combined)
