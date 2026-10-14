# =============================================================================
# combined_resolution_figure.R — combined resolution figure
# =============================================================================
# WHAT IT DOES: Draws the sampling-sensitivity figure of a population: left, the slope
#   of the regression of the PPI-summary-statistic correlation on the number of
#   sampling points (heatmap, from R/analysis/resolution_analysis.R); right, the effect
#   of the number of sampling points on every PPI from the mixed model
#   PPI ~ n + summary statistics + (1 | genotype), fitted here on the resolution tables
#   (forest plot).
# REQUIRES:     POPULATION_DIR/analyses/regression_summary_stats/ (regression_sumstats_dataprep.R)
#               and POPULATION_DIR/analyses/resolution/ (resolution_analysis.R).
# PRODUCES:     POPULATION_DIR/plots/combined_resolution_figure.pdf and the mixed-model
#               results POPULATION_DIR/analyses/resolution/mixed_model.csv
# HOW TO RUN:   setwd("~/PlastQuant")
#               POPULATION_DIR <- here::here("output", "populations", "ppi_comparison")
#               source(here::here("R","figures","combined_resolution_figure.R"))
# -----------------------------------------------------------------------------
# PARAMETERS — set before sourcing
#   POPULATION_DIR  path, REQUIRED   population folder with the resolution analysis      [COMMON]
# =============================================================================
options(warn = -1)  # silence warnings even if the project .Rprofile was not loaded

library(dplyr)
library(ggplot2)
library(patchwork)
source(here::here("R", "pipeline", "load_pipeline.R"))

if (!exists("POPULATION_DIR")) stop("set POPULATION_DIR (a population folder) first")
res_dir <- file.path(POPULATION_DIR, "analyses", "resolution")
if (!file.exists(file.path(res_dir, "correlation_regression.csv"))) {
  stop("no resolution analysis in ", POPULATION_DIR, " (run R/analysis/resolution_analysis.R first)")
}
model_df <- read.csv(file.path(res_dir, "correlation_regression.csv"), stringsAsFactors = FALSE)

# Mixed model on the PPI and summary-statistic tables of all sampling designs
tables <- read_regression_tables(file.path(POPULATION_DIR, "analyses", "regression_summary_stats"), range = "full")
out <- sample_size_mixed_model(tables)
write.csv(out, file.path(res_dir, "mixed_model.csv"), row.names = FALSE)

# ── LEFT panel: slope heatmap ─────────────────────────────────────────────────

slope_heat <- model_df %>%
  filter(!is.na(slope)) %>%
  mutate(
    p_label = ifelse(is.na(p_value), "",
                ifelse(p_value < 0.001, "***",
                  ifelse(p_value < 0.01, "**",
                    ifelse(p_value < 0.05, "*", "")))),
    PS   = factor(PPI,  levels = rev(sort(unique(PPI)))),
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

# ── RIGHT panel: confounder model forest plot ─────────────────────────────────

plot_tbl <- out %>%
  filter(!is.na(estimate)) %>%
  arrange(estimate) %>%
  mutate(variable = factor(PPI, levels = PPI))

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

out_path <- population_subdir(POPULATION_DIR, "plots")
ggsave(
  filename = file.path(out_path, "combined_resolution_figure.pdf"),
  plot     = combined,
  width    = 16,
  height   = 20,
  units    = "cm",
  dpi      = 900,
  device   = "pdf"
)
message("Figure written to ", file.path(out_path, "combined_resolution_figure.pdf"))
