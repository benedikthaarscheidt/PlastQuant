# =============================================================================
# 05_evaluate_gwas.R — GWAS accuracy evaluation, per reaction-norm shape
# =============================================================================
# WHAT IT DOES: Compares the GWAS hits of one resolution of a population against the
#   population's causal SNPs.
#   Every shape is evaluated against its own causal SNPs only, so the power of a
#   parameter is always the power within one shape (a hit on a causal SNP of another
#   shape is a false positive). Draws one power heatmap per shape plus an overview.
# REQUIRES:     A resolution folder with GWAS results (<REPLICATE_HANDLING>/gwas/gwas_results.csv) of a
#               population with a truth table (<population>/data/truth.csv), e.g. from a scenario_*.R driver
#               or run_population_gwas() on own data with known causal SNPs.
# PRODUCES:     RESOLUTION_DIR/<REPLICATE_HANDLING>/gwas/gwas_accuracy_summary.csv and
#               RESOLUTION_DIR/<REPLICATE_HANDLING>/plots/gwas_power_heatmap.pdf.
# HOW TO RUN:   setwd("~/PlastQuant")
#               RESOLUTION_DIR <- here::here("output", "populations", "scenario_20_5", "n_6")
#               source(here::here("R", "gwas", "05_evaluate_gwas.R"))
# -----------------------------------------------------------------------------
# PARAMETERS — set these before sourcing (guarded defaults)
#   RESOLUTION_DIR  path, REQUIRED              resolution folder to evaluate         [COMMON]
#   REPLICATE_HANDLING  "mean_before" (default) or "replicates"   which PPIs of the resolution
#   ALPHA           numeric, default 0.05       Bonferroni level: P < ALPHA / n_SNPs
# =============================================================================
options(warn = -1)  # silence warnings even if the project .Rprofile was not loaded
library(ggplot2)
source(here::here("R", "pipeline", "load_pipeline.R"))

if (!exists("RESOLUTION_DIR")) stop("set RESOLUTION_DIR (a resolution folder, e.g. output/populations/<name>/n_6) first")
if (!exists("REPLICATE_HANDLING")) REPLICATE_HANDLING <- "mean_before"   # PPIs of mean_before/ or replicates/
PPI_SET_DIR <- file.path(RESOLUTION_DIR, REPLICATE_HANDLING)
if (!exists("ALPHA")) ALPHA <- 0.05

# -----------------------------------------------------------------------------
# 1. Load the population's GWAS results and truth
# -----------------------------------------------------------------------------
message("Evaluating GWAS of: ", RESOLUTION_DIR)
gwas_res <- read.csv(file.path(PPI_SET_DIR, "gwas", "gwas_results.csv"), stringsAsFactors = FALSE)
population_dir <- dirname(normalizePath(RESOLUTION_DIR))
truth    <- read.csv(file.path(population_dir, "data", "truth.csv"), stringsAsFactors = FALSE)
names(truth)[names(truth) == "Trait"] <- "Shape"
n_snps   <- length(unique(gwas_res$SNP))
# Simulation settings exist only for simulated runs (own data has no config.rds)
cfg_file <- file.path(population_dir, "data", "config.rds")
cfg      <- if (file.exists(cfg_file)) readRDS(cfg_file) else NULL

# -----------------------------------------------------------------------------
# 2. Power per shape, parameter and PPI (Bonferroni)
# -----------------------------------------------------------------------------
power <- evaluate_power(gwas_res, truth, n_snps, alpha = ALPHA)
power <- power[, c("Shape", "Score", "Parameter", "True_SNPs", "Found", "Power")]
write.csv(power, file.path(PPI_SET_DIR, "gwas", "gwas_accuracy_summary.csv"), row.names = FALSE)

subtitle_base <- paste0(
  if (!is.null(cfg)) paste0("h2 = ", cfg$h2, ", ", cfg$Mc, " causal SNPs per parameter, ",
                            cfg$n_per_shape, " genotypes per shape, "),
  n_snps, " SNPs  (Bonferroni P < ", format(ALPHA / n_snps, digits = 2, scientific = TRUE), ")"
)

# -----------------------------------------------------------------------------
# 3. Heatmaps: overview of all shapes, then one page per shape
# -----------------------------------------------------------------------------
make_heatmap <- function(df, title, subtitle) {
  df <- df[df$Parameter != "False_Positives", ]
  ggplot(df, aes(x = Parameter, y = Score, fill = Power)) +
    geom_tile(color = "white") +
    scale_fill_gradient(low = "white", high = "red", limits = c(0, 1), na.value = "grey90") +
    geom_text(aes(label = ifelse(is.na(Power), "", sprintf("%.2f", Power))), color = "black", size = 3) +
    labs(title = title, subtitle = subtitle, x = "Causal parameter (truth of this shape)",
         y = "Plasticity index") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

plots_dir <- population_subdir(PPI_SET_DIR, "plots")
pdf(file.path(plots_dir, "gwas_power_heatmap.pdf"), width = 10, height = 8)
print(make_heatmap(power, "GWAS power per shape", subtitle_base) +
        facet_grid(. ~ Shape, scales = "free_x", space = "free_x"))
for (s in unique(power$Shape)) {
  fp <- power[power$Shape == s & power$Parameter == "False_Positives", ]
  print(make_heatmap(power[power$Shape == s, ], paste0("GWAS power: ", s),
                     paste0(subtitle_base, "\nFalse positives (all PPIs): ", sum(fp$Found))))
}
dev.off()

message("Summary CSV : ", file.path(PPI_SET_DIR, "gwas", "gwas_accuracy_summary.csv"))
message("Heatmap PDF : ", file.path(plots_dir, "gwas_power_heatmap.pdf"))
