# =============================================================================
# 05_evaluate_gwas.R — GWAS accuracy evaluation, per reaction-norm shape
# =============================================================================
# WHAT IT DOES: Compares the GWAS hits of one run against the simulated causal SNPs.
#   Every shape is evaluated against its own causal SNPs only, so the power of a
#   parameter is always the power within one shape (a hit on a causal SNP of another
#   shape is a false positive). Draws one power heatmap per shape plus an overview.
# REQUIRES:     A run folder written by run_simulation_scenario() / run_pipeline()
#               (gwas_results.csv, truth.csv, config.rds), e.g. by a scenario_*.R driver.
# PRODUCES:     RUN_DIR/gwas_accuracy_summary.csv and RUN_DIR/plots/gwas_power_heatmap.pdf.
# HOW TO RUN:   setwd("~/PlastQuant")
#               OUTPUT_BASE <- here::here("output", "scenario_20_5")   # newest run is used
#               source(here::here("R", "05_evaluate_gwas.R"))
# -----------------------------------------------------------------------------
# PARAMETERS — set these before sourcing (guarded defaults)
#   RUN_DIR      path, default newest run_* in OUTPUT_BASE    run folder to evaluate  [COMMON]
#   OUTPUT_BASE  path                                         scenario folder (if no RUN_DIR)
#   ALPHA        numeric, default 0.05                        Bonferroni level: P < ALPHA / n_SNPs
# =============================================================================
options(warn = -1)  # silence warnings even if the project .Rprofile was not loaded
library(ggplot2)
source(here::here("R", "pipeline", "load_pipeline.R"))

if (!exists("RUN_DIR")) {
  if (!exists("OUTPUT_BASE")) stop("set RUN_DIR (a run folder) or OUTPUT_BASE (a scenario folder) first")
  RUN_DIR <- latest_run_dir(OUTPUT_BASE)
}
if (!exists("ALPHA")) ALPHA <- 0.05

# -----------------------------------------------------------------------------
# 1. Load the run
# -----------------------------------------------------------------------------
message("Evaluating GWAS of: ", RUN_DIR)
gwas_res <- read.csv(file.path(RUN_DIR, "gwas_results.csv"), stringsAsFactors = FALSE)
truth    <- read.csv(file.path(RUN_DIR, "truth.csv"), stringsAsFactors = FALSE)
cfg      <- readRDS(file.path(RUN_DIR, "config.rds"))
n_snps   <- length(unique(gwas_res$SNP))

# -----------------------------------------------------------------------------
# 2. Power per shape, parameter and PPI (Bonferroni)
# -----------------------------------------------------------------------------
power <- evaluate_power(gwas_res, truth, n_snps, alpha = ALPHA)
power <- power[, c("Shape", "Score", "Parameter", "True_SNPs", "Found", "Power")]
write.csv(power, file.path(RUN_DIR, "gwas_accuracy_summary.csv"), row.names = FALSE)

subtitle_base <- paste0(
  "h2 = ", cfg$h2, ", ", cfg$Mc, " causal SNPs per parameter, ",
  cfg$n_per_shape, " genotypes per shape, ", n_snps, " SNPs",
  "  (Bonferroni P < ", format(ALPHA / n_snps, digits = 2, scientific = TRUE), ")"
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

plots_dir <- file.path(RUN_DIR, "plots")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
pdf(file.path(plots_dir, "gwas_power_heatmap.pdf"), width = 10, height = 8)
print(make_heatmap(power, "GWAS power per shape", subtitle_base) +
        facet_grid(. ~ Shape, scales = "free_x", space = "free_x"))
for (s in unique(power$Shape)) {
  fp <- power[power$Shape == s & power$Parameter == "False_Positives", ]
  print(make_heatmap(power[power$Shape == s, ], paste0("GWAS power: ", s),
                     paste0(subtitle_base, "\nFalse positives (all PPIs): ", sum(fp$Found))))
}
dev.off()

message("Summary CSV : ", file.path(RUN_DIR, "gwas_accuracy_summary.csv"))
message("Heatmap PDF : ", file.path(plots_dir, "gwas_power_heatmap.pdf"))
