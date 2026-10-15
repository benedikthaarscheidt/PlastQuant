# =============================================================================
# combine_resolution_gwas.R — GWAS power over the resolutions of one population
# =============================================================================
# WHAT IT DOES: Collects gwas_power.csv of every resolution n_<n> of a population (one
#   replicate handling) into one table and draws the power of every PPI against the
#   number of measurement points, per shape and causal parameter.
# REQUIRES:     the resolutions with their GWAS (resolution_gwas.R resolution ...)
# PRODUCES:     POPULATION_DIR/analyses/resolution_gwas/<handling>/gwas_power_by_resolution.csv
#               POPULATION_DIR/analyses/resolution_gwas/<handling>/gwas_power_by_resolution.pdf
# HOW TO RUN:   from the repo root
#   Rscript --vanilla R/gwas/scenarios/combine_resolution_gwas.R <population folder> [mean_before|replicates]
#   or in R: POPULATION_DIR <- "..."; source(here::here("R", "gwas", "scenarios", "combine_resolution_gwas.R"))
# -----------------------------------------------------------------------------
# PARAMETERS — command-line arguments, or set these before sourcing
#   POPULATION_DIR      path                             population folder           [COMMON]
#   REPLICATE_HANDLING  "mean_before" (default) or "replicates"
# =============================================================================
options(warn = -1)  # silence warnings even if the project .Rprofile was not loaded
source(here::here("R", "pipeline", "load_pipeline.R"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) POPULATION_DIR <- args[1]
if (length(args) >= 2) REPLICATE_HANDLING <- args[2]
if (!exists("POPULATION_DIR")) stop("set POPULATION_DIR (or pass the population folder as argument)")
if (!exists("REPLICATE_HANDLING")) REPLICATE_HANDLING <- "mean_before"

res_dirs <- resolution_dirs(POPULATION_DIR)
power_files <- file.path(res_dirs, REPLICATE_HANDLING, "gwas", "gwas_power.csv")
done <- file.exists(power_files)
if (!any(done)) stop("no resolution of ", POPULATION_DIR, " has gwas/gwas_power.csv for ", REPLICATE_HANDLING)
if (!all(done)) message("without GWAS (skipped): ", paste(basename(res_dirs[!done]), collapse = ", "))

power <- do.call(rbind, lapply(which(done), function(k) {
  d <- read.csv(power_files[k], stringsAsFactors = FALSE)
  cbind(N_points = as.integer(sub("^n_", "", basename(res_dirs[k]))), d)
}))
power <- power[order(power$Shape, power$Parameter, power$Score, power$N_points), ]
rownames(power) <- NULL

out_dir <- file.path(population_subdir(POPULATION_DIR, "analyses"), "resolution_gwas", REPLICATE_HANDLING)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(power, file.path(out_dir, "gwas_power_by_resolution.csv"), row.names = FALSE)

# One page per shape: power of every PPI against n, one panel per causal parameter;
# last page: false positives of every PPI against n, one panel per shape
grDevices::pdf(file.path(out_dir, "gwas_power_by_resolution.pdf"), width = 12, height = 8)
for (s in unique(power$Shape)) {
  d <- power[power$Shape == s & power$Parameter != "False_Positives", ]
  print(ggplot2::ggplot(d, ggplot2::aes(x = N_points, y = Power, colour = Score)) +
          ggplot2::geom_line() + ggplot2::geom_point(size = 1) +
          ggplot2::facet_wrap(~ Parameter) +
          ggplot2::scale_x_log10(breaks = sort(unique(power$N_points))) +
          ggplot2::coord_cartesian(ylim = c(0, 1)) +
          ggplot2::labs(title = paste0("GWAS power against the resolution: ", s),
                        subtitle = paste0(basename(POPULATION_DIR), ", ", REPLICATE_HANDLING),
                        x = "Measurement points n", y = "Power (found / causal SNPs)", colour = "PPI") +
          ggplot2::theme_minimal())
}
fp <- power[power$Parameter == "False_Positives", ]
print(ggplot2::ggplot(fp, ggplot2::aes(x = N_points, y = Found, colour = Score)) +
        ggplot2::geom_line() + ggplot2::geom_point(size = 1) +
        ggplot2::facet_wrap(~ Shape) +
        ggplot2::scale_x_log10(breaks = sort(unique(power$N_points))) +
        ggplot2::labs(title = "False positives against the resolution", x = "Measurement points n",
                      y = "Significant non-causal SNPs", colour = "PPI") +
        ggplot2::theme_minimal())
invisible(grDevices::dev.off())
message("GWAS power over ", sum(done), " resolutions written to ", out_dir)
