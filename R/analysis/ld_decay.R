# =============================================================================
# ld_decay.R — LD between the simulated SNPs against their distance
# =============================================================================
# WHAT IT DOES: Mean r^2 between the dosages of SNPs d apart, for d = 1..MAX_DISTANCE,
#   from the SNP matrix of a population, with the LD window of the simulation (the
#   largest distance with mean r^2 >= ld_r2) marked. Shows how strong ld_strength is.
# REQUIRES:     a simulated population (data/snp_matrix.rds, data/config.rds).
# PRODUCES:     POPULATION_DIR/plots/ld_decay.pdf
#               POPULATION_DIR/analyses/ld_decay.csv
# HOW TO RUN:   setwd("~/PlastQuant")
#               POPULATION_DIR <- here::here("output", "populations", "My_population")
#               source(here::here("R", "analysis", "ld_decay.R"))
# -----------------------------------------------------------------------------
# PARAMETERS — set these before sourcing (guarded defaults)
#   POPULATION_DIR  path, REQUIRED   population folder
#   MAX_DISTANCE    integer (50)     largest distance in SNPs
# =============================================================================
source(here::here("R", "pipeline", "load_pipeline.R"))

if (!exists("POPULATION_DIR")) stop("set POPULATION_DIR (a population folder) first")
if (!exists("MAX_DISTANCE")) MAX_DISTANCE <- 50

X <- readRDS(file.path(POPULATION_DIR, "data", "snp_matrix.rds"))
cfg <- readRDS(file.path(POPULATION_DIR, "data", "config.rds"))
ld_strength <- if (is.null(cfg$ld_strength)) 0 else cfg$ld_strength
r2 <- if (is.null(cfg$ld_r2)) 0.1 else cfg$ld_r2
window <- if (is.null(cfg$ld_window)) 0 else cfg$ld_window

decay <- ld_decay(X, seq_len(min(MAX_DISTANCE, ncol(X) - 2)))
write.csv(decay, file.path(population_subdir(POPULATION_DIR, "analyses"), "ld_decay.csv"), row.names = FALSE)

p <- ggplot2::ggplot(decay, ggplot2::aes(x = Distance, y = R2)) +
  ggplot2::geom_line() +
  ggplot2::geom_point(size = 1) +
  ggplot2::geom_hline(yintercept = r2, linetype = "dashed") +
  ggplot2::geom_vline(xintercept = window, linetype = "dotted") +
  ggplot2::labs(title = "LD decay of the simulated SNPs",
                subtitle = paste0(basename(POPULATION_DIR), ", ld_strength = ", ld_strength,
                                  ", LD window = ", window, " SNPs (mean r^2 >= ", r2, ")"),
                x = "Distance (SNPs)", y = "Mean r^2") +
  ggplot2::theme_minimal(base_size = 10)
grDevices::pdf(file.path(population_subdir(POPULATION_DIR, "plots"), "ld_decay.pdf"), width = 7, height = 4.5)
print(p)
invisible(grDevices::dev.off())
message("LD decay written to ", file.path(POPULATION_DIR, "plots", "ld_decay.pdf"))
