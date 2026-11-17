# =============================================================================
# scenario_20_5.R — heritability scenario: h2 = 0.2, 5 causal SNPs per parameter
# =============================================================================
# Simulates the population (R/sim) into output/populations/scenario_20_5/ (data/,
# variance/), measures it at 6 points (n_6/), computes all PPIs from the replicate-averaged
# reaction norms (n_6/mean_before/ppis/) and runs the GWAS with the power per shape and
# its heatmaps (n_6/mean_before/gwas/, n_6/mean_before/plots/).
# All other settings (variance components, study space, ...) are the defaults of
# sim_config() in R/sim/simulation_config.R.
# HOW TO RUN: Rscript --vanilla R/gwas/scenarios/scenario_20_5.R   (from the repo root)
# -----------------------------------------------------------------------------
# PARAMETERS — set these before sourcing (guarded defaults)
#   H2  the heritability setting, default 0.6
#   M_C  the number of causal SNPs, default 5
#   MAF  the MAF to prefilter the original SNPs, default 0.05
# =============================================================================
options(warn = -1)  # silence warnings even if the project .Rprofile was not loaded

if (!exists("H2")) H2 <- 0.6
if (!exists("M_C")) M_C <- 5
if (!exists("MAF")) MAF <- 0.05
source(here::here("R", "pipeline", "load_pipeline.R"))
geno <- new.env()
load(here::here("data", "A-Datasets", "Genotyping_matrix_Negro_et_al_2019.Rdata"), envir = geno)
map <- readRDS(here::here("data", "A-Datasets", "SNPs_DL_windows_Negro_et_al_2019.rds"))
map <- map[map$SNP.name %in% colnames(geno$mat2), ]
map <- map[order(map$Chromosome, map$Position), ]            # genomic order
snps <- prune_snps(geno$mat2[, map$SNP.name], r2 = 0.1, window = 50, chrom = map$Chromosome, MAF = MAF)
#snps <- snps[1:200,]
rm(geno, map)
scenario_name = paste("scenario", as.character(H2*100), as.character(M_C), as.character(MAF*100),sep = "_")
POPULATION_DIR <- here::here("output", "populations", scenario_name)

#simulate_population(sim_config(h2 = 0.2, Mc = 5), POPULATION_DIR)   # true reaction norms -> data/, variance/
simulate_population_from_snps(sim_config(h2 = H2, Mc = M_C), snps, POPULATION_DIR)   # -> data/, variance/RESOLUTION_DIR <- sample_resolution(POPULATION_DIR, 6)    # measured at 6 points -> n_6/data/
RESOLUTION_DIR <- sample_resolution(POPULATION_DIR, 6)    # measured at 6 points -> n_6/data/                                         # -> n_6/data/
compute_population_ppis(RESOLUTION_DIR)                   # PPIs -> n_6/mean_before/ppis/
run_population_gwas(RESOLUTION_DIR)                       # GWAS, power, heatmaps -> n_6/mean_before/gwas/, plots/

