# =============================================================================
# scenario_h2_Mc.R — heritability scenario on the maize SNP panel: h2, Mc, MAF and seeds set below
# =============================================================================
# Prunes the maize panel (Negro et al. 2019) once to tag SNPs with MAF >= MAF; then, for
# every seed in SEEDS, simulates the population on it (R/sim) into
# output/populations/scenario_<h2>_<Mc>_<MAF>_seed<seed>/ (data/, variance/), measures it
# at 6 points (n_6/), computes all PPIs from the replicate-averaged reaction norms
# (n_6/mean_before/ppis/) and runs the GWAS with the power per shape, its heatmaps and the
# MAF of the hits (n_6/mean_before/gwas/, n_6/mean_before/plots/).
# All other settings (variance components, study space, ...) are the defaults of
# sim_config() in R/sim/simulation_config.R.
# HOW TO RUN (from the repo root):
#   Rscript --vanilla R/gwas/scenarios/scenario_h2_Mc.R            # defaults below
#   Rscript --vanilla -e "H2 <- 0.4; M_C <- 10; SEEDS <- 1:10; source('R/gwas/scenarios/scenario_h2_Mc.R')"
# -----------------------------------------------------------------------------
# PARAMETERS — set these before sourcing (guarded defaults)
#   H2               the heritability setting, default 0.6
#   M_C              the number of causal SNPs, default 5
#   MAF              minimum minor allele frequency of the SNPs, for the pruning and the
#                    GWAS, default 0.05
#   SEEDS            seeds of the simulation (causal SNPs, effects, plants, measurement
#                    errors), one population per seed, default 1:5
#   N_GWAS_SNPS      number of SNPs the GWAS tests: always the causal SNPs, plus a random
#                    draw from all other SNPs of the panel (not pruned, MAF >= MAF) up to
#                    this number (per seed); Inf or more than available: all, default Inf
#   TRUE_PARAMETERS  TRUE: the true parameters of the reaction norms are GWAS traits
#                    next to the PPIs (true_L, true_a, ...), default TRUE
# =============================================================================
options(warn = -1)  # silence warnings even if the project .Rprofile was not loaded

if (!exists("H2")) H2 <- 0.6
if (!exists("M_C")) M_C <- 5
if (!exists("MAF")) MAF <- 0.05
if (!exists("SEEDS")) SEEDS <- 1:5
if (!exists("TRUE_PARAMETERS")) TRUE_PARAMETERS <- TRUE
if (!exists("N_GWAS_SNPS")) N_GWAS_SNPS <- Inf
source(here::here("R", "pipeline", "load_pipeline.R"))
geno <- new.env()
load(here::here("data", "A-Datasets", "Genotyping_matrix_Negro_et_al_2019.Rdata"), envir = geno)
map <- readRDS(here::here("data", "A-Datasets", "SNPs_DL_windows_Negro_et_al_2019.rds"))
map <- map[map$SNP.name %in% colnames(geno$mat2), ]
map <- map[order(map$Chromosome, map$Position), ]            # genomic order
panel <- geno$mat2[, map$SNP.name]
rm(geno)
snps <- prune_snps(panel, r2 = 0.1, window = 50, chrom = map$Chromosome, MAF = MAF)
#snps <- snps[1:200,]
rm(map)   # panel: all SNPs of the same lines, from which the GWAS draws its SNPs (MAF filter in the GWAS)
for (seed in SEEDS) {
  scenario_name = paste("scenario", as.character(H2*100), as.character(M_C), as.character(MAF*100), paste0("seed", seed), sep = "_")
  POPULATION_DIR <- here::here("output", "populations", scenario_name)
  message("== ", scenario_name)

  #simulate_population(sim_config(h2 = 0.2, Mc = 5), POPULATION_DIR)   # true reaction norms -> data/, variance/
  simulate_population_from_snps(sim_config(h2 = H2, Mc = M_C, seed = seed), snps, POPULATION_DIR)   # -> data/, variance/
  RESOLUTION_DIR <- sample_resolution(POPULATION_DIR, 6)    # measured at 6 points -> n_6/data/
  compute_population_ppis(RESOLUTION_DIR)                   # PPIs -> n_6/mean_before/ppis/
  run_population_gwas(RESOLUTION_DIR, gwas_snps = panel, min_maf = MAF, true_parameters = TRUE_PARAMETERS,
                      n_gwas_snps = N_GWAS_SNPS)            # GWAS, power, heatmaps, MAF of the hits -> n_6/mean_before/gwas/, plots/
}
