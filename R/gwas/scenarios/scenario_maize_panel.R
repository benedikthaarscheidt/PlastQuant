# =============================================================================
# scenario_maize_panel.R — simulated reaction norms on the real maize SNPs
# =============================================================================
# WHAT IT DOES: Takes the SNPs of the maize panel (Negro et al. 2019) as genotypes and
#   simulates the reaction norms on them (simulate_population_from_snps(), snp_panel.R),
#   so that LD, allele frequencies and relatedness are those of real maize lines. Then
#   measures the population at N_POINTS points, computes the PPIs and runs the GWAS with
#   the power per shape, as for a simulated population.
# REQUIRES:     data/A-Datasets/Genotyping_matrix_Negro_et_al_2019.Rdata (object mat2:
#               lines x SNPs, dosages 0/1/2)
# PRODUCES:     POPULATION_DIR/data/ (incl. causal_tags.csv), variance/, n_<N_POINTS>/data/,
#               n_<N_POINTS>/<REPLICATE_HANDLING>/ppis/, gwas/, plots/
# HOW TO RUN:   from the repo root, from the terminal (the GWAS runs on N_CORES cores):
#               OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 Rscript --vanilla R/gwas/scenarios/scenario_maize_panel.R
# -----------------------------------------------------------------------------
# PARAMETERS — set these before sourcing (guarded defaults)
#   POPULATION_DIR      path, default output/populations/maize_panel_<h2*100>_<Mc>
#   H2                  numeric (0.8)          narrow-sense heritability on the latent scale
#   MC                  integer (5)            causal SNPs per parameter
#   SEED                integer (1)            seed of the genetics
#   SNP_SUBSET          integer, NULL (all)    use only the first SNPs (quick checks)
#   N_POINTS            integer (6)            resolution
#   REPLICATE_HANDLING  "mean_before"          replicate handling of the PPIs
#   N_CORES             integer (1)            cores for the GWAS
# All other settings are the defaults of sim_config() in R/sim/simulation_config.R.
# =============================================================================
options(warn = -1)  # silence warnings even if the project .Rprofile was not loaded
source(here::here("R", "pipeline", "load_pipeline.R"))

if (!exists("H2")) H2 <- 0.8
if (!exists("MC")) MC <- 5
if (!exists("SEED")) SEED <- 1
if (!exists("SNP_SUBSET")) SNP_SUBSET <- NULL
if (!exists("N_POINTS")) N_POINTS <- 6
if (!exists("REPLICATE_HANDLING")) REPLICATE_HANDLING <- "mean_before"
if (!exists("N_CORES")) N_CORES <- 1
if (!exists("POPULATION_DIR")) {
  POPULATION_DIR <- here::here("output", "populations", paste0("maize_panel_", round(H2 * 100), "_", MC))
}

# SNPs of all lines of the panel
geno <- new.env()
load(here::here("data", "A-Datasets", "Genotyping_matrix_Negro_et_al_2019.Rdata"), envir = geno)
snps <- geno$mat2
rm(geno)
if (!is.null(SNP_SUBSET)) snps <- snps[, seq_len(min(SNP_SUBSET, ncol(snps))), drop = FALSE]
message("Maize panel: ", nrow(snps), " lines, ", ncol(snps), " SNPs")

simulate_population_from_snps(sim_config(h2 = H2, Mc = MC, seed = SEED), snps, POPULATION_DIR)   # -> data/, variance/
rm(snps)
RESOLUTION_DIR <- sample_resolution(POPULATION_DIR, N_POINTS)                  # -> n_<N_POINTS>/data/
compute_population_ppis(RESOLUTION_DIR, REPLICATE_HANDLING)                     # -> <handling>/ppis/
run_population_gwas(RESOLUTION_DIR, REPLICATE_HANDLING, n_cores = N_CORES)      # -> <handling>/gwas/, plots/
