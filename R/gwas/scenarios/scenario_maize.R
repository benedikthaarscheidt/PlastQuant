# =============================================================================
# scenario_maize.R — the empirical maize data through the pipeline
# =============================================================================
# WHAT IT DOES: Loads the maize data (244 varieties, traits biomass, leaf area and WUE in
#   8 environments: seasons S12, S13, S16, W13 under drought (WD) and well watered (WW))
#   and the SNPs of the varieties, writes them as a population folder (own data), and
#   runs the same steps as for simulated data: PPIs, GWAS, and the summary of the GWAS
#   hits (06_analyze_maize_gwas_results.R).
#   Environment values: the environments are not controlled, so every environment gets
#   an environment index as in Finlay-Wilkinson: for each trait its own index, the
#   environment mean of that trait as given in the FW/AMMI files of the data
#   (ENV_Biomass, ENV_Leaf_area, ENV_WUE). The environments of a trait are ordered
#   along its index.
# REQUIRES:     data/A-Datasets/: Biomass_RDPI_CV_Ratio.txt, Leaf_area_RDPI_CV_Ratio.txt,
#               WUE_RDPI_CV_Ratio.txt (trait values), Biomass_FW_AMMI.txt,
#               Leaf_area_FW_AMMI.txt, WUE_FW_AMMI.txt (environment indices),
#               Genotyping_matrix_Negro_et_al_2019.Rdata
#               (genotyping data: https://forgemia.inra.fr/gqe-base/djabali-drought-plasticity-qtls-article1)
# PRODUCES:     output/populations/maize/: data/ (SNPs), n_8/data/ (the three traits),
#               n_8/ppis/, n_8/gwas/, and the GWAS summary of 06_analyze_maize_gwas_results.R
# HOW TO RUN:   Rscript --vanilla R/gwas/scenarios/scenario_maize.R   (from the repo root)
# -----------------------------------------------------------------------------
# PARAMETERS — set these before sourcing (guarded defaults)
#   POPULATION_DIR  path, default output/populations/maize   population folder        [COMMON]
#   SNP_SUBSET      integer, default NULL (all 978134 SNPs)   use only the first SNPs (quick checks)
# =============================================================================
options(warn = -1)  # silence warnings even if the project .Rprofile was not loaded
source(here::here("R", "pipeline", "load_pipeline.R"))

if (!exists("POPULATION_DIR")) POPULATION_DIR <- here::here("output", "populations", "maize")
if (!exists("SNP_SUBSET")) SNP_SUBSET <- NULL
data_dir <- here::here("data", "A-Datasets")

# -----------------------------------------------------------------------------
# 1. Load the phenotypes: one file per trait, one column per environment
# -----------------------------------------------------------------------------
environment_names <- c("S12_WD", "S13_WD", "S16_WD", "W13_WD", "S12_WW", "S13_WW", "S16_WW", "W13_WW")
trait_files  <- c(biomass = "Biomass_RDPI_CV_Ratio.txt", leafArea = "Leaf_area_RDPI_CV_Ratio.txt",
                  WUE = "WUE_RDPI_CV_Ratio.txt")
trait_prefix <- c(biomass = "Biomass", leafArea = "Leaf_area", WUE = "WUE")
raw <- lapply(trait_files, function(f) read.delim(file.path(data_dir, f), check.names = FALSE, stringsAsFactors = FALSE))

# Environment index of every trait: its environment means from the FW/AMMI files
# (column ENV_<trait>, one value per environment ENV = <treatment>_<season>). The value
# an environment's rows share is taken, so a single mistyped row cannot change it.
index_files <- c(biomass = "Biomass_FW_AMMI.txt", leafArea = "Leaf_area_FW_AMMI.txt", WUE = "WUE_FW_AMMI.txt")
environments <- lapply(setNames(names(index_files), names(index_files)), function(trait) {
  d <- read.delim(file.path(data_dir, index_files[[trait]]), stringsAsFactors = FALSE)
  index_col <- paste0("ENV_", trait_prefix[[trait]])
  vapply(environment_names, function(e) {
    parts <- strsplit(e, "_")[[1]]                                   # S12_WD -> ENV WD_S12
    v <- d[[index_col]][d$ENV == paste(parts[2], parts[1], sep = "_")]
    as.numeric(names(which.max(table(v))))
  }, numeric(1))
})

# Trait matrices (Genotype = variety, one column per environment) -> long tables
phenotypes <- lapply(setNames(names(trait_files), names(trait_files)), function(trait) {
  d <- raw[[trait]]
  tab <- data.frame(Genotype = as.character(d$Variety), d[paste0(trait_prefix[[trait]], "_", environment_names)],
                    check.names = FALSE, stringsAsFactors = FALSE)
  names(tab)[-1] <- environment_names
  trait_matrix_to_long(tab, environments[[trait]])
})

# -----------------------------------------------------------------------------
# 2. Load the SNPs of the phenotyped varieties (dosages 0/1/2)
# -----------------------------------------------------------------------------
geno <- new.env()
load(file.path(data_dir, "Genotyping_matrix_Negro_et_al_2019.Rdata"), envir = geno)   # object mat2: lines x SNPs
snp_matrix <- geno$mat2[as.character(raw$biomass$Variety), , drop = FALSE]
rm(geno)
if (!is.null(SNP_SUBSET)) snp_matrix <- snp_matrix[, seq_len(min(SNP_SUBSET, ncol(snp_matrix))), drop = FALSE]
message("Maize: ", nrow(snp_matrix), " varieties, ", ncol(snp_matrix), " SNPs, ", length(environment_names), " environments")

# -----------------------------------------------------------------------------
# 3. Population folder: SNPs in data/, the measurements as the resolution n_8
# -----------------------------------------------------------------------------
write_population(pipeline_input(phenotypes, environments, snp_matrix), POPULATION_DIR, snp_format = "rds")
rm(snp_matrix)
RESOLUTION_DIR <- resolution_path(POPULATION_DIR, length(environment_names))

# -----------------------------------------------------------------------------
# 4. PPIs, GWAS and the summary of the GWAS hits
# -----------------------------------------------------------------------------
compute_population_ppis(RESOLUTION_DIR)      # -> n_8/ppis/
run_population_gwas(RESOLUTION_DIR)          # -> n_8/gwas/gwas_results.csv (no truth: no power)
source(here::here("R", "gwas", "06_analyze_maize_gwas_results.R"))
