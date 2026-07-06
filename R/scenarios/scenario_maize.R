# =============================================================================
# scenario_maize.R — scenario driver: empirical maize dataset
# =============================================================================
# WHAT IT DOES: Sets the parameters for the empirical maize analysis and runs the GWAS
#   accuracy evaluation (05_evaluate_gwas.R). The maize score computation
#   (03_plasticity_scores_maize.R) is available but sourced separately.
# REQUIRES:     None to launch — self-contained driver. Uses data/maize_scores_list.rds
#               and the maize datasets under data/A-Datasets/.
# PRODUCES:     Maize GWAS evaluation under OUTPUT_BASE (output/scenario_maize/).
# HOW TO RUN:   setwd("~/PlastQuant"); source(here::here("R","scenarios","scenario_maize.R"))
# -----------------------------------------------------------------------------
# PARAMETERS (edit here — set UNCONDITIONALLY below so they override script defaults)
#   OUTPUT_BASE           output location for the maize run
#   DO_NOT_PLOT_MAHATTAN  = TRUE  skip Manhattan plots
# =============================================================================
options(warn = -1)  # silence warnings even if the project .Rprofile was not loaded

OUTPUT_BASE = here::here("output", "scenario_maize")
DO_NOT_PLOT_MAHATTAN = TRUE

#source(here::here("R", "03_plasticity_scores_maize.R"))

source(here::here("R", "05_evaluate_gwas.R"))


