# =============================================================================
# scenario_40_20.R — scenario driver: HERITABILITY_5TH=2, CAUSAL_SNP_NUM=20
# =============================================================================
# WHAT IT DOES: Sets the parameters for this scenario, then runs the plasticity-score
#   computation (03_plasticity_scores.R) and the GWAS accuracy evaluation
#   (05_evaluate_gwas.R) end-to-end.
# REQUIRES:     None — self-contained driver. Run from the project root.
# PRODUCES:     Scores + GWAS evaluation under OUTPUT_BASE (output/scenario_40_20/).
# HOW TO RUN:   setwd("~/PlastQuant"); source(here::here("R","scenarios","scenario_40_20.R"))
# -----------------------------------------------------------------------------
# PARAMETERS (edit here — set UNCONDITIONALLY below so they override script defaults)
#   OUTPUT_BASE           output location for this scenario
#   HERITABILITY_5TH      = 2    heritability setting for the 5th trait          [COMMON]
#   CAUSAL_SNP_NUM        = 20    causal SNPs per parameter                       [COMMON]
#   DO_NOT_PLOT_MAHATTAN  = TRUE  skip Manhattan plots
# =============================================================================
options(warn = -1)  # silence warnings even if the project .Rprofile was not loaded

OUTPUT_BASE = here::here("output", "scenario_40_20")
HERITABILITY_5TH = 2
CAUSAL_SNP_NUM = 20
DO_NOT_PLOT_MAHATTAN = TRUE

source(here::here("R", "03_plasticity_scores.R"))

source(here::here("R", "05_evaluate_gwas.R"))
