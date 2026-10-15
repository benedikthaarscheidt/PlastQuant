# =============================================================================
# plasticity_score_realisation_plot_maize.R — plasticity-score realisations (maize)
# =============================================================================
# WHAT IT DOES: Plots the realised plasticity scores for the empirical maize dataset.
#   Sources analysis/rank_corr_biplot.R.
# REQUIRES:     The maize PPIs of scenario_maize.R (RESOLUTION_DIR, default
#               output/populations/maize/n_8); analysis/rank_corr_biplot.R.
# PRODUCES:     Maize plasticity-score realisation figure(s).
# HOW TO RUN:   setwd("~/PlastQuant"); source(here::here("R", "analysis", "plasticity_score_realisation_plot_maize.R"))
# -----------------------------------------------------------------------------
# PARAMETERS (edit the named assignment near the top of the script)
#   form_ranges  list  the reaction-norm ranges/forms plotted                              [COMMON]
#   RESOLUTION_DIR  path, default output/populations/maize/n_8   maize resolution folder
# =============================================================================
options(warn = -1)  # silence warnings even if the project .Rprofile was not loaded

if (!exists("RESOLUTION_DIR")) RESOLUTION_DIR <- here::here("output", "populations", "maize", "n_8")
form_ranges <- list(
  biomass     = 1:244,
  leafArea   = 245:488,
  WUE = 489:732
)

#source(here::here("R", "analysis", "plasticity_score_realisation_plot.R"))
source(here::here("R", "analysis", "rank_corr_biplot.R"))
