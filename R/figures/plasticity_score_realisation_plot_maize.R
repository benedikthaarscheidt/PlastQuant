# =============================================================================
# plasticity_score_realisation_plot_maize.R — plasticity-score realisations (maize)
# =============================================================================
# WHAT IT DOES: Plots the realised plasticity scores for the empirical maize dataset.
#   Sources figures/rank_corr_biplot.R.
# REQUIRES:     Maize score outputs (data/maize_scores_list.rds and the maize pipeline
#               outputs); figures/rank_corr_biplot.R.
# PRODUCES:     Maize plasticity-score realisation figure(s).
# HOW TO RUN:   setwd("~/PlastQuant"); source(here::here("R","figures","plasticity_score_realisation_plot_maize.R"))
# -----------------------------------------------------------------------------
# PARAMETERS (edit the named assignment near the top of the script)
#   form_ranges  list  the reaction-norm ranges/forms plotted                              [COMMON]
#   the input score path is set in the script body.
# =============================================================================
options(warn = -1)  # silence warnings even if the project .Rprofile was not loaded

form_ranges <- list(
  biomass     = 1:244,
  leafArea   = 245:488,
  WUE = 489:732
)

#source(here::here("R", "figures", "plasticity_score_realisation_plot.R"))
source(here::here("R", "figures", "rank_corr_biplot.R"))
