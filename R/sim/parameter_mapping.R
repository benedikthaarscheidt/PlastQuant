# =============================================================================
# parameter_mapping.R — latent value -> value of a reaction-norm parameter
# =============================================================================

# L = exp(sigma_L * ell); bounded parameters: Q_theta(Phi(ell)) on their study space.
# Works on vectors and matrices; all maps are increasing in ell.
latent_to_parameter <- function(ell, theta, shape, cfg) {
  if (theta == "L") return(exp(cfg$sigma_L * ell))
  space <- cfg$study_space[[shape]][[theta]]
  if (is.null(space)) stop("no study space for parameter ", theta, " of shape ", shape)
  u <- latent_to_unit(ell)
  if (space$type == "uniform") return(space$min + (space$max - space$min) * u)
  space$min * (space$max / space$min)^u
}
