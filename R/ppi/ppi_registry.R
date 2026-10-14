# =============================================================================
# ppi_registry.R — which plasticity indices (PPIs) are computed, and with what
# =============================================================================
# Every PPI is one entry. `inputs` maps the arguments of the function to what the
# pipeline provides:
#   level "genotype"   (called once per genotype; its replicate-mean reaction norm)
#       trait      numeric vector, trait value per environment
#       env        numeric vector, environment values
#       env_index  1, 2, ..., number of environments
#       covariate  vector of 1s (one per environment)
#   level "population" (called once per shape)
#       matrix     genotypes x environments matrix of replicate means
#       env        numeric vector, environment values
# `args` holds fixed extra arguments. `extract` turns the result into one number per
# genotype (default: the first element). To add a PPI, append an entry with
# register_ppi(); to leave PPIs out, use the `exclude` argument of compute_ppis().
# =============================================================================

ppi_entry <- function(name, fun, inputs = c(trait_values = "trait"), args = list(),
                      level = "genotype", extract = NULL) {
  if (!level %in% c("genotype", "population")) stop("level must be 'genotype' or 'population'")
  list(name = name, fun = fun, inputs = inputs, args = args, level = level, extract = extract)
}

register_ppi <- function(registry, entry) {
  registry[[entry$name]] <- entry
  validate_ppi_registry(registry[entry$name])
  registry
}

validate_ppi_registry <- function(registry) {
  keys <- list(genotype = c("trait", "env", "env_index", "covariate"), population = c("matrix", "env"))
  for (e in registry) {
    bad_key <- setdiff(e$inputs, keys[[e$level]])
    if (length(bad_key)) stop("PPI ", e$name, ": unknown input ", paste(bad_key, collapse = ", "))
    formal_names <- names(formals(e$fun))
    if (!"..." %in% formal_names) {
      bad_arg <- setdiff(c(names(e$inputs), names(e$args)), formal_names)
      if (length(bad_arg)) stop("PPI ", e$name, ": function has no argument ", paste(bad_arg, collapse = ", "))
    }
  }
  invisible(TRUE)
}

# The PPIs of the study, with the inputs they were computed with so far
default_ppi_registry <- function() {
  pp <- asNamespace("ppindices")
  trait_only <- c(
    CV_t = "calculate_CVt", RN = "calculate_reaction_norm_slope", D_slope = "calculate_D_slope",
    RC = "calculate_RC", RDPI = "calculate_rdpi", ESPI = "calculate_ESPI",
    ESPIID = "calculate_espiid", PSI = "calculate_PSI", RPI = "calculate_RPI", PQ = "calculate_PQ",
    PR = "calculate_PR", NRW = "calculate_NRW", ESP = "calculate_ESP", CEV = "calculate_CEV",
    PFI = "calculate_PFI", APC = "calculate_APC", SI = "calculate_SI", RSI = "calculate_RSI",
    EVS = "calculate_EVS"
  )
  reg <- lapply(names(trait_only), function(n) ppi_entry(n, get(trait_only[[n]], envir = pp)))
  names(reg) <- names(trait_only)
  reg$PPi   <- ppi_entry("PPi", pp$calculate_Phenotypic_Plasticity_Index, inputs = c(data = "trait"))
  reg$RNN   <- ppi_entry("RNN", pp$calculate_reaction_norm_non_linear, args = list(degree = 3))
  reg$gPi   <- ppi_entry("gPi", pp$calculate_grand_plasticity,
                         inputs = c(trait_values = "trait", env_data = "env", covariate_data = "covariate"))
  reg$PPF   <- ppi_entry("PPF", pp$calculate_PPF,
                         inputs = c(trait_values = "trait", env_groups = "env", covariate_values = "covariate"))
  reg$PImd  <- ppi_entry("PImd", pp$calculate_PImd, inputs = c(trait_values = "trait", env = "env"))
  reg$PILSM <- ppi_entry("PILSM", pp$calculate_PILSM,
                         inputs = c(trait_values = "trait", env = "env", covariates = "covariate"))
  reg$RTR   <- ppi_entry("RTR", pp$calculate_RTR, inputs = c(trait_values = "trait", env_values = "env"))
  reg$PIR   <- ppi_entry("PIR", pp$calculate_PIR, inputs = c(trait_values = "trait", env_values = "env"))
  reg$FW    <- ppi_entry("FW", pp$calculate_finlay_wilkinson, inputs = c(Y = "matrix"),
                         level = "population", extract = function(res) setNames(res$beta, res$genotype))
  validate_ppi_registry(reg)
  reg
}
