#this files contains the following indices/functions:
#generate_synthetic_data,
#Phenotypic Stability Index (calculate_PSI),
#Relative Plasticity Index (calculate_RPI) - tested,
#Plasticity Quotient (calculate_PQ) - tested,
#Phenotypic Range (PR) (calculate_PR) - tested,
#Norm of reaction width (calculate_NRW) - tested,
#Environment-Specific Plasticity (ESP) (calculate_ESP) - tested,
#Calculate Plasticity Differential (PD) (calculate_PD) - tested,
#Calculate Fitness Plasticity Index (FPI) (calculate_FPI) - tested,
#Calculate Transplant Plasticity Score (TPS)(calculate_TPS) - tested,



##################### datasets for testing

#df_test1 = data.frame(Column1 = c(rep(4, 10), rep(2, 10)), Column2 = c(rep(10, 10), rep(1, 10)))
#df_test2 = data.frame(Column0 = c(rep(1, 10), rep(2, 10)),Column1 = c(rep(2, 10), rep(1, 10)), Column2 = c(rep(2, 10), rep(4, 10)))
#df_test3=data.frame(Column0 = c(rep(3, 10), rep(2, 10),rep(1,10)),Column1 = c(rep(2, 10), rep(1, 15),rep(3,5)), Column2 = c(rep(2, 10), rep(4, 10),rep(3,10)))
#df_test4 = data.frame(
#  Column0 = c(rep(3, 10), rep(2, 10), rep(1, 10)),   # Response variable
#  Column1 = c(rep(2, 10), rep(3, 10), rep(1, 10)),   # Control (2) and Treatment (3)
#  Column2 = c(rep(1, 10), rep(1, 10), rep(1, 10))    # Covariate (matches values of Column0)
#)
#
#df_test5 = data.frame(
#  Column0 = c(rep(3, 10), rep(2, 10), rep(1, 10)),   # Response variable
#  Column1 = c(rep(3, 10), rep(2, 10), rep(1, 10)),   # Control (2) and Treatment (3)
#  Column2 = c(rep(3, 10), rep(2, 10), rep(1, 10))    # Covariate (matches values of Column0)
#)
#
#df_test_simple = data.frame(
#  Column1 = c(2, 2, 3, 3),      # Control (2) and Treatment (3)
#  Column0 = c(10, 12, 20, 22),  # Response variable (trait)
#  Column2 = c(1, 1, 2, 2)       # Covariate (e.g., biomass)
#)
#
#df_test6=data.frame(
#  Column0 = c(rep(3, 10), rep(2, 10), rep(1, 10)),
#  Column1 = c(rep(6, 10), rep(4, 10), rep(2, 10)),
#  Column2 = c(rep(1, 5),rep(2,5), rep(2, 10), rep(3, 10)),
#  Column3 = c(rep(1, 10), rep(1, 10), rep(1, 10))
#)



########################################

#' @title Calculate the Phenotypic Stability Index (PSI) Using Linear Regression
#'
#' @description
#' Computes the Phenotypic Stability Index (PSI) for a trait measured across environments.
#' PSI is defined as \(1 / (1 + |β|)\), where \(β\) is the slope from regressing trait values
#' on environmental values. A smaller absolute slope (i.e., more stable trait) yields a PSI closer to 1.
#'
#' @param trait_values
#'   Numeric vector of trait measurements. Length must be ≥ 2.
#' @param env_values
#'   (Optional) Numeric vector of environmental values corresponding to `trait_values`.
#'   If `NULL` (default), environments are assumed equidistant: `1, 2, …, length(trait_values)`.
#'
#' @return
#'   A single numeric PSI value:
#'   \deqn{\mathrm{PSI} = \frac{1}{1 + |\beta|},}
#'   where \(\beta\) is the regression coefficient from \(\mathrm{trait\_values} \sim \mathrm{env\_values}\).
#'
#' @details
#' If `env_values` has near-zero variance, a warning is issued since the regression slope
#' (and hence PSI) may be unstable.
#'
#' @examples
#' # Example trait in three environments, two replicates each
#' trait_values <- c(100, 110, 120, 105, 115, 125)
#' env_values   <- c(1, 2, 3,   1,   2,   3)
#'
#' # Calculate PSI with explicit env_values
#' PSI_result <- calculate_PSI(trait_values, env_values)
#' print(PSI_result)
#'
#' # Calculate PSI assuming equidistant environments
#' PSI_no_env <- calculate_PSI(trait_values)
#' print(PSI_no_env)
#'
#' @export
calculate_PSI <- function(trait_values, env_values = NULL) {
  # Input validation
  if (!is.numeric(trait_values)) {
    stop("`trait_values` must be a numeric vector.")
  }

  n_values <- length(trait_values)
  if (n_values < 2) {
    stop("At least two trait values are required to compute PSI.")
  }

  # Default equidistant environments
  if (is.null(env_values)) {
    env_values <- seq_len(n_values)
  }

  if (!is.numeric(env_values)) {
    stop("`env_values` must be numeric.")
  }
  if (length(env_values) != n_values) {
    stop("`trait_values` and `env_values` must have the same length.")
  }

  # Warn if env_values have near-zero variance
  if (sd(env_values) < .Machine$double.eps) {
    warning("Environmental values have near-zero variance; PSI may be unreliable.")
  }

  # Fit linear model and extract slope
  beta <- unname(coef(lm(trait_values ~ env_values))[2])

  # Compute stability index
  stability_score <- 1 / (1 + abs(beta))
  return(stability_score)
}





#external_light = rep(c(0.4, 0.6, 0.8), each = 100)
#PSI = calculate_PSI(synthetic_data1, trait_cols = 5, external_env_factor = external_light)
#print(PSI)

##################################


#' @title Calculate the Relative Plasticity Index (RPI)
#'
#' @description
#' Computes the Relative Plasticity Index (RPI) of a trait across environments.
#' For two environments, RPI is defined as the mean of \eqn{|x_{i,1} - x_{i,2}| / (x_{i,1} + x_{i,2})}
#' over paired observations. When no specific pair is given, RPI is averaged over all unique
#' environment pairs.
#'
#' @param trait_values
#'   A numeric vector of trait measurements.
#' @param env_values
#'   (Optional) A vector of the same length as \code{trait_values} giving the environment
#'   label or value for each observation. If \code{NULL} (default), environments are assumed
#'   equidistant (1, 2, …, length(trait_values)).
#' @param env1
#'   (Optional) A single environment label or value. Must appear in \code{env_values}.
#'   If provided together with \code{env2}, RPI is computed only for that pair.
#' @param env2
#'   (Optional) A single environment label or value. Must appear in \code{env_values}.
#'   If provided together with \code{env1}, RPI is computed only for that pair.
#'
#' @return
#' A single numeric value:
#' \itemize{
#'   \item If \code{env1} and \code{env2} are provided, returns the mean
#'     of \eqn{|x_{i,env1} - x_{i,env2}| / (x_{i,env1} + x_{i,env2})} over paired observations.
#'   \item Otherwise, returns the grand mean of RPI across all unique environment pairs.
#' }
#' Returns \code{NA} (with a warning) if fewer than two unique environments are present.
#'
#' @details
#' 1. Validates that \code{trait_values} is numeric and that \code{env_values} length matches.
#' 2. If fewer than two unique environments exist, returns \code{NA}.
#' 3. If \code{env1} and \code{env2} are specified, pairs up observations (up to the smaller group size)
#'    and computes RPI for that pair.
#' 4. Otherwise, iterates over every combination of two distinct environments, computes pairwise RPI,
#'    and returns the average.
#'
#' @examples
#' # Two-environment example
#' trait <- c(5, 6, 7, 8, 9, 10)
#' env   <- rep(c("A", "B"), each = 3)
#' calculate_RPI(trait, env, env1 = "A", env2 = "B")
#'
#' # Multi-environment example: averaged over all pairs
#' trait <- c(2, 4, 6, 8, 10, 12)
#' env   <- rep(1:3, each = 2)
#' calculate_RPI(trait, env)
#'
#' # Equidistant environments assumed if env_values = NULL
#' calculate_RPI(c(1, 2, 3, 4))
#'
#' @export
calculate_RPI = function(trait_values, env_values = NULL, env1 = NULL, env2 = NULL) {
  # Ensure trait_values is numeric
  if (!is.numeric(trait_values)) stop("trait_values must be a numeric vector.")

  n = length(trait_values)

  # Assume equidistant environments if not provided
  if (is.null(env_values)) {
    env_values = seq_len(n)
  }

  # Ensure env_values has correct length
  if (length(env_values) != n) stop("env_values must have the same length as trait_values.")

  unique_envs = unique(env_values)

  # If only one environment is present, return NA
  if (length(unique_envs) < 2) {
    warning("Only one unique environment present. RPI cannot be calculated.")
    return(NA)
  }

  # If specific environments are given, calculate RPI for that pair
  if (!is.null(env1) && !is.null(env2)) {
    if (!(env1 %in% unique_envs) || !(env2 %in% unique_envs)) stop("env1 and env2 must be in env_values.")

    idx1 = which(env_values == env1)
    idx2 = which(env_values == env2)

    # Ensure equal sample size
    min_len = min(length(idx1), length(idx2))
    idx1 = idx1[seq_len(min_len)]
    idx2 = idx2[seq_len(min_len)]

    RPI_values = abs(trait_values[idx1] - trait_values[idx2]) / (trait_values[idx1] + trait_values[idx2])
    return(mean(RPI_values, na.rm = TRUE))
  }

  # Otherwise, compute RPI for all pairs
  env_combinations = combn(unique_envs, 2, simplify = TRUE)
  RPI_results = numeric(ncol(env_combinations))

  for (i in seq_len(ncol(env_combinations))) {
    env1_current = env_combinations[1, i]
    env2_current = env_combinations[2, i]

    idx1 = which(env_values == env1_current)
    idx2 = which(env_values == env2_current)

    min_len = min(length(idx1), length(idx2))
    idx1 = idx1[seq_len(min_len)]
    idx2 = idx2[seq_len(min_len)]

    RPI_values = abs(trait_values[idx1] - trait_values[idx2]) / (trait_values[idx1] + trait_values[idx2])
    RPI_results[i] = mean(RPI_values, na.rm = TRUE)
  }
  RPI_grand = sum(RPI_results) / length(RPI_results)
  names(RPI_results) = apply(env_combinations, 2, paste, collapse = "-")
  return(RPI_grand)
}
calculate_RPI = function(trait_values, env_values = NULL, env1 = NULL, env2 = NULL) {
  # Ensure trait_values is numeric
  if (!is.numeric(trait_values)) stop("trait_values must be a numeric vector.")

  n = length(trait_values)

  # Assume equidistant environments if not provided
  if (is.null(env_values)) {
    env_values = seq_len(n)
  }

  # Ensure env_values has correct length
  if (length(env_values) != n) stop("env_values must have the same length as trait_values.")

  unique_envs = unique(env_values)

  # If only one environment is present, return NA
  if (length(unique_envs) < 2) {
    warning("Only one unique environment present. RPI cannot be calculated.")
    return(NA)
  }

  # If specific environments are given, calculate RPI for that pair
  if (!is.null(env1) && !is.null(env2)) {
    if (!(env1 %in% unique_envs) || !(env2 %in% unique_envs)) stop("env1 and env2 must be in env_values.")

    idx1 = which(env_values == env1)
    idx2 = which(env_values == env2)

    # Ensure equal sample size
    min_len = min(length(idx1), length(idx2))
    idx1 = idx1[seq_len(min_len)]
    idx2 = idx2[seq_len(min_len)]

    RPI_values = abs(trait_values[idx1] - trait_values[idx2]) / (trait_values[idx1] + trait_values[idx2])
    return(mean(RPI_values, na.rm = TRUE))
  }

  # Otherwise, compute RPI for all pairs
  env_combinations = combn(unique_envs, 2, simplify = TRUE)
  RPI_results = numeric(ncol(env_combinations))

  for (i in seq_len(ncol(env_combinations))) {
    env1_current = env_combinations[1, i]
    env2_current = env_combinations[2, i]

    idx1 = which(env_values == env1_current)
    idx2 = which(env_values == env2_current)

    min_len = min(length(idx1), length(idx2))
    idx1 = idx1[seq_len(min_len)]
    idx2 = idx2[seq_len(min_len)]

    RPI_values = abs(trait_values[idx1] - trait_values[idx2]) / (trait_values[idx1] + trait_values[idx2])
    RPI_results[i] = mean(RPI_values, na.rm = TRUE)
  }
  RPI_grand=sum(RPI_results)/length(RPI_results)
  names(RPI_results) = apply(env_combinations, 2, paste, collapse = "-")
  return(RPI_grand)
}



## test - passed on synthetic dataset

#calculate_RPI(df_test2,trait_col=c(2,3),env_col=1)

######################################

#' @title Calculate Plasticity Quotient (PQ)
#' @description Computes the Plasticity Quotient (PQ) based on the range of trait values
#' and the corresponding environmental factor range.
#'
#' @param trait_values A numeric vector representing the measured trait values across environments.
#' @param env_values (Optional) A numeric vector representing the environmental conditions
#' in which the traits were measured. If NULL, equidistant environments are assumed.
#'
#' @details
#' The Plasticity Quotient (PQ) is calculated as:
#'
#' \deqn{PQ = \frac{\max(trait\_values) - \min(trait\_values)}{\max(env\_values) - \min(env\_values)}}
#'
#' If `env_values` is not provided, it is assumed that the environments are equidistant (i.e.,
#' `env_values = 1, 2, 3, ..., n` where `n` is the number of trait values).
#'
#' If the environmental range is zero (i.e., all values are identical), a warning is issued and `NA` is returned.
#'
#' @return A numeric value representing the Plasticity Quotient (PQ).
#'
#' @examples
#' # Example 1: Using explicitly defined environmental values
#' trait_values = c(10, 12, 15, 17, 20)
#' env_values = c(5, 10, 15, 20, 25)
#' calculate_PQ(trait_values, env_values)
#'
#' # Example 2: Using equidistant environments (default behavior)
#' trait_values = c(10, 15, 20, 25, 30)
#' calculate_PQ(trait_values)
#'
#' # Example 3: Handling zero environment range (returns NA with warning)
#' trait_values = c(10, 15, 20, 25, 30)
#' env_values = rep(5, 5)
#' calculate_PQ(trait_values, env_values)
#'
#' @export
calculate_PQ = function(trait_values, env_values = NULL) {
  # Ensure trait_values is numeric
  if (!is.numeric(trait_values)) stop("trait_values must be a numeric vector.")

  n = length(trait_values)

  # If no explicit env_values are given, assume equidistant environments
  if (is.null(env_values)) {
    env_values = seq_len(n)
  }

  # Ensure env_values is numeric and of correct length
  if (!is.numeric(env_values)) stop("env_values must be a numeric vector.")
  if (length(env_values) != n) stop("trait_values and env_values must have the same length.")

  # Compute the range of trait values
  range_trait = abs(max(trait_values) - min(trait_values))

  # Compute the range of environment values (used for standardization)
  range_env = abs(max(env_values) - min(env_values))

  # Ensure env range is non-zero to avoid division errors
  if (range_env == 0) {
    warning("Environment range is zero. PQ calculation is not meaningful.")
    return(NA)
  }

  # Compute the Plasticity Quotient (PQ)
  PQ_value = range_trait / range_env

  return(PQ_value)
}


## test - passed on synthetic dataframe

#calculate_PQ(df_test6,1,3,2)

###########################################

#' @title Calculate Phenotypic Range (PR)
#' @description Computes the Phenotypic Range (PR) for a trait across different environments.
#'
#' @param trait_values A numeric vector representing the measured trait values.
#' @param env_values (Optional) A numeric vector representing the environmental conditions.
#' If NULL, equidistant environments are assumed.
#' @param across A logical value. If TRUE, calculates the PR across all environments combined.
#'
#' @details
#' The Phenotypic Range (PR) is defined as:
#'
#' \deqn{PR = \max(trait\_values) - \min(trait\_values)}
#'
#' If `across = TRUE`, PR is calculated across **all** environments.
#'
#' If `env_values` is **not provided**, equidistant environments are assumed.
#'
#' @return A numeric value representing the Phenotypic Range (PR).
#'
#' @examples
#' # Example 1: Using explicitly defined environmental values
#' trait_values = c(10, 12, 15, 17, 20)
#' env_values = c(5, 10, 15, 20, 25)
#' calculate_PR(trait_values, env_values)
#'
#' # Example 2: Using equidistant environments (default behavior)
#' trait_values = c(10, 15, 20, 25, 30)
#' calculate_PR(trait_values)
#'
#' # Example 3: Calculating PR across all environments
#' calculate_PR(trait_values, env_values, across = TRUE)
#'
#' @export
calculate_PR = function(trait_values, env_values = NULL, across = TRUE) {
  # Ensure trait_values is numeric
  if (!is.numeric(trait_values)) stop("trait_values must be a numeric vector.")

  n = length(trait_values)

  # If no explicit env_values are given, assume equidistant environments
  if (is.null(env_values)) {
    env_values = seq_len(n)
  }

  # Ensure env_values is numeric and matches trait_values in length
  if (!is.numeric(env_values)) stop("env_values must be a numeric vector.")
  if (length(env_values) != n) stop("trait_values and env_values must have the same length.")

  # If across = TRUE, compute PR across all environments
  if (across) {
    PR_value = max(trait_values, na.rm = TRUE) - min(trait_values, na.rm = TRUE)
    return(PR_value)
  }

  # Compute PR within each environment
  unique_envs = unique(env_values)
  PR_values = numeric(length(unique_envs))

  for (i in seq_along(unique_envs)) {
    env_mask = env_values == unique_envs[i]
    env_trait_values = trait_values[env_mask]
    PR_values[i] = max(env_trait_values, na.rm = TRUE) - min(env_trait_values, na.rm = TRUE)
  }

  return(PR_values)
}



## test - passed on synthetic dataset


#calculate_PR(df_test6,1,3)
#calculate_PR(df_test6,1,trait_cols=c(2,3),across=T)






##########################

#' @title Calculate Norm of Reaction Width (NRW)
#' @description Computes the Norm of Reaction Width (NRW) for a given trait across multiple environments.
#'
#' @param trait_values A numeric vector representing measured trait values.
#' @param env_values (Optional) A numeric vector representing environments. If NULL, assumes equidistant environments.
#' @param group_values (Optional) A vector indicating genotypes or other grouping factors. If provided, NRW is computed per group.
#' @param across Logical. If TRUE, calculates NRW **across all environments**.
#'
#' @details
#' NRW is calculated as:
#' \deqn{NRW = \max(\bar{z}_e) - \min(\bar{z}_e)}
#' where **\(\bar{z}_e\)** is the mean trait value per environment.
#'
#' If `across = TRUE`, NRW is computed **for all environments**.
#'
#' If an environment has **only one measurement**, NRW is set to `NA`.
#'
#' @return A numeric value or named vector representing NRW.
#'
#' @examples
#' trait_values = c(10, 12, 15, 17, 20, 25, 30)
#' env_values = c(1, 1, 2, 2, 3, 3, 3)
#' calculate_NRW(trait_values, env_values)
#'
#' @export
calculate_NRW = function(trait_values, env_values = NULL, group_values = NULL, across = FALSE) {
  if (!is.numeric(trait_values)) stop("trait_values must be numeric.")

  # Assume equidistant environments if not provided
  if (is.null(env_values)) env_values = seq_along(trait_values)
  if (length(env_values) != length(trait_values)) stop("trait_values and env_values must be the same length.")

  # Compute NRW across all environments
  if (across) {
    return(max(trait_values, na.rm = TRUE) - min(trait_values, na.rm = TRUE))
  }

  # Compute NRW per group (e.g., genotype)
  if (!is.null(group_values)) {
    unique_groups = unique(group_values)
    NRW_values = setNames(numeric(length(unique_groups)), unique_groups)

    for (g in unique_groups) {
      mask = group_values == g
      mean_per_env = tapply(trait_values[mask], env_values[mask], mean, na.rm = TRUE)
      NRW_values[g] = max(mean_per_env, na.rm = TRUE) - min(mean_per_env, na.rm = TRUE)
    }

    return(NRW_values)
  }

  # Default: Calculate NRW across environments
  mean_per_env = tapply(trait_values, env_values, mean, na.rm = TRUE)
  NRW_value = max(mean_per_env, na.rm = TRUE) - min(mean_per_env, na.rm = TRUE)

  return(NRW_value)
}


## test - passed on synthetic dataset


#calculate_NRW(df_test6,trait_cols=3)

##############################


#' @title Calculate Environmental Sensitivity Performance (ESP)
#' @description Computes the Environmental Sensitivity Performance (ESP) for multiple traits across environments.
#'
#' @param trait_values A numeric vector representing measured trait values.
#' @param env_values (Optional) A numeric vector representing environments. If NULL, assumes equidistant environments.
#' @param env_subset (Optional) A vector specifying which environments to include in the calculation. If NULL, all unique environments are used.
#'
#' @details
#' ESP is calculated as:
#' \deqn{ESP = \frac{\bar{z}_e - \bar{z}}{\bar{z}}}
#' where **\(\bar{z}_e\)** is the mean trait value in environment \( e \) and **\(\bar{z}\)** is the overall mean across all environments.
#'
#' If `env_values` is NULL, environments are assumed to be **equidistant**.
#'
#' @return A named numeric vector containing ESP values for each environment.
#'
#' @examples
#' trait_values = c(10, 12, 15, 17, 20, 25, 30)
#' env_values = c(1, 1, 2, 2, 3, 3, 3)
#' calculate_ESP(trait_values, env_values)
#'
#' @export
calculate_ESP = function(trait_values, env_values = NULL, env_subset = NULL) {
  if (!is.numeric(trait_values)) stop("trait_values must be numeric.")

  # Assume equidistant environments if not provided
  if (is.null(env_values)) env_values = seq_along(trait_values)
  if (length(env_values) != length(trait_values)) stop("trait_values and env_values must be the same length.")

  # Use all unique environments unless specific ones are given
  unique_envs = unique(env_values)
  if (!is.null(env_subset)) unique_envs = unique_envs[unique_envs %in% env_subset]

  # Compute mean trait value across all environments
  mean_trait_all = mean(trait_values, na.rm = TRUE)

  # Compute ESP per environment
  ESP_values = setNames(numeric(length(unique_envs)), unique_envs)

  for (e in unique_envs) {
    mean_trait_env = mean(trait_values[env_values == e], na.rm = TRUE)
    ESP_values[as.character(e)] = (mean_trait_env - mean_trait_all) / mean_trait_all
  }

  y=abs(ESP_values)
  ESP_grand=sum(y)
  return(ESP_grand)
}



## test - passed on synthetic dataset

#calculate_ESP(df_test6,1,trait_cols=c(2,3),env=NULL)

######################### 01.10.2024

#' @title Calculate Generalized Plasticity Differential (PD)
#' @description Computes the Plasticity Differential (PD) for traits across environments.
#' If control vs. stress environments are explicitly given, it calculates PD based on those.
#' Otherwise, PD is computed using pairwise environment comparisons.
#'
#' @param trait_values A numeric vector of measured trait values.
#' @param env_values Optional. A numeric or categorical vector representing environments.
#'   If NULL, assumes equidistant environments.
#' @param control_stress_vector Optional. A vector indicating control vs. stress conditions (e.g., "Control"/"Stress" or `0`/`1`).
#'   If provided, PD is calculated as |Trait_Stress - Trait_Control|.
#' @param method A string specifying the calculation approach:
#'   - `"pairwise"`: Computes mean absolute difference across all environment pairs (default if `control_stress_vector` is missing).
#'   - `"reference"`: Uses the lowest mean trait environment as reference.
#'   - `"variability"`: Uses the range (max - min) of trait values as PD.
#'
#' @return A single numeric value representing the computed PD.
#'
#' @examples
#' trait_values = c(100, 110, 120, 130, 140, 150)
#' control_stress = c(0, 0, 1, 1, 1, 1)  # Numeric encoding of Control (0) and Stress (1)
#'
#' # Explicit control vs. stress comparison
#' calculate_general_PD(trait_values, control_stress_vector = control_stress)
#'
#' # Pairwise comparisons with equidistant environments (default)
#' calculate_general_PD(trait_values)
#'
#' # Variability-based PD
#' calculate_general_PD(trait_values, method = "variability")
#'
#' @export
calculate_general_PD = function(trait_values, env_values = NULL, control_stress_vector = NULL, method = "pairwise") {
  if (!is.numeric(trait_values)) stop("trait_values must be numeric.")

  num_values = length(trait_values)
  if (num_values < 2) stop("At least two trait values are required.")

  # If no env_values are provided, assume equidistant environments
  if (is.null(env_values)) {
    env_values = seq_len(num_values)
  }

  if (length(env_values) != num_values) stop("trait_values and env_values must have the same length.")

  unique_envs = unique(env_values)
  num_envs = length(unique_envs)

  if (num_envs < 2) stop("At least two distinct environments are required.")

  # If a control vs. stress vector is provided, use it directly
  if (!is.null(control_stress_vector)) {
    if (length(control_stress_vector) != num_values) {
      stop("control_stress_vector must have the same length as trait_values.")
    }

    # Automatically detect whether the vector is categorical or numeric
    if (all(control_stress_vector %in% c("Control", "Stress"))) {
      control_values = trait_values[control_stress_vector == "Control"]
      stress_values = trait_values[control_stress_vector == "Stress"]
    } else if (all(control_stress_vector %in% c(0, 1))) {
      control_values = trait_values[control_stress_vector == 0]
      stress_values = trait_values[control_stress_vector == 1]
    } else {
      stop("control_stress_vector must contain 'Control'/'Stress' labels or numeric 0/1 values.")
    }

    if (length(control_values) != length(stress_values)) {
      stop("Control and stress groups must have the same number of values.")
    }

    # Compute PD as mean absolute difference
    return(mean(abs(stress_values - control_values), na.rm = TRUE))
  }

  # If no control vs. stress is provided, use the selected method
  if (method == "pairwise") {
    pd_values = c()
    env_combinations = combn(unique_envs, 2)
    for (i in 1:ncol(env_combinations)) {
      env1 = env_combinations[1, i]
      env2 = env_combinations[2, i]

      values1 = trait_values[env_values == env1]
      values2 = trait_values[env_values == env2]

      if (length(values1) != length(values2)) next  # Skip if unequal sample sizes

      pd_values = c(pd_values, mean(abs(values1 - values2), na.rm = TRUE))
    }
    return(mean(pd_values, na.rm = TRUE))

  } else if (method == "reference") {
    env_means = tapply(trait_values, env_values, mean, na.rm = TRUE)
    reference_env = names(which.min(env_means))

    pd_values = c()
    for (env in unique_envs) {
      if (env == reference_env) next
      values1 = trait_values[env_values == reference_env]
      values2 = trait_values[env_values == env]
      if (length(values1) != length(values2)) next  # Skip if unequal sample sizes
      pd_values = c(pd_values, mean(abs(values1 - values2), na.rm = TRUE))
    }
    return(mean(pd_values, na.rm = TRUE))

  } else if (method == "variability") {
    return(max(trait_values, na.rm = TRUE) - min(trait_values, na.rm = TRUE))

  } else {
    stop("Invalid method. Choose 'pairwise', 'reference', or 'variability'.")
  }
}


## test - passed on synthetic dataset

#calculate_PD(df_test6,1,c(2,3),3,1)
##########################


#' @title Calculate Fitness Plasticity Index (FPI)
#' @description This function calculates the Fitness Plasticity Index (FPI) for a given trait
#' across control and stress conditions. If no control-stress mapping is provided, FPI is calculated across all environment pairs.
#'
#' @param trait_values A numeric vector containing trait values measured across different environments.
#' @param env_values (Optional) A numeric vector specifying the environment associated with each trait measurement.
#'                   If NULL, environments are assumed to be equidistant.
#' @param control_stress (Optional) A binary vector (0/1 or "Control"/"Stress") specifying which measurements belong to control and stress.
#'                       If NULL, pairwise FPI is calculated across all environments.
#'
#' @return A numeric value (mean FPI) if control-stress mapping is provided, otherwise a named vector of mean FPI values across all pairs.
#'
#' @examples
#' trait_values = c(10, 12, 15, 17, 20, 22, 25, 27, 30, 32)
#'
#' # Compute FPI across all environment pairs with equidistant environments
#' result = calculate_FPI(trait_values)
#' print(result)
#'
#' # Specify explicit environments
#' env_values = seq_along(trait_values)
#' control_stress = c(0, 0, 1, 1, 1, 1, 0, 0, 1, 1) # 0 = Control, 1 = Stress
#' result_with_control = calculate_FPI(trait_values, env_values, control_stress)
#' print(result_with_control)
#'
#' @export
calculate_FPI = function(trait_values, env_values = NULL, control_stress = NULL) {
  # Ensure trait_values is numeric
  if (!is.numeric(trait_values)) stop("trait_values must be a numeric vector.")

  n = length(trait_values)

  # Assume equidistant environments if not provided
  if (is.null(env_values)) {
    env_values = seq_len(n)
  }

  # Ensure env_values has correct length
  if (length(env_values) != n) stop("env_values must have the same length as trait_values.")

  unique_envs = unique(env_values)

  # If only one environment is present, return NA
  if (length(unique_envs) < 2) {
    warning("Only one unique environment present. FPI cannot be calculated.")
    return(NA)
  }

  # If control-stress mapping is provided
  if (!is.null(control_stress)) {
    if (length(control_stress) != n) stop("control_stress must have the same length as trait_values.")

    # Ensure binary format (0/1 or Control/Stress)
    if (!all(control_stress %in% c(0, 1, "Control", "Stress"))) {
      stop("control_stress must be a vector of 0/1 or 'Control'/'Stress'.")
    }

    # Convert to 0/1 if needed
    if (is.character(control_stress)) {
      control_stress = ifelse(control_stress == "Control", 0, 1)
    }

    control_values = trait_values[control_stress == 0]
    stress_values = trait_values[control_stress == 1]

    # Ensure equal sample size
    min_len = min(length(control_values), length(stress_values))
    control_values = control_values[seq_len(min_len)]
    stress_values = stress_values[seq_len(min_len)]

    # Compute FPI
    FPI_values = (stress_values - control_values) / control_values
    return(mean(FPI_values, na.rm = TRUE))
  }

  # Otherwise, compute FPI for all environment pairs
  env_combinations = combn(unique_envs, 2, simplify = TRUE)
  FPI_results = numeric(ncol(env_combinations))

  for (i in seq_len(ncol(env_combinations))) {
    env1 = env_combinations[1, i]
    env2 = env_combinations[2, i]

    idx1 = which(env_values == env1)
    idx2 = which(env_values == env2)

    min_len = min(length(idx1), length(idx2))
    idx1 = idx1[seq_len(min_len)]
    idx2 = idx2[seq_len(min_len)]

    FPI_values = (trait_values[idx2] - trait_values[idx1]) / trait_values[idx1]
    FPI_results[i] = mean(FPI_values, na.rm = TRUE)
  }

  FPI_grand=mean(FPI_results)
  names(FPI_results) = apply(env_combinations, 2, paste, collapse = "-")
  return(FPI_grand)
}



## test - passed on synthetic dataset

#calculate_FPI(df_test6,1,c(2,3),1,2)


#################################


#' @title Calculate Transplant Plasticity Score (TPS)
#' @description This function calculates the Transplant Plasticity Score (TPS) for a given trait when an organism is transplanted to a different environment.
#' The TPS quantifies the relative change in a trait between the native environment and the transplanted environment.
#'
#' @param trait_values A numeric vector containing trait values measured across different environments.
#' @param env_values A numeric vector specifying the environment associated with each trait measurement.
#'                   If NULL, environments are assumed to be equidistant.
#' @param native_env A value indicating the native environment. Must be present in `env_values`.
#' @param transplanted_env A value indicating the transplanted environment. Must be present in `env_values`.
#'
#' @return A numeric value representing the mean TPS for the trait.
#' @examples
#' trait_values = c(50, 60, 70, 80)
#' env_values = c(1, 1, 2, 2)  # Native = 1, Transplanted = 2
#'
#' # Calculate TPS between native (1) and transplanted (2)
#' tps_result = calculate_TPS(trait_values, env_values, native_env = 1, transplanted_env = 2)
#' print(tps_result)
#'
#' @export
calculate_TPS = function(trait_values, env_values, native_env, transplanted_env) {
  # Ensure trait_values is numeric
  if (!is.numeric(trait_values)) stop("trait_values must be a numeric vector.")

  # Ensure env_values is provided and has correct length
  if (is.null(env_values)) stop("env_values must be provided.")
  if (length(env_values) != length(trait_values)) stop("env_values must have the same length as trait_values.")

  # Check if specified native and transplanted environments exist in env_values
  unique_envs = unique(env_values)
  if (!(native_env %in% unique_envs) || !(transplanted_env %in% unique_envs)) {
    stop("native_env and transplanted_env must be in env_values.")
  }

  # Extract trait values for native and transplanted environments
  native_data = trait_values[env_values == native_env]
  transplanted_data = trait_values[env_values == transplanted_env]

  # Ensure equal sample size
  min_len = min(length(native_data), length(transplanted_data))
  native_data = native_data[seq_len(min_len)]
  transplanted_data = transplanted_data[seq_len(min_len)]

  # Compute TPS
  TPS_values = (transplanted_data - native_data) / native_data
  return(mean(TPS_values, na.rm = TRUE))
}
