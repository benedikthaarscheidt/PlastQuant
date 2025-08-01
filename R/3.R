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
  # Validate trait_values
  if (!is.numeric(trait_values) || !is.vector(trait_values)) {
    stop("`trait_values` must be a numeric vector.")
  }
  # Drop NA trait observations
  keep <- !is.na(trait_values)
  if (sum(keep) < 2) {
    stop("At least two non-NA trait values are required to compute PSI.")
  }
  trait <- trait_values[keep]
  n     <- length(trait)

  # Handle env_values
  if (is.null(env_values)) {
    env <- seq_len(n)
  } else {
    if (!is.numeric(env_values)) {
      if (is.factor(env_values)) {
        env_vals <- as.numeric(levels(env_values))[as.integer(env_values)]
      } else if (is.character(env_values)) {
        env_vals <- suppressWarnings(as.numeric(env_values))
        if (any(is.na(env_vals) & !is.na(env_values))) {
          stop("Cannot coerce character `env_values` to numeric.")
        }
      } else {
        stop("`env_values` must be numeric, factor, or character.")
      }
      env <- env_vals[keep]
    } else {
      if (length(env_values) != length(trait_values)) {
        stop("`env_values` must have the same length as `trait_values`.")
      }
      env <- env_values[keep]
    }
  }

  # Check environmental range
  if (abs(diff(range(env, na.rm = TRUE))) < .Machine$double.eps) {
    warning("Environmental values have zero or near-zero range; PSI cannot be calculated.")
    return(NA_real_)
  }

  # Fit linear model: trait ~ env
  fit <- lm(trait ~ env)
  coef_env <- unname(coef(fit)[ "env" ])
  if (is.na(coef_env)) {
    warning("Model slope is NA; PSI cannot be calculated.")
    return(NA_real_)
  }

  # Compute stability score
  1 / (1 + abs(coef_env))
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
calculate_RPI <- function(trait_values, env_values = NULL, env1 = NULL, env2 = NULL) {
  # Validate inputs
  if (!is.numeric(trait_values) || !is.vector(trait_values)) {
    stop("`trait_values` must be a numeric vector.")
  }
  if (!is.null(env_values) && length(env_values) != length(trait_values)) {
    stop("`env_values` must have the same length as `trait_values`.")
  }

  # Drop NA pairs
  if (is.null(env_values)) {
    keep <- !is.na(trait_values)
  } else {
    keep <- !is.na(trait_values) & !is.na(env_values)
    env_values <- env_values[keep]
  }
  trait_values <- trait_values[keep]
  n <- length(trait_values)
  if (n < 2) {
    warning("Not enough non-NA observations; RPI cannot be calculated.")
    return(NA_real_)
  }

  # Prepare environment factor
  if (is.null(env_values)) {
    env_f <- factor(seq_len(n))
  } else {
    if (!(is.numeric(env_values) || is.factor(env_values) || is.character(env_values))) {
      stop("`env_values` must be numeric, factor, or character.")
    }
    env_f <- factor(env_values)
  }
  levels_env <- levels(env_f)
  if (length(levels_env) < 2) {
    warning("Less than two unique environments; RPI cannot be calculated.")
    return(NA_real_)
  }

  # Helper to compute RPI for one pair
  pair_rpi <- function(e1, e2) {
    idx1 <- which(env_f == e1)
    idx2 <- which(env_f == e2)
    min_len <- min(length(idx1), length(idx2))
    if (min_len == 0) stop("No observations found for specified environments.")
    idx1 <- idx1[seq_len(min_len)]
    idx2 <- idx2[seq_len(min_len)]
    denom <- trait_values[idx1] + trait_values[idx2]
    r <- ifelse(denom == 0, NA_real_, abs(trait_values[idx1] - trait_values[idx2]) / denom)
    mean(r, na.rm = TRUE)
  }

  # Specific pair
  if (!is.null(env1) || !is.null(env2)) {
    if (is.null(env1) || is.null(env2)) {
      stop("Both `env1` and `env2` must be provided together.")
    }
    if (!(env1 %in% levels_env) || !(env2 %in% levels_env)) {
      stop("`env1` and `env2` must exist in `env_values`.")
    }
    return(pair_rpi(env1, env2))
  }

  # All pairs
  combos <- combn(levels_env, 2, simplify = FALSE)
  rpi_vals <- vapply(combos, function(pr) pair_rpi(pr[1], pr[2]), numeric(1))
  mean(rpi_vals, na.rm = TRUE)
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
calculate_PQ <- function(trait_values, env_values = NULL) {
  # Validate trait_values
  if (!is.numeric(trait_values) || !is.vector(trait_values)) {
    stop("`trait_values` must be a numeric vector.")
  }
  # Drop NA trait observations
  keep <- !is.na(trait_values)
  if (sum(keep) < 2) {
    stop("At least two non-NA trait values are required to compute PQ.")
  }
  trait <- trait_values[keep]

  n <- length(trait)
  # Handle env_values
  if (is.null(env_values)) {
    env <- seq_len(n)
  } else {
    if (! (is.numeric(env_values) || is.factor(env_values) || is.character(env_values))) {
      stop("`env_values` must be numeric, factor, or character.")
    }
    if (length(env_values) != length(trait_values)) {
      stop("`env_values` must have the same length as `trait_values`.")
    }
    # drop corresponding env NA
    env_raw <- env_values[keep]
    if (is.factor(env_raw)) {
      env <- as.numeric(levels(env_raw))[as.integer(env_raw)]
    } else if (is.character(env_raw)) {
      env <- suppressWarnings(as.numeric(env_raw))
      if (any(is.na(env) & !is.na(env_raw))) stop("Cannot coerce character `env_values` to numeric.")
    } else {
      env <- env_raw
    }
  }
  # Compute trait range and env range
  trait_range <- max(trait, na.rm=TRUE) - min(trait, na.rm=TRUE)
  env_range   <- max(env, na.rm=TRUE) - min(env, na.rm=TRUE)
  if (env_range == 0) {
    warning("Environmental range is zero; PQ cannot be calculated.")
    return(NA_real_)
  }
  trait_range / env_range
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
calculate_PR <- function(trait_values, env_values = NULL, across = TRUE) {
  # Coerce logical trait_values to numeric
  if (!is.numeric(trait_values)) {
    if (is.logical(trait_values)) {
      trait_values <- as.numeric(trait_values)
    } else {
      stop("`trait_values` must be a numeric vector.")
    }
  }
  if (!is.vector(trait_values)) {
    stop("`trait_values` must be a numeric vector.")
  }
  # Validate across
  if (!is.logical(across) || length(across) != 1 || is.na(across)) {
    stop("`across` must be a single TRUE or FALSE value.")
  }

  # Validate env_values length
  n <- length(trait_values)
  if (!is.null(env_values) && length(env_values) != n) {
    stop("`env_values` must have the same length as `trait_values`.")
  }

  # Drop NA observations
  if (is.null(env_values)) {
    keep <- !is.na(trait_values)
  } else {
    keep <- !is.na(trait_values) & !is.na(env_values)
    env_values <- env_values[keep]
  }
  trait <- trait_values[keep]
  if (length(trait) < 2) {
    warning("Not enough non-NA observations; PR cannot be calculated.")
    return(NA_real_)
  }

  # Global range branch
  if (across) {
    return(max(trait, na.rm=TRUE) - min(trait, na.rm=TRUE))
  }

  # Within-environment branch
  # Build environment factor
  if (is.null(env_values)) {
    env_f <- factor(seq_along(trait))
  } else if (is.factor(env_values) || is.numeric(env_values) || is.character(env_values)) {
    env_f <- factor(env_values)
  } else {
    stop("`env_values` must be numeric, factor, or character.")
  }
  levels_env <- levels(env_f)
  if (length(levels_env) == 0) {
    warning("No valid environments; PR cannot be calculated.")
    return(NA_real_)
  }

  # Compute PR per environment
  pr_vals <- vapply(levels_env, function(e) {
    vals <- trait[env_f == e]
    if (length(vals) < 1) return(NA_real_)
    max(vals, na.rm=TRUE) - min(vals, na.rm=TRUE)
  }, numeric(1))
  names(pr_vals) <- levels_env
  pr_vals
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
calculate_NRW <- function(trait_values, env_values = NULL, group_values = NULL, across = FALSE) {
  # Validate trait_values
  if (missing(trait_values) || !is.numeric(trait_values)) {
    stop("'trait_values' must be a numeric vector.")
  }
  n <- length(trait_values)
  if (n < 1) stop("'trait_values' must have at least one element.")

  # Validate env_values or assign default
  if (is.null(env_values)) {
    env_values <- seq_len(n)
  } else {
    if (!(is.numeric(env_values) || is.factor(env_values) || is.character(env_values))) {
      stop("'env_values' must be numeric, factor, or character.")
    }
    if (length(env_values) != n) {
      stop("'env_values' and 'trait_values' must be the same length.")
    }
  }

  # Validate across flag
  if (!is.logical(across) || length(across) != 1 || is.na(across)) {
    stop("'across' must be a single TRUE or FALSE.")
  }

  # Validate group_values
  if (!is.null(group_values)) {
    if (length(group_values) != n) {
      stop("'group_values' must be the same length as 'trait_values'.")
    }
  }

  # Filter out entries where trait or env is NA
  valid_idx <- !is.na(trait_values) & !is.na(env_values)
  if (!any(valid_idx)) {
    warning("No valid data points after removing NAs; returning NA.")
    return(NA_real_)
  }
  trait <- trait_values[valid_idx]
  env   <- env_values[valid_idx]
  grp   <- if (!is.null(group_values)) group_values[valid_idx] else NULL

  # Global range across all values
  if (across) {
    if (all(is.na(trait))) {
      warning("All trait_values are NA; returning NA.")
      return(NA_real_)
    }
    return(max(trait, na.rm = TRUE) - min(trait, na.rm = TRUE))
  }

  # Range per group if provided
  if (!is.null(grp)) {
    groups <- unique(grp)
    result <- setNames(numeric(length(groups)), groups)
    for (g in groups) {
      mask <- grp == g
      traits_g <- trait[mask]
      envs_g   <- env[mask]
      if (all(is.na(traits_g))) {
        warning(sprintf("All trait_values are NA for group %s; returning NA for this group.", g))
        result[g] <- NA_real_
        next
      }
      means <- tapply(traits_g, envs_g, mean, na.rm = TRUE)
      if (length(means) < 1 || all(is.na(means))) {
        warning(sprintf("No valid means for group %s; returning NA.", g))
        result[g] <- NA_real_
        next
      }
      result[g] <- max(means, na.rm = TRUE) - min(means, na.rm = TRUE)
    }
    return(result)
  }

  # Default: per-environment range
  means <- tapply(trait, env, mean, na.rm = TRUE)
  if (all(is.na(means))) {
    warning("No valid means per environment; returning NA.")
    return(NA_real_)
  }
  max(means, na.rm = TRUE) - min(means, na.rm = TRUE)
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
calculate_ESP <- function(trait_values, env_values = NULL, env_subset = NULL) {
  # Validate trait_values
  if (missing(trait_values) || !is.numeric(trait_values)) {
    stop("'trait_values' must be a numeric vector.")
  }
  n <- length(trait_values)
  if (n < 1) stop("'trait_values' must have at least one element.")

  # Validate env_values or assign default
  if (is.null(env_values)) {
    env_values <- seq_len(n)
  } else {
    if (!(is.numeric(env_values) || is.factor(env_values) || is.character(env_values))) {
      stop("'env_values' must be numeric, factor, or character.")
    }
    if (length(env_values) != n) {
      stop("'env_values' and 'trait_values' must be the same length.")
    }
  }

  # Filter out entries where trait or env is NA
  valid_idx <- !is.na(trait_values) & !is.na(env_values)
  if (!any(valid_idx)) {
    warning("No valid data points after removing NAs; returning NA.")
    return(NA_real_)
  }
  trait <- trait_values[valid_idx]
  env   <- env_values[valid_idx]

  # Determine environments to use
  unique_envs <- unique(env)
  if (!is.null(env_subset)) {
    if (!(is.numeric(env_subset) || is.factor(env_subset) || is.character(env_subset))) {
      stop("'env_subset' must be numeric, factor, or character.")
    }
    selected <- unique_envs %in% env_subset
    if (!any(selected)) {
      warning("None of the 'env_subset' values match 'env_values'; returning NA.")
      return(NA_real_)
    }
    unique_envs <- unique_envs[selected]
  }

  # Compute overall mean
  if (all(is.na(trait))) {
    warning("All 'trait_values' are NA; returning NA.")
    return(NA_real_)
  }
  mean_all <- mean(trait, na.rm = TRUE)

  # Compute ESP per environment
  ESP_values <- setNames(numeric(length(unique_envs)), unique_envs)
  for (e in unique_envs) {
    subset_vals <- trait[env == e]
    if (all(is.na(subset_vals))) {
      warning(sprintf("All 'trait_values' are NA for environment %s; assigning NA.", e))
      ESP_values[as.character(e)] <- NA_real_
    } else {
      mean_env <- mean(subset_vals, na.rm = TRUE)
      ESP_values[as.character(e)] <- (mean_env - mean_all) / mean_all
    }
  }

  # Aggregate
  abs_vals <- abs(ESP_values)
  if (all(is.na(abs_vals))) {
    warning("No valid ESP values; returning NA.")
    return(NA_real_)
  }
  if (any(is.na(abs_vals))) {
    warning("Some environments had NA ESP and were excluded from the sum.")
  }
  ESP_grand <- sum(abs_vals, na.rm = TRUE)
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
calculate_general_PD <- function(
    trait_values,
    env_values = NULL,
    control_stress_vector = NULL,
    method = "pairwise"
) {
  # Validate trait_values
  if (!is.numeric(trait_values)) {
    stop("'trait_values' must be numeric.")
  }
  n <- length(trait_values)
  if (n < 2) {
    stop("At least two trait values are required.")
  }

  # Setup env_values (default equidistant sequence)
  if (is.null(env_values)) {
    env_values <- seq_len(n)
  }
  if (length(env_values) != n) {
    stop("'env_values' must be the same length as 'trait_values'.")
  }
  env <- as.factor(env_values)

  # Validate control_stress_vector length
  if (!is.null(control_stress_vector) && length(control_stress_vector) != n) {
    stop("'control_stress_vector' must be the same length as 'trait_values'.")
  }

  # Filter NAs
  valid <- !is.na(trait_values) & !is.na(env)
  if (sum(valid) < 2) {
    warning("Not enough valid data points; returning NA.")
    return(NA_real_)
  }
  trait <- trait_values[valid]
  env   <- as.factor(env_values[valid])

  # Control vs Stress method
  if (!is.null(control_stress_vector)) {
    csv <- control_stress_vector[valid]
    if (all(csv %in% c(0,1))) {
      ctrl <- trait[csv == 0]
      strs <- trait[csv == 1]
    } else if (all(csv %in% c("Control", "Stress"))) {
      ctrl <- trait[csv == "Control"]
      strs <- trait[csv == "Stress"]
    } else {
      stop("'control_stress_vector' must be 0/1 or 'Control'/'Stress'.")
    }
    if (length(ctrl) != length(strs)) {
      stop("Control and Stress groups must have equal length.")
    }
    return(mean(abs(strs - ctrl), na.rm = TRUE))
  }

  # Validate method
  if (!method %in% c("pairwise", "reference", "variability")) {
    stop("'method' must be 'pairwise', 'reference', or 'variability'.")
  }

  # Variability method: max - min
  if (method == "variability") {
    return(max(trait, na.rm = TRUE) - min(trait, na.rm = TRUE))
  }

  # Reference method: compare environment means to lowest-mean environment
  if (method == "reference") {
    env_means <- tapply(trait, env, mean, na.rm = TRUE)
    ref_env <- names(which.min(env_means))
    diffs <- abs(env_means - env_means[ref_env])
    diffs <- diffs[names(diffs) != ref_env]
    return(mean(diffs, na.rm = TRUE))
  }

  # Pairwise method: average PD across all env pairs
  if (method == "pairwise") {
    levels_env <- levels(env)
    combos <- utils::combn(levels_env, 2, simplify = FALSE)
    pd_vals <- numeric(0)
    for (pair in combos) {
      v1 <- trait[env == pair[1]]
      v2 <- trait[env == pair[2]]
      v1 <- v1[!is.na(v1)]; v2 <- v2[!is.na(v2)]
      if (length(v1) < 1 || length(v2) < 1) {
        warning(sprintf("No valid data for '%s' vs '%s'; skipped.", pair[1], pair[2]))
        next
      }
      if (length(v1) != length(v2)) {
        warning(sprintf("Unequal sample sizes for '%s' vs '%s'; pairing with minimum length.", pair[1], pair[2]))
        m <- min(length(v1), length(v2))
        v1 <- head(v1, m); v2 <- head(v2, m)
      }
      pd_vals <- c(pd_vals, mean(abs(v1 - v2), na.rm = TRUE))
    }
    if (length(pd_vals) == 0) {
      warning("No PD values computed; returning NA.")
      return(NA_real_)
    }
    return(mean(pd_vals, na.rm = TRUE))
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
calculate_FPI <- function(trait_values,
                          env_values     = NULL,
                          control_stress = NULL,
                          summary_fun    = mean,
                          na.rm          = TRUE,
                          drop_zero      = TRUE) {

  # 1. Basic checks
  if (!is.numeric(trait_values)) {
    stop("`trait_values` must be a numeric vector.")
  }
  n <- length(trait_values)
  if (n < 2L) {
    stop("Need at least two trait values to compute FPI.")
  }

  # 2. env_values default & check
  if (is.null(env_values)) {
    env_values <- seq_len(n)
  }
  if (length(env_values) != n) {
    stop("`env_values` must have same length as `trait_values`.")
  }
  # allow factors/chars → force to factor
  if (is.factor(env_values) || is.character(env_values)) {
    env_values <- as.factor(env_values)
  } else if (!is.numeric(env_values) && !is.integer(env_values)) {
    stop("`env_values` must be numeric, integer, factor, or character.")
  }

  # 3. control_stress branch
  if (!is.null(control_stress)) {
    if (length(control_stress) != n) {
      stop("`control_stress` must have same length as `trait_values`.")
    }
    # logical → 0/1
    if (is.logical(control_stress)) {
      control_stress <- ifelse(control_stress, 1L, 0L)
    }
    # factor/char → "Control"/"Stress"
    if (is.factor(control_stress) || is.character(control_stress)) {
      cs_c <- as.character(control_stress)
      if (!all(cs_c %in% c("Control", "Stress"))) {
        stop("If character/factor, `control_stress` must be 'Control'/'Stress'.")
      }
      control_stress <- ifelse(cs_c == "Control", 0L, 1L)
    }
    if (!all(control_stress %in% c(0L, 1L))) {
      stop("`control_stress` must be logical, 0/1, or 'Control'/'Stress'.")
    }

    # split
    control_idx <- which(control_stress == 0L)
    stress_idx  <- which(control_stress == 1L)
    if (length(control_idx) < 1 || length(stress_idx) < 1) {
      stop("Need at least one Control and one Stress observation.")
    }

    # drop zeros (ignore NAs)
    if (drop_zero && any(trait_values[control_idx] == 0, na.rm = TRUE)) {
      zeros <- sum(trait_values[control_idx] == 0, na.rm = TRUE)
      unit  <- if (zeros == 1) "value" else "values"
      warning(sprintf(
        "Dropping %d zero control %s to avoid division by zero.",
        zeros, unit
      ))
      keep        <- !is.na(trait_values[control_idx]) & trait_values[control_idx] != 0
      control_idx <- control_idx[keep]
    }

    # pair up
    m            <- min(length(control_idx), length(stress_idx))
    control_idx  <- control_idx[seq_len(m)]
    stress_idx   <- stress_idx[seq_len(m)]

    FPI_vec  <- (trait_values[stress_idx] - trait_values[control_idx]) /
      trait_values[control_idx]
    pairwise <- setNames(FPI_vec,
                         paste0("Control–Stress_", seq_along(FPI_vec)))
    overall  <- summary_fun(FPI_vec, na.rm = na.rm)

    return(list(pairwise = pairwise, overall = overall))
  }

  # 4. Pairwise env/env
  uniq_envs <- unique(env_values)
  if (length(uniq_envs) < 2L) {
    warning("Only one unique environment; returning NA.")
    return(list(pairwise = NA_real_, overall = NA_real_))
  }

  combos  <- utils::combn(uniq_envs, 2L, simplify = FALSE)
  out_vec <- numeric(length(combos))
  names(out_vec) <- vapply(combos,
                           function(x) paste0(x[1], "–", x[2]),
                           character(1))

  for (nm in names(out_vec)) {
    evs  <- strsplit(nm, "–", fixed = TRUE)[[1]]
    idx1 <- which(env_values == evs[1])
    idx2 <- which(env_values == evs[2])

    m    <- min(length(idx1), length(idx2))
    idx1 <- idx1[seq_len(m)]
    idx2 <- idx2[seq_len(m)]

    if (drop_zero && any(trait_values[idx1] == 0, na.rm = TRUE)) {
      zeros <- sum(trait_values[idx1] == 0, na.rm = TRUE)
      unit  <- if (zeros == 1) "value" else "values"
      warning(sprintf(
        "Dropping %d zero baseline %s for pair %s",
        zeros, unit, nm
      ))
      keep <- !is.na(trait_values[idx1]) & trait_values[idx1] != 0
      idx1 <- idx1[keep]
      idx2 <- idx2[keep]
    }

    vals     <- (trait_values[idx2] - trait_values[idx1]) / trait_values[idx1]
    out_vec[nm] <- unname(summary_fun(vals, na.rm = na.rm))
  }

  overall <- summary_fun(out_vec, na.rm = na.rm)
  return(list(pairwise = out_vec, overall = overall))
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
calculate_TPS <- function(trait_values,
                          env_values,
                          native_env,
                          transplanted_env,
                          summary_fun = mean,
                          na.rm       = TRUE,
                          drop_zero   = TRUE) {

  # 1. Basic checks
  if (!is.numeric(trait_values)) {
    stop("`trait_values` must be a numeric vector.")
  }
  n <- length(trait_values)
  if (n < 1L) {
    stop("`trait_values` must have length >= 1.")
  }
  if (missing(env_values) || is.null(env_values)) {
    stop("`env_values` must be provided.")
  }
  if (length(env_values) != n) {
    stop("`env_values` must have the same length as `trait_values`.")
  }
  # allow factor/char
  if (is.factor(env_values) || is.character(env_values)) {
    env_values <- as.factor(env_values)
  } else if (!is.numeric(env_values) && !is.integer(env_values)) {
    stop("`env_values` must be numeric, integer, factor, or character.")
  }

  # 2. native/transplanted must exist
  uniq <- unique(env_values)
  if (!(native_env %in% uniq)) {
    stop("`native_env` not found in `env_values`.")
  }
  if (!(transplanted_env %in% uniq)) {
    stop("`transplanted_env` not found in `env_values`.")
  }

  # 3. Pull out and pair
  native_idx      <- which(env_values == native_env)
  transplanted_idx<- which(env_values == transplanted_env)
  if (length(native_idx) < 1 || length(transplanted_idx) < 1) {
    stop("Need at least one observation in each env.")
  }
  m               <- min(length(native_idx), length(transplanted_idx))
  native_data     <- trait_values[native_idx[seq_len(m)]]
  transplanted_data <- trait_values[transplanted_idx[seq_len(m)]]

  # 4. Drop zeros if requested
  if (drop_zero && any(native_data == 0, na.rm = TRUE)) {
    zeros <- sum(native_data == 0, na.rm = TRUE)
    unit  <- if (zeros == 1) "value" else "values"
    warning(sprintf(
      "Dropping %d zero native %s to avoid division by zero.",
      zeros, unit
    ))
    keep         <- !is.na(native_data) & native_data != 0
    native_data  <- native_data[keep]
    transplanted_data <- transplanted_data[keep]
  }

  # 5. If nothing left, warn & return NA
  if (length(native_data) < 1) {
    warning("No valid pairs remain after zero-dropping; returning NA.")
    return(list(raw = numeric(0), overall = NA_real_))
  }

  # 6. Compute
  tps      <- (transplanted_data - native_data) / native_data
  raw      <- unname(tps)
  overall  <- summary_fun(raw, na.rm = na.rm)

  list(raw = raw, overall = overall)
}
