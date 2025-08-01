#this files contains the following indices/functions:
#generate_synthetic_data,
#check_and_install_packages,
#RDPI	(rdpi_calculation) - tested,
#RDPIs (rdpi_mean_calculation) - tested,
#ESPI (calculate_ESPI) - tested,
#ESPIid (espiid_calculation) - tested,
#evwpi_calculation (idea from Benedikt)


##### datasets for testing
#
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
############################




#' @title Calculate Relative Distance Plasticity Index (RDPI)
#'
#' @description
#' Calculates the Relative Distance Plasticity Index (RDPI) for trait values across different
#' environmental conditions. The calculation follows Valladares et al. (2006) methodology.
#'
#' @param trait_values Numeric vector containing trait measurements (already averaged across replicates)
#' @param env_values Optional vector of environmental conditions. If NULL, equidistant steps are assumed
#'
#' @return Numeric value representing the RDPI score
#'
#' @examples
#' # With equidistant environments
#' traits = c(10, 12, 15, 18, 20)
#' rdpi = calculate_rdpi(traits)
#'
#' # With specified environments
#' traits = c(10, 12, 15, 18, 20)
#' envs = c(1, 2, 4, 6, 8)  # Non-equidistant environments
#' rdpi = calculate_rdpi(traits, envs)
#'
#' @export
calculate_rdpi <- function(trait_values, env_values = NULL) {
  # Validate trait_values
  if (!is.numeric(trait_values)) {
    stop("trait_values must be a numeric vector.")
  }
  n <- length(trait_values)
  if (n < 2) {
    stop("trait_values must contain at least two observations.")
  }

  # Handle env_values
  if (is.null(env_values)) {
    env_values <- seq_len(n)
  }
  if (length(env_values) != n) {
    stop("trait_values and env_values must have the same length.")
  }
  env <- factor(env_values)
  if (nlevels(env) < 2) {
    stop("At least two environment levels are required to calculate RDPI.")
  }

  # Compute environment means
  means <- tapply(trait_values, env, function(x) {
    if (all(is.na(x))) {
      NA_real_
    } else {
      mean(x, na.rm = TRUE)
    }
  })
  # Drop NA means
  means <- means[!is.na(means)]
  if (length(means) < 2) {
    warning("Not enough valid environment means; RDPI cannot be calculated.")
    return(NA_real_)
  }

  # Generate all unique pairs
  env_levels <- names(means)
  pairs <- utils::combn(env_levels, 2, simplify = FALSE)

  # Compute RDPI per pair
  rdpi_vals <- vapply(pairs, function(pr) {
    m1 <- means[ pr[1] ]
    m2 <- means[ pr[2] ]
    abs_diff <- abs(m2 - m1)
    sum_vals <- m1 + m2
    if (sum_vals == 0) {
      0
    } else {
      abs_diff / sum_vals
    }
  }, numeric(1))

  # Return mean RDPI across pairs
  mean(rdpi_vals, na.rm = TRUE)
}



## Example usage with synthetic data
#external_light = rep(c(0.4, 0.6, 0.8), each = 100)
#external_water = sample(rep(c("Low", "High"), each = 150))
#
## Calculate RDPI with factors not in the dataframe
#RDPI= rdpi_calculation(synthetic_data1, sp = NULL, trait_cols = c(3, 4), factors = 2, factors_not_in_dataframe = list(external_light, external_water))
#
#
#### test - passed on synthetic dataset and crosschecked with data from package
#
#
#df_test20 = data.frame(
#  Column0 = as.factor(c(rep(1, 20), rep(2, 20))),  # Two groups/species, 20 observations each
#  Column1 = as.factor(c(rep(1, 10), rep(2, 10), rep(1, 10), rep(2, 10))),  # Two environments per group
#  Column2 = c(rep(2, 10), rep(1, 10), rep(2, 10), rep(1, 10))  # Numeric trait values
#)
##rdpi_calculation(df_test20,sp=1, trait_cols=3,factors = 2)


#####################################################








#' @title Environmental Sensitivity Performance Index (ESPI)
#'
#' @description
#' Calculates the Environmental Sensitivity Performance Index (ESPI) for a single trait
#' across environmental conditions. ESPI quantifies the sensitivity of a trait by
#' comparing the difference between its maximum and minimum mean values to the range of the
#' environmental gradient.
#'
#' ESPI is defined as:
#' \deqn{ESPI = \frac{\max(\bar{z}_e) - \min(\bar{z}_e)}{\lvert \max(e) - \min(e)\rvert},}
#' where \(\bar{z}_e\) is the mean trait value in environment \(e\) and \(e\) is the environmental value.
#'
#' @param trait_values
#'   Numeric vector of trait measurements. At least two observations are required.
#' @param env_values
#'   Optional numeric, factor, or character vector of environmental values corresponding to `trait_values`.
#'   If `NULL`, environments are assumed equidistant (1, 2, …, length(trait_values)).
#'
#' @return
#'   Numeric scalar: ESPI value. Returns `NA` with a warning if ESPI cannot be calculated (e.g., zero range or missing data).
#'
#' @examples
#' trait <- c(10, 15, 20, 25, 30)
#' env   <- c(1, 2, 3, 4, 5)
#' calculate_ESPI(trait, env)
#'
#' # Without specifying env_values (assumes equidistant)
#' calculate_ESPI(trait)
#'
#' @export
calculate_ESPI <- function(trait_values, env_values = NULL) {
  # Validate trait_values
  if (!is.numeric(trait_values) || !is.vector(trait_values)) {
    stop("`trait_values` must be a numeric vector.")
  }
  n <- length(trait_values)
  if (n < 2) {
    stop("`trait_values` must contain at least two observations.")
  }

  # Handle env_values
  if (is.null(env_values)) {
    env_values <- seq_len(n)
  }
  if (length(env_values) != n) {
    stop("`env_values` must have the same length as `trait_values`.")
  }

  # Coerce env_values to numeric
  if (!is.numeric(env_values)) {
    if (is.factor(env_values)) {
      env_num <- as.numeric(levels(env_values))[as.integer(env_values)]
    } else if (is.character(env_values)) {
      env_num <- suppressWarnings(as.numeric(env_values))
      if (any(is.na(env_num) & !is.na(env_values))) {
        stop("Cannot coerce character `env_values` to numeric.")
      }
    } else {
      stop("`env_values` must be numeric, factor, or character.")
    }
  } else {
    env_num <- env_values
  }

  # Drop NA observations
  keep <- !is.na(trait_values) & !is.na(env_num)
  if (!any(keep)) {
    warning("All observations are NA; ESPI cannot be calculated.")
    return(NA_real_)
  }
  trait <- trait_values[keep]
  env   <- env_num[keep]

  # Compute environment means
  gm <- tapply(trait, env, mean, na.rm = TRUE)
  gm <- gm[!is.na(gm)]
  if (length(gm) < 1) {
    warning("No valid environment means; ESPI cannot be calculated.")
    return(NA_real_)
  }
  max_mean <- max(gm)
  min_mean <- min(gm)

  # Compute environmental range
  env_range <- max(env) - min(env)
  if (env_range == 0) {
    warning("Environmental range is zero; ESPI cannot be calculated.")
    return(NA_real_)
  }

  # Calculate ESPI
  (max_mean - min_mean) / env_range
}

## test - passed on synthetic dataset


#df_test2.2 = data.frame(Column0 = c(rep(1, 10), rep(2, 10)),Column1 = c(rep(2, 5),rep(4,5) ,rep(1, 10)), Column2 = c(rep(2, 10), rep(4, 10)))
#calculate_ESPI(df_test2.2,trait_cols = c(2,3),env_col = 1)


#####################################


#' @title Environmental Sensitivity Performance Index for Individual Differences (ESPIID)
#'
#' @description
#' Computes the Environmental Sensitivity Performance Index for Individual Differences (ESPIID)
#' for a single trait across an environmental gradient. ESPIID quantifies the phenotypic distance
#' between individuals in different environments, normalized by the environmental distance.
#'
#' Specifically, for each pair of environments \(i,j\):
#' \deqn{\mathrm{ESPIID}_{ij} = \frac{\mathrm{Dist}(z_i, z_j)}{|e_i - e_j|},}
#' where \(\mathrm{Dist}(z_i,z_j)\) is the mean (or median) absolute difference between all individuals
#' in environments \(i\) and \(j\), and \(e_i,e_j\) are the numeric environment indices.
#'
#' @param trait_values
#'   A numeric vector of trait measurements. Must have length ≥ 2.
#' @param use_median
#'   Logical; if `TRUE`, uses the median of all pairwise absolute differences when computing distances.
#'   If `FALSE` (default), uses the mean.
#' @param aggregate
#'   Logical; if `TRUE` (default), returns the overall ESPIID as the mean across all environment pairs.
#'   If `FALSE`, returns a named numeric vector of ESPIID values for each environment pair.
#'
#' @return
#' If `aggregate = TRUE`, a single numeric ESPIID value. If `aggregate = FALSE`, a named numeric vector
#' where names are of the form `"i-j"` for each environment pair.
#'
#' @examples
#' # Equidistant environments, mean-based, aggregated
#' vals <- c(10, 12, 15, 18, 20, 25)
#' calculate_espiid(vals)
#'
#' # Use median instead of mean
#' calculate_espiid(vals, use_median = TRUE)
#'
#' # Return all pairwise ESPIID values
#' calculate_espiid(vals, aggregate = FALSE)
#'
#' @export
calculate_espiid <- function(trait_values,
                             use_median = FALSE,
                             aggregate  = TRUE) {
  # Validate inputs
  if (!is.numeric(trait_values) || !is.vector(trait_values)) {
    stop("`trait_values` must be a numeric vector.")
  }
  if (!is.logical(use_median) || length(use_median) != 1 || is.na(use_median)) {
    stop("`use_median` must be a single TRUE or FALSE value.")
  }
  if (!is.logical(aggregate)  || length(aggregate)  != 1 || is.na(aggregate)) {
    stop("`aggregate` must be a single TRUE or FALSE value.")
  }
  # Drop NA observations
  keep <- !is.na(trait_values)
  if (sum(keep) < 2) {
    warning("At least two non-NA trait values are required; ESPIID cannot be calculated.")
    return(NA_real_)
  }
  trait <- trait_values[keep]
  n     <- length(trait)

  # Equidistant environments (1...n)
  env <- seq_len(n)

  # Prepare output
  pair_names <- combn(n, 2, function(idx) paste(idx, collapse = "-"))
  n_pairs    <- length(pair_names)
  espiid_vec <- numeric(n_pairs)
  names(espiid_vec) <- pair_names

  # Compute ESPIID for each pair
  idx <- 1L
  for (i in seq_len(n - 1L)) {
    for (j in (i + 1L):n) {
      diff_val <- abs(trait[j] - trait[i])
      # for single observations, mean == median == diff_val
      summary_diff <- if (use_median) {
        diff_val
      } else {
        diff_val
      }
      env_dist <- abs(j - i)
      espiid_vec[idx] <- if (env_dist > 0) summary_diff / env_dist else NA_real_
      idx <- idx + 1L
    }
  }

  # Aggregate or return full vector
  if (aggregate) {
    mean(espiid_vec, na.rm = TRUE)
  } else {
    espiid_vec
  }
}

#l=list(external_light)
#
#external_factors=list(external_light,external_water)
#
### test - passed on synthetic dataset
#
#espiid_calculation(df_test2.2, trait_cols=2,factors=1)




