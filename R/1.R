#this files contains the following indices/functions:
#generate_synthetic_data,
#impute,
#coefficient-of variation total (calculate_CVt) - tested,
#slope of norm reaction (calculate_reaction_norm_slope) - tested,
#slope of plastic response (D) (calculate_D_slope)- tested,
#response coefficient (RC) (calculate_RC)- tested,
#Standard deviation of means (CVm) (calculate_CVm)- tested,
#Standard deviation of medians (CVmd)(calculate_CVmd)- tested,
#Grand plasticity (calculate_grand_plasticity)- tested,
#combine_factors- tested,
#Phenotypic Plasticity Index (calculate_PPF)- tested,
#Phenotypic Plasticity Index (calculate_Phenotypic_Plasticity_Index)- tested,
#PImd (calculate_PImd)- tested,
#PILSM (calculate_PILSM)- tested,
#RTR (calculate_RTR)- tested,
#PIR (calculate_PIR) - tested


################################

#' @import ggplot2
#' @import dplyr
#' @import purrr
#' @import MASS
#' @import emmeans
#' @import gridExtra
#' @import stats
#' @importFrom roxygen2
#' @importFrom reshape2
NULL



################################




#' @title Coefficient of Variation for Plasticity (CVt)
#'
#' @description
#' Calculates the Coefficient of Variation (CVt) for a single numeric vector of trait values.
#' Returns `NA` if fewer than two observations, or throws an error if there are any `NA`s.
#'
#' @param trait_values A numeric vector containing the trait values for a single genotype.
#' @return A numeric value representing CVt (sd/mean), or `NA` if length < 2.
#'
#' @examples
#' # Example usage
#' trait_values <- c(100.2, 92, 83.5)
#' cvt <- calculate_CVt(trait_values)
#' print(cvt)
#'
#' @export
calculate_CVt = function(trait_values) {
  if (length(trait_values) < 2) {
    return(NA)  # Avoid division by zero for single values
  }

  if (any(is.na(trait_values))) {
    stop("The input vector contains missing values. Please handle missing data before calling this function.")
  }
  if(mean(trait_values)==0){
    return(NA)
  }
  # Calculate and return the CVt
  return(sd(trait_values) / mean(trait_values))
}




### test - passed in by hand calculation with this dataset

#df_test1 = data.frame(genotype=c(rep(1,10),rep(1,10)),Column1 = c(rep(4, 10), rep(2, 10)), Column2 = c(rep(10, 10), rep(1, 10)))
#
##traitwise
#cv=calculate_CVt(df_test1,genotype_col = 1, trait_col = 2)
#print(sd(df_test1[,2])/mean(df_test1[,2]))
##print(sd(df_test1[,2])/mean(df_test1[,2]))
#print(cv)
##total
#df_test1_flattened=as.numeric(unlist(df_test1))
#calculate_CVt(df_test1,traitwise = F)
#print(sd(df_test1_flattened)/mean(df_test1_flattened))



################################



#' @title Slope of a Linear Reaction Norm
#'
#' @description
#' Fits a simple linear regression of trait values on an environmental gradient
#' and returns the estimated slope coefficient. If no explicit `environments`
#' vector is provided, environments are assumed to be sequential (1, 2, …).
#'
#' @param trait_values A numeric vector of trait measurements taken across environments.
#'   Must have length ≥ 2 and contain no missing values.
#' @param environments Optional numeric vector of the same length as `trait_values`,
#'   indicating the environmental values or indices. If `NULL` (default),
#'   uses `seq_along(trait_values)` as the environmental gradient.
#'
#' @return A single numeric value: the estimated slope of the regression line
#'   (i.e., the change in trait per unit change in environment). Returns `NA`
#'   if `length(trait_values) < 2`.
#'
#' @examples
#' # Default sequential environments:
#' trait <- c(5, 6, 7, 8, 9)
#' calculate_reaction_norm_slope(trait)
#'
#' # Custom environment values:
#' env   <- c(10, 20, 30, 40, 50)
#' trait <- c(2.1, 2.8, 3.5, 4.2, 4.9)
#' calculate_reaction_norm_slope(trait, environments = env)
#'
#' # Too few observations yields NA:
#' calculate_reaction_norm_slope(42)
#'
#' @export
calculate_reaction_norm_slope = function(trait_values, environments=NULL) {
  if (length(trait_values) < 2) {
    return(NA)
  }
  if(any(is.na(trait_values))) {
    stop("The input vector contains missing values. Please handle missing data before calling this function.")
  }
  if(mean(trait_values) == 0) {
    return(NA)
  }
  if(!is.null(environments) && length(environments) != length(trait_values)) {
    stop("If environments are provided, they must match the length of trait_values.")
  }
  if(!is.null(environments) && any(is.na(environments))){
    stop("If environments are provided, they cannot have any NAs")
  }

  if(is.null(environments)){
    environments = seq_along(trait_values)
  }

  lm_fit = lm(trait_values ~ environments)

  return(unname(coef(lm_fit)["environments"]))
}


##############################

#' @title Nonlinear Reaction Norm Nonlinearity Score
#'
#' @description
#' Fits a raw‐polynomial regression of trait values on an environmental gradient
#' and returns a single nonlinearity score, defined as the sum of the absolute
#' values of all polynomial coefficients beyond the intercept. This provides an
#' index of how much the reaction norm deviates from linearity.
#'
#' @param trait_values A numeric vector of trait measurements for one genotype
#'   across environments. Must have at least `degree + 1` values and contain no `NA`.
#' @param degree An integer ≥ 1 specifying the degree of the polynomial model
#'   to fit (default: 2).
#' @param environments Optional numeric vector of the same length as `trait_values`,
#'   giving the environmental values. If `NULL`, uses sequential indices.
#'
#' @return A single numeric nonlinearity score (sum of absolute non‐intercept
#'   coefficients), or `NA` if there are insufficient data or model fails.
#'
#' @examples
#' vals <- c(10, 12, 15, 20, 25)
#' calculate_reaction_norm_non_linear(vals)
#' env <- c(0, 5, 10, 15, 20)
#' calculate_reaction_norm_non_linear(vals, degree = 3, environments = env)
#' calculate_reaction_norm_non_linear(vals, degree = 4)
#'
#' @export
calculate_reaction_norm_non_linear <- function(trait_values,
                                               degree = 2,
                                               environments = NULL) {
  # Input validation
  if (!is.numeric(trait_values)) {
    stop("`trait_values` must be a numeric vector.")
  }
  if (any(is.na(trait_values))) {
    stop("`trait_values` contains missing values; please handle NAs first.")
  }
  if (!is.numeric(degree) || length(degree) != 1 || degree < 1 || degree != as.integer(degree)) {
    stop("`degree` must be a single integer ≥ 1.")
  }
  if (!is.null(environments)) {
    if (!is.numeric(environments) || length(environments) != length(trait_values)) {
      stop("`environments` must be a numeric vector the same length as `trait_values`.")
    }
    if (any(is.na(environments))) {
      stop("`environments` contains missing values; please handle NAs first.")
    }
  }

  # Default environments if not provided
  if (is.null(environments)) {
    environments <- seq_along(trait_values)
  }

  # Ensure enough data points
  if (length(trait_values) < degree + 1) {
    warning("Not enough data points for the specified degree; returning NA.")
    return(NA_real_)
  }

  # Fit the polynomial model
  model <- tryCatch({
    lm(trait_values ~ poly(environments, degree, raw = TRUE))
  }, error = function(e) {
    warning("Model fitting failed; returning NA.")
    return(NULL)
  })
  if (is.null(model)) {
    return(NA_real_)
  }

  # Extract and sum non-intercept coefficients
  coefs <- coef(model)
  nonlinearity_score <- sum(abs(coefs[-1]))

  return(nonlinearity_score)
}



#test - passed with synthetic dataset

#df_test2 = data.frame(Column0 = c(rep(1, 10), rep(2, 10)),Column1 = c(rep(2, 10), rep(1, 10)), Column2 = c(rep(2, 10), rep(4, 10)))
#calculate_reaction_norm_slope(df_test2,env_col = 1,trait_cols = c(2,3),plot = T)

################################


#NOTE: if the resource availability is an actual measurement of a metabolite which is being used by the plant then the grouping of the plants into high vs low resource availability should be done by clustering.
# Sadly I am missing the insight into the common practices in the field.


#' @title D‑Slope Plasticity Metric
#'
#' @description
#' Computes the “D slope” plasticity metric, defined as the difference between the mean of the
#' upper‑fraction and the mean of the lower‑fraction of a sorted vector of trait values.
#' This quantifies the scope of plastic response by comparing extreme ends of the trait distribution.
#'
#' @param trait_values A numeric vector of trait values. Must contain at least one value.
#' @param lower_fraction A numeric between 0 and 1 (exclusive) indicating the fraction of the
#'   lowest sorted values to include in the “lower” group. Default is 0.2 (lowest 20%).
#' @param upper_fraction A numeric between 0 and 1 (inclusive) indicating the fraction of the
#'   highest sorted values to include in the “upper” group. Default is 0.8 (top 20%).
#'
#' @return A single numeric:
#'   \code{mean(upper_values) - mean(lower_values)}.
#'   If \code{lower_fraction * length(trait_values)} < 1 or \code{(1 - upper_fraction) * length(trait_values)} < 1,
#'   one group will be empty and its mean will be \code{NA}, yielding \code{NA}.
#'
#' @examples
#' # Default: bottom 20% vs top 20%
#' vals <- c(10, 12, 20, 22, 25, 28, 30, 32, 35, 40)
#' calculate_D_slope(vals)
#'
#' # Custom fractions: bottom 10% vs top 10%
#' calculate_D_slope(vals, lower_fraction = 0.1, upper_fraction = 0.9)
#'
#' # If you want only the extremes (5% each):
#' calculate_D_slope(vals, lower_fraction = 0.05, upper_fraction = 0.95)
#'
#' @export
calculate_D_slope = function(trait_values, lower_fraction = 0.2, upper_fraction = 0.8) {
  # Ensure trait values are sorted and valid
  sorted_values = sort(trait_values, na.last = TRUE)
  if (length(sorted_values) < 2) {
    return(NA)  # Not enough data to calculate D slope
  }
  if (any(is.na(sorted_values))) {
    stop("The input vector contains missing values. Please handle missing data before calling this function.")
  }

  if(lower_fraction<0 || lower_fraction>1 || upper_fraction<0 || upper_fraction>1) {
    stop("`lower_fraction` and `upper_fraction` must be between 0 and 1.")
  }
  if (lower_fraction > upper_fraction) {
    stop("`lower_fraction` must be less than `upper_fraction`.")
  }

  # Calculate boundaries
  lower_boundary_index = floor(length(sorted_values) * lower_fraction)
  upper_boundary_index = ceiling(length(sorted_values) * upper_fraction)

  # Extract lower and upper segments
  lower_values = sorted_values[1:lower_boundary_index]
  upper_values = sorted_values[(upper_boundary_index+1):length(sorted_values)]

  # Calculate D slope
  D_slope = mean(upper_values, na.rm = TRUE) - mean(lower_values, na.rm = TRUE)

  return(D_slope)
}

#test - passed with synthetic

#calculate_D_slope(df_test2,env_col = 1,trait_cols = 2)
#df_test3=data.frame(Column0 = c(rep(3, 10), rep(2, 10),rep(1,10)),Column1 = c(rep(2, 10), rep(1, 15),rep(3,5)), Column2 = c(rep(2, 10), rep(4, 10),rep(3,10)))
#calculate_D_slope(df_test3,env_col = 1,trait_cols = 3)
#mean(df_test3[1:15,3])-mean(df_test3[16:nrow(df_test3),3])




################################



#' @title Calculate the Response Coefficient (RC)
#' @description
#' Computes the Response Coefficient (RC), defined as the ratio of the mean of the highest‐value
#' fraction of a sorted trait vector to the mean of the lowest‐value fraction. This metric
#' quantifies the proportional change between high and low trait values.
#'
#' @param trait_values Numeric vector of trait measurements.
#'   Must be length ≥ 2 and non‐missing.
#' @param lower_fraction Numeric scalar in (0, 1). Fraction of observations to include in the “low” group.
#'   Defaults to 0.5 (the lowest 50 %).
#' @param upper_fraction Numeric scalar in (0, 1]. Fraction of observations to include in the “high” group.
#'   Defaults to 0.5 (the highest 50 %).
#'
#' @return
#' A single numeric value:
#' \deqn{\mathrm{RC} = \frac{\bar{x}_{\mathrm{high}}}{\bar{x}_{\mathrm{low}}}}
#' where \(\bar{x}_{\mathrm{high}}\) and \(\bar{x}_{\mathrm{low}}\) are the means of the top and
#' bottom fractions of the sorted trait vector, respectively.
#'
#' @examples
#' trait_values <- c(10, 12, 20, 22, 25, 28, 30, 32, 35, 40)
#'
#' # Default (50% low vs. 50% high)
#' calculate_RC(trait_values)
#'
#' # 20% low vs. 20% high
#' calculate_RC(trait_values, lower_fraction = 0.2, upper_fraction = 0.2)
#'
#' @export
calculate_RC = function(trait_values, lower_fraction = 0.5, upper_fraction = 0.5) {
  # Ensure trait values are sorted
  sorted_values = sort(trait_values, na.last = TRUE)
  if (length(sorted_values) < 2) {
    return(NA)  # Not enough data to calculate RC
  }
  if (any(is.na(sorted_values))) {
    stop("The input vector contains missing values. Please handle missing data before calling this function.")
  }
  if (lower_fraction < 0 || lower_fraction > 1 || upper_fraction < 0 || upper_fraction > 1) {
    stop("`lower_fraction` and `upper_fraction` must be between 0 and 1.")
  }
  if (lower_fraction > upper_fraction) {
    stop("`lower_fraction` must be less than `upper_fraction`.")
  }

  # Calculate indices for splitting
  lower_boundary_index = floor(length(sorted_values) * lower_fraction)
  upper_boundary_index = ceiling(length(sorted_values) * (upper_fraction))

  # Extract lower and upper segments
  lower_values = sorted_values[1:lower_boundary_index]
  upper_values = sorted_values[(upper_boundary_index + 1):length(sorted_values)]
  # Calculate means
  mean_low = mean(lower_values, na.rm = TRUE)
  mean_high = mean(upper_values, na.rm = TRUE)
  if (mean_low == 0) {
    return(NA)  # Avoid division by zero
  }
  if (mean_high == 0) {
    return(NA)  # Avoid division by zero
  }
  # Calculate Response Coefficient
  RC = mean_high / mean_low

  return(RC)
}
# test - passed with synthetic dataset

#calculate_RC(df_test3,env_col = 1,trait_cols = 3)
#
#mean(df_test3[1:15,3])/mean(df_test3[16:nrow(df_test3),3])


################################

#' @title Coefficient of Variation of Group Means (CVm)
#'
#' @description
#' Computes the Coefficient of Variation of Means (CVm), defined as the standard deviation
#' of group means divided by their grand mean. You can supply either:
#'
#' 1. A **list** of numeric vectors (each element is one group), or
#' 2. A **numeric vector** of trait measurements plus a matching **group_labels** vector.
#'
#' @param trait_values
#'   Either a **list** of numeric vectors (one per group), or a **numeric vector** of trait measurements.
#' @param group_labels
#'   (Optional) A vector of the same length as `trait_values` giving the group label for each observation.
#'   Required when `trait_values` is a numeric vector.
#'
#' @return
#' A single numeric value:
#' \deqn{\mathrm{CVm} = \frac{\mathrm{sd}(\bar{x}_g)}{\mathrm{mean}(\bar{x}_g)},}
#' where \(\bar{x}_g\) are the means of each group.
#'
#' @examples
#' # As a list of groups:
#' trait_groups <- list(
#'   Env1 = c(10, 12, 11),
#'   Env2 = c(20, 22, 21),
#'   Env3 = c(15, 18, 17)
#' )
#' calculate_CVm(trait_values = trait_groups)
#'
#' # As a vector with labels:
#' traits <- c(10,12,11, 20,22,21, 15,18,17)
#' labels <- rep(c("Env1","Env2","Env3"), each = 3)
#' calculate_CVm(trait_values = traits, group_labels = labels)
#'
#' @export
calculate_CVm <- function(trait_values, group_labels = NULL) {
  if (is.list(trait_values)) {
    # List-of-groups mode
    group_means <- sapply(trait_values, mean, na.rm = TRUE)
  } else {
    # Vector+labels mode
    if (is.null(group_labels)) {
      stop("`group_labels` must be provided when `trait_values` is a numeric vector.")
    }
    if (length(group_labels) != length(trait_values)) {
      stop("`group_labels` and `trait_values` must have the same length.")
    }
    group_means <- tapply(trait_values, group_labels, mean, na.rm = TRUE)
  }

  sd_means   <- sd(group_means, na.rm = TRUE)
  mean_means <- mean(group_means, na.rm = TRUE)
  if (mean_means == 0) {
    return(NA)  # Avoid division by zero
  }
  if (sd_means == 0) {
    return(0)  # If all group means are the same, CVm is 0
  }
  return(sd_means / mean_means)
}

## test - passed with synthetic dataset

#calculate_CVm(df_test3,env_col = 1,trait_col = 3)
#
#sd(c(mean(df_test3[1:10,3]),mean(df_test3[11:20,3]),mean(df_test3[21:30,3])))/mean(mean(df_test3[,3]))


###############################


#' @title Coefficient of Variation of Group Medians (CVmd)
#'
#' @description
#' Computes the Coefficient of Variation of Medians (CVmd), defined as the standard deviation
#' of group medians divided by their grand mean. You can supply either:
#'
#' 1. A **list** of numeric vectors (each element is one group), or
#' 2. A **numeric vector** of trait values plus a matching **group_labels** vector.
#'
#' @param trait_values
#'   Either a **list** of numeric vectors (one per group), or a **numeric vector** of trait measurements.
#' @param group_labels
#'   (Optional) A vector of the same length as `trait_values` giving the group label for each observation.
#'   Required when `trait_values` is a numeric vector.
#'
#' @return
#' A single numeric value:
#' \deqn{\mathrm{CVmd} = \frac{\mathrm{sd}(\tilde{x}_g)}{\mathrm{mean}(\tilde{x}_g)},}
#' where \(\tilde{x}_g\) are the medians of each group.
#'
#' @examples
#' # As a list of groups:
#' trait_groups <- list(
#'   Env1 = c(10, 12, 11),
#'   Env2 = c(20, 22, 21),
#'   Env3 = c(15, 18, 17)
#' )
#' calculate_CVmd(trait_values = trait_groups)
#'
#' # As a vector with labels:
#' traits <- c(10,12,11, 20,22,21, 15,18,17)
#' labels <- rep(c("Env1","Env2","Env3"), each = 3)
#' calculate_CVmd(trait_values = traits, group_labels = labels)
#'
#' @export
calculate_CVmd <- function(trait_values, group_labels = NULL) {
  if (is.list(trait_values)) {
    # List-of-groups mode
    medians <- sapply(trait_values, median, na.rm = TRUE)
  } else {
    # Vector+labels mode
    if (is.null(group_labels)) {
      stop("`group_labels` must be provided when `trait_values` is a numeric vector.")
    }
    if (length(group_labels) != length(trait_values)) {
      stop("`group_labels` and `trait_values` must have the same length.")
    }
    medians <- tapply(trait_values, group_labels, median, na.rm = TRUE)
  }

  sd_med   <- sd(medians, na.rm = TRUE)
  mean_med <- mean(medians, na.rm = TRUE)
  if (mean_med == 0) {
    return(NA)  #
  }
  if (sd_med == 0) {
    return(0)  # If all group medians are the same, CVmd is 0
  }
  sd_med / mean_med
}


# test - passed on synthetic dataset

#calculate_CVmd(df_test3,env_col = 1,trait_col = 2)
#
#sd(c(median(df_test3[1:10,2]),median(df_test3[11:20,2]),median(df_test3[21:30,2])))/mean(c(median(df_test3[1:10,2]),median(df_test3[11:20,2]),median(df_test3[21:30,2])))

###############################

#' @title Calculate Grand Plasticity (Pi) for a Single Trait
#'
#' @description
#' Fits a linear model of a single trait on environment (and optional covariate),
#' obtains least‐square means for each environment via **emmeans**, and then
#' computes the coefficient of variation (CV) of those means.  Optionally you
#' may exclude a control environment from the calculation.
#'
#' @param trait_values
#'   Numeric vector of trait measurements for one genotype across environments.
#' @param env_data
#'   Factor or vector (same length as `trait_values`) indicating the environment
#'   of each observation.
#' @param covariate_data
#'   Optional numeric vector (same length) to adjust for a covariate in the model.
#' @param control_env
#'   Optional level of `env_data` to treat as the control; that level's emmean will
#'   be excluded when computing the CV.
#'
#' @return
#'   A single numeric value:
#'   \eqn{\mathrm{Pi} = \frac{\mathrm{sd}(\hat\mu_e)}{\mathrm{mean}(\hat\mu_e)},}
#'   where \(\hat\mu_e\) are the emmeans for each environment (minus any control_env).
#'
#' @examples
#' trait_values   <- c(10, 12, 20, 22)
#' env_data       <- factor(c(1,1,2,2))
#' covariate_data <- c(5, 6, 7, 8)
#'
#' # Without covariate, include all envs
#' calculate_grand_plasticity(trait_values, env_data)
#'
#' # With covariate, exclude "Control"
#' calculate_grand_plasticity(
#'   trait_values   = trait_values,
#'   env_data       = env_data,
#'   covariate_data = covariate_data,
#'   control_env    = "Control"
#' )
#'
#' @importFrom stats     lm
#' @importFrom emmeans   emmeans
#' @export
calculate_grand_plasticity = function(trait_values, env_data, covariate_data = NULL, control_env = NULL) {

  env_data = factor(env_data)

  if (!is.null(covariate_data)) {
    model = lm(trait_values ~ covariate_data + env_data)
  } else {
    model = lm(trait_values ~ env_data)
  }
  adjusted_means = emmeans::emmeans(model, ~ env_data)
  adjusted_means_summary = summary(adjusted_means)

  if (!is.null(control_env)) {
    treatment_means = adjusted_means_summary$emmean[adjusted_means_summary$env_data != control_env]
  } else {
    treatment_means = adjusted_means_summary$emmean
  }

  if (length(treatment_means) == 0) {
    stop("Could not extract treatment means. Check your input data and control environment.")
  }

  sd_means = sd(treatment_means, na.rm = TRUE)
  grand_mean = mean(treatment_means, na.rm = TRUE)

  grand_plasticity = sd_means / grand_mean

  return(grand_plasticity)
}

## test - passed on synthetic dataset
#
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
#  Column1 = c(1, 1, 2, 2),      # Control (2) and Treatment (3)
#  Column0 = c(10, 12, 20, 22),  # Response variable (trait)
#  Column2 = c(1, 1, 2, 3)       # Covariate (e.g., biomass)
#)



#calculate_grand_plasticity(trait_values=df_test4[,2], env_data=df_test4[,1], covariate_data = NULL, control_env = NULL) this should produce 0.5
#calculate_grand_plasticity(df_test5,env_col = 1,trait_col = 2,covariate_col = 3,control_env = 2)

#calculate_grand_plasticity(df_test_simple,env_col = 1,trait_col = 2,covariate_col = 3,control_env = 2)

##########################################




## Example usage with synthetic data
#external_light = rep(c(0.4, 0.6, 0.8), each = 100)
#external_water = sample(rep(c("Low", "High"), each = 150))
#
#combined_factors= combine_factors(synthetic_data1, factors_not_in_dataframe = list(external_light,external_water), factors=1)

##########################################

##' @title Calculate the Phenotypic Plasticity Factor (PPF) for One or Multiple Traits
##'
##' @description Calculates the Phenotypic Plasticity Factor (PPF) for a single genotype
##' across specified environment pairs or all combinations. Models trait ~ env (+ covariates)
##' via linear regression and extracts least-square means using emmeans.
##' @param trait_values Numeric vector of trait measurements.
##' @param env_groups   Vector (factor/character) of environment labels for each measurement.
##' @param covariate_values Optional numeric vector or data.frame of covariates.
##' @param env_pairs    Optional list of length-2 character vectors. Defaults to all pairs.
##' @return Numeric scalar: mean PPF across all specified pairs (percent).
##' @examples
##' tv <- c(rnorm(5,10), rnorm(5,12)); env <- rep(c("A","B"), each=5)
##' calculate_PPF(tv, env)
calculate_PPF <- function(trait_values, env_groups,
                          covariate_values = NULL,
                          env_pairs        = NULL) {
  if (!requireNamespace("emmeans", quietly = TRUE)) {
    stop("Package 'emmeans' is required. Please install it first.")
  }

  # Prepare grouping factor
  env <- factor(env_groups)
  if (nlevels(env) < 2) {
    stop("At least two environments are required to calculate PPF.")
  }

  # Determine environment pairs
  if (is.null(env_pairs)) {
    env_pairs <- utils::combn(levels(env), 2, simplify = FALSE)
  } else if (is.vector(env_pairs) && length(env_pairs) == 2 && !is.list(env_pairs)) {
    env_pairs <- list(env_pairs)
  }

  PPF_values <- numeric(length(env_pairs))
  names(PPF_values) <- vapply(env_pairs, paste, "", collapse = "-")

  if (interactive()) {
    pb <- txtProgressBar(min = 0, max = length(env_pairs), style = 3)
    on.exit(close(pb), add = TRUE)
  }

  for (i in seq_along(env_pairs)) {
    if (interactive()) setTxtProgressBar(pb, i)
    pair <- env_pairs[[i]]

    sel    <- env %in% pair
    env_fac <- factor(env[sel], levels = pair)

    # Build data frame including covariates if present
    if (!is.null(covariate_values)) {
      if (is.data.frame(covariate_values)) {
        cov_df <- covariate_values[sel, , drop = FALSE]
      } else {
        cov_df <- data.frame(covariate = covariate_values[sel])
      }
      df <- cbind(trait = trait_values[sel], env = env_fac, cov_df)
      # Fit full model
      full_fml <- as.formula(paste("trait ~ env +", paste(names(cov_df), collapse = " + ")))
      fit_full <- tryCatch(lm(full_fml, data = df), error = function(e) e)
      # Check for collinearity (NA coefficients)
      if (inherits(fit_full, "error") || any(is.na(coef(fit_full)))) {
        warning(sprintf(
          "Collinearity detected for environments %s; using unadjusted model for PPF.",
          paste(pair, collapse = " vs ")
        ), call. = FALSE)
        # fallback to unadjusted
        df <- df[, c("trait", "env")]
        fit <- lm(trait ~ env, data = df)
      } else {
        fit <- fit_full
      }
    } else {
      df  <- data.frame(trait = trait_values[sel], env = env_fac)
      fit <- lm(trait ~ env, data = df)
    }

    # Compute LSMs
    emm <- emmeans::emmeans(fit, ~ env)
    emms <- as.data.frame(emm)
    m1 <- emms$emmean[emms$env == pair[1]]
    m2 <- emms$emmean[emms$env == pair[2]]

    # Compute PPF, guard against zero
    PPF_values[i] <- if (is.na(m1) || m1 == 0) 0 else 100 * abs(m1 - m2) / m1
  }

  return(mean(PPF_values))
}




#
#synthetic_data2=combine_factors(synthetic_data1,factors=NULL, factors_not_in_dataframe=list(external_water))
#
### test - passed on synthetic dataset
#
#df_test6 = data.frame(
#  Column1 = c(rep(3, 15), rep(2, 15)),   # Response variable
#  Column2 = c(rep(4, 15), rep(2, 15)),   # Control (2) and Treatment (3)
#  Column3 = c(rep(3, 10), rep(2, 10), rep(1, 10))    # Covariate (matches values of Column0)
#)
#df_test6$Column1 = as.factor(df_test6$Column1)
#calculate_PPF(df_test6,env_col = 1, trait_cols = 2, env_pairs = list(2,3))
#
#model = lm(Column2 ~ Column1 , data = df_test6)
#lsmeans_env = emmeans(model, ~ Column1)
#summary_lsmeans = summary(lsmeans_env)
## Extract only the LSMeans (adjusted means)
#lsmeans_values = summary_lsmeans$emmean
#100*abs((lsmeans_values[[1]]-lsmeans_values[[2]])/lsmeans_values[[1]])
#
#
#model = lm(Column3 ~ Column1 , data = df_test6)
#lsmeans_env = emmeans(model, ~ Column1)
#summary_lsmeans = summary(lsmeans_env)
## Extract only the LSMeans (adjusted means)
#lsmeans_values = summary_lsmeans$emmean
#100*abs((lsmeans_values[[1]]-lsmeans_values[[2]])/lsmeans_values[[1]])


###########################################


#' @title Phenotypic Plasticity Index (Pi)
#'
#' @description
#' Calculates the Phenotypic Plasticity Index (Pi), defined as the difference between
#' the maximum and minimum values of a trait divided by the maximum value.  Pi may be
#' computed on a single numeric vector or on specified columns of a data frame or matrix.
#'
#' @param data
#'   A numeric vector, or a data frame / matrix where each column is a trait.
#' @param trait_cols
#'   (Optional) A vector of column names or indices indicating which columns of `data`
#'   to analyze.  If `NULL` (the default), `data` must be a numeric vector and Pi is
#'   calculated for that vector.
#'
#' @return
#'   - If `trait_cols = NULL`, a single numeric value: Pi for the input vector.
#'   - Otherwise, a named numeric vector of Pi values, one per specified column.
#'
#' @examples
#' # Single-vector mode
#' x <- c(5, 8, 12, 7, 10)
#' calculate_Phenotypic_Plasticity_Index(x)
#'
#' # Data-frame mode
#' df <- data.frame(
#'   Height = c(10, 15, 12, 14, 18),
#'   Width  = c(20, 25, 22, 23, 27),
#'   Depth  = c(30, 35, 32, 34, 37)
#' )
#' # Compute Pi for Height and Width only
#' calculate_Phenotypic_Plasticity_Index(df, trait_cols = c("Height", "Width"))
#'
#' @export
calculate_Phenotypic_Plasticity_Index <- function(data, trait_cols = NULL) {
  # Single-vector mode
  if (is.null(trait_cols)) {
    if (!is.numeric(data)) {
      stop("When trait_cols is NULL, 'data' must be a numeric vector.")
    }
    max_value <- max(data, na.rm = TRUE)
    min_value <- min(data, na.rm = TRUE)
    return((max_value - min_value) / max_value)
  }

  # Multi-column mode: check input type
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("When trait_cols is provided, 'data' must be a data.frame or matrix.")
  }

  # Resolve trait_cols to names
  if (is.numeric(trait_cols)) {
    if (any(trait_cols < 1 | trait_cols > ncol(data))) {
      stop("Numeric trait_cols indices out of range.")
    }
    col_names <- colnames(data)[trait_cols]
  } else if (is.character(trait_cols)) {
    if (!all(trait_cols %in% colnames(data))) {
      missing <- setdiff(trait_cols, colnames(data))
      stop(sprintf("Columns not found: %s", paste(missing, collapse = ", ")))
    }
    col_names <- trait_cols
  } else {
    stop("trait_cols must be NULL, numeric, or character.")
  }

  # Compute Pi per column
  Pi_values <- numeric(length(col_names))
  names(Pi_values) <- col_names
  for (i in seq_along(col_names)) {
    trait_data <- data[[col_names[i]]]
    if (!is.numeric(trait_data)) {
      stop(sprintf("Column '%s' is not numeric.", col_names[i]))
    }
    max_value <- max(trait_data, na.rm = TRUE)
    min_value <- min(trait_data, na.rm = TRUE)
    Pi_values[i] <- (max_value - min_value) / max_value
  }

  # Return scalar for single trait
  if (length(Pi_values) == 1) {
    return(unname(Pi_values))
  }
  return(Pi_values)
}

##test - passed on a synthetic dataset (look at the dataset for confirmation)

#calculate_Phenotypic_Plasticity_Index(df_test6,trait_col = 2)




####################################

#' @title Proportional Inter‑Median Difference (PImd)
#'
#' @description
#' Calculates the Proportional Inter‑Median Difference (PImd), defined as the difference
#' between the maximum and minimum medians of a trait across environmental groups,
#' divided by the maximum median.  This provides a measure of relative variability
#' of the trait across the specified conditions.
#'
#' @param trait_values
#'   Numeric vector of trait measurements.
#' @param env
#'   Vector (numeric, character, or factor) indicating the environment for each measurement.
#'
#' @return
#'   A single numeric value:
#'   \deqn{\frac{\max(\tilde{x}_g) - \min(\tilde{x}_g)}{\max(\tilde{x}_g)}},
#'   where \(\tilde{x}_g\) are the medians of each environment group.
#'
#' @examples
#' trait_values <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
#' env <- factor(c('A','A','B','B','C','C','A','B','C','A'))
#' calculate_PImd(trait_values, env)
#'
#' @export
calculate_PImd <- function(trait_values, env = NULL) {
  # Input checks
  if (!is.numeric(trait_values)) {
    stop("trait_values must be a numeric vector.")
  }
  n <- length(trait_values)
  # Handle missing env
  if (is.null(env)) {
    env <- seq_len(n)
  }
  if (all(is.na(trait_values))) {
    warning("All trait_values are NA; PImd cannot be calculated.")
    return(NA_real_)
  }
  if (length(env) != n) {
    stop("Length of env (", length(env), ") must match length of trait_values (", n, ").")
  }
  # Convert to factor
  env <- factor(env)
  levels_env <- levels(env)
  # Compute medians for each level
  medians <- vapply(levels_env, function(l) {
    median(trait_values[env == l], na.rm = TRUE)
  }, numeric(1))
  # If all medians are NA, cannot compute
  if (all(is.na(medians))) {
    warning("All medians are NA; PImd cannot be calculated.")
    return(NA_real_)
  }
  max_med <- max(medians, na.rm = TRUE)
  min_med <- min(medians, na.rm = TRUE)
  # Div/0 guard
  if (max_med == 0) {
    warning("Maximum median is zero; PImd cannot be calculated.", call. = FALSE)
    return(NA_real_)
  }
  # Compute PImd
  return((max_med - min_med) / max_med)
}


##test - passed on synthetic dataset (check  dataset for confirmation)

#calculate_PImd(df_test6,trait_cols = c(2,3),env_col = 1)

###############################################



#' @title Proportional Inter‑Least‑Square Mean Difference (PILSM)
#'
#' @description
#' Computes the Proportional Inter‑Least‑Square Mean Difference (PILSM), defined as the
#' difference between the maximum and minimum least‑square means (LSMs) of a trait across
#' environments, divided by the maximum LSM.
#'
#' @param trait_values
#'   Numeric vector of trait measurements.
#' @param env
#'   Vector (numeric, character, or factor) indicating the environment for each observation.
#' @param covariates
#'   Optional numeric vector or data frame of covariate values to include in the model.
#'
#' @return
#'   A single numeric value:
#'   \deqn{\frac{\max(\hat{\mu}_e) - \min(\hat{\mu}_e)}{\max(\hat{\mu}_e)}},
#'   where \(\hat{\mu}_e\) are the estimated least‑square means for each environment.
#'
#' @examples
#' trait_values <- c(
#'   rep(10:14, each = 1),
#'   rep(20:24, each = 1),
#'   rep(30:34, each = 1)
#' )
#' env <- factor(rep(c("Env1", "Env2", "Env3"), each = 5))
#' covariates <- data.frame(SoilQuality = rep(c(3, 4, 5), each = 5))
#'
#' # Without covariates
#' calculate_PILSM(trait_values, env)
#'
#' # With covariates
#' calculate_PILSM(trait_values, env, covariates)
#'
#' @export
calculate_PILSM <- function(trait_values, env = NULL, covariates = NULL) {
  # Validate trait_values
  if (!is.numeric(trait_values)) {
    stop("trait_values must be a numeric vector.")
  }
  n <- length(trait_values)

  # Handle env
  if (is.null(env)) {
    env <- seq_len(n)
  }
  if (length(env) != n) {
    stop(sprintf("Length of env (%d) must match length of trait_values (%d).", length(env), n))
  }
  env <- factor(env)

  # Handle covariates
  if (!is.null(covariates)) {
    if (is.vector(covariates) && !is.data.frame(covariates)) {
      covariates <- data.frame(Covariate = covariates)
    }
    if (!is.data.frame(covariates)) {
      stop("covariates must be a vector or data.frame.")
    }
    if (nrow(covariates) != n) {
      stop(sprintf("Number of rows in covariates (%d) must match length of trait_values (%d).", nrow(covariates), n))
    }
  }

  # Build model data
  df <- data.frame(trait = trait_values, env = env)
  if (!is.null(covariates)) df <- cbind(df, covariates)

  # Fit model
  fml <- if (is.null(covariates)) {
    trait ~ env
  } else {
    paste("trait ~ env +", paste(names(covariates), collapse = " + "))
  }
  fit <- lm(as.formula(fml), data = df)

  # Compute least-square means
  if (!requireNamespace("emmeans", quietly = TRUE)) {
    stop("Package 'emmeans' is required for PILSM calculation.")
  }
  emm_tbl <- as.data.frame(emmeans::emmeans(fit, ~ env))
  lsms    <- emm_tbl$emmean

  # Remove NA emmeans
  if (all(is.na(lsms))) {
    warning("All LSMs are NA; PILSM cannot be calculated.")
    return(NA_real_)
  }
  max_lsm <- max(lsms, na.rm = TRUE)
  min_lsm <- min(lsms, na.rm = TRUE)

  if (max_lsm == 0) {
    warning("Maximum LSM is zero; PILSM cannot be calculated.")
    return(NA_real_)
  }

  (max_lsm - min_lsm) / max_lsm
}



## test - passed on synthetic dataset


#calculate_PILSM(df_test6,trait_col=c(2,3),env_col=1,plot = F)
#
#model = lm(Column2 ~ Column1 , data = df_test6)
#lsmeans_env = emmeans(model, ~ Column1)
#summary_lsmeans = summary(lsmeans_env)
## Extract only the LSMeans (adjusted means)
#lsmeans_values = summary_lsmeans$emmean
#(max(lsmeans_values)-min(lsmeans_values))/max(lsmeans_values)

################################################


#' @title Relative Trait Response (RTR)
#'
#' @description
#' Computes the Relative Trait Response (RTR) for a trait measured along an environmental gradient.
#' RTR is defined as the difference between the mean trait value at the high end of the gradient
#' and the mean trait value at the low end, normalized by the absolute maximum trait value:
#' \deqn{\mathrm{RTR} = \frac{\bar{z}_{\mathrm{high}} - \bar{z}_{\mathrm{low}}}{\max|z|}}.
#'
#' @param trait_values
#'   Numeric vector of trait measurements.
#' @param env_values
#'   Numeric vector of the same length as `trait_values` giving the environmental values.
#' @param env_low
#'   A fraction (0 < frac < 1) specifying the low‐end cutoff.
#'   Defaults to 0.2 (the 20th percentile).
#' @param env_high
#'   A fraction (0 < frac < 1) specifying the high‐end cutoff.
#'   Defaults to 0.2 (the 80th percentile).
#' @return
#'   A single numeric value: the Relative Trait Response.
#'
#' @examples
#' trait_values <- c(
#'   10,12,11,13,14,15,13,14,12,13,
#'   20,22,21,23,24,25,23,24,22,23
#' )
#' env_values <- seq(1, 10, length.out = 20)
#'
#' # Using 20% tails
#' calculate_RTR(trait_values, env_values)
#'
#' # Using explicit thresholds
#' calculate_RTR(trait_values, env_values, env_low = 2, env_high = 9)
#'
#' @export
calculate_RTR <- function(trait_values, env_values,
                          env_low = 0.2, env_high = 0.2) {


  if (!is.numeric(env_low) || length(env_low) != 1 || env_low < 0 || env_low > 1) {
    warning("env_low must be a single numeric value between 0 and 1.")
    return(NA)
  }
  if (!is.numeric(env_high) || length(env_high) != 1 || env_high < 0 || env_high > 1) {
    warning("env_high must be a single numeric value between 0 and 1.")
    return(NA)
  }

  if (env_low > env_high) {
    warning("env_low cannot be higher than env_high.")
    return(NA)
  }
  if(all(is.na(trait_values))){
    warning("All observations have NA in trait_values or env_values; RTR cannot be calculated.")
    return(NA_real_)
  }
  if (!is.numeric(trait_values) || !is.numeric(env_values)) {
    stop("trait_values and env_values must be numeric vectors.")

  }
  n <- length(trait_values)
  if (length(env_values) != n) {
    stop("trait_values and env_values must have the same length.")
  }
  # Handle missing data
  keep <- !(is.na(trait_values) | is.na(env_values))
  if (!any(keep)) {
    warning("All observations have NA in trait_values or env_values; RTR cannot be calculated.")
    return(NA_real_)
  }
  trait <- trait_values[keep]
  env    <- env_values[keep]



  if (env_low >= 0 && env_low <= 1) {
    low_thresh <- stats::quantile(env, probs = env_low, na.rm = TRUE)
  } else {
    low_thresh <- env_low
  }
  # Determine high threshold
  if (!is.numeric(env_high) || length(env_high) != 1) {
    stop("env_high must be a single numeric value.")
  }
  if (env_high >= 0 && env_high <= 1) {
    high_thresh <- stats::quantile(env, probs = 1 - env_high, na.rm = TRUE)
  } else {
    high_thresh <- env_high
  }

  # Subset observations
  low_idx  <- which(env <= low_thresh)
  high_idx <- which(env >= high_thresh)
  if (length(low_idx) == 0 || length(high_idx) == 0) {
    warning("No observations found in low or high env thresholds; RTR cannot be calculated.")
    return(NA_real_)
  }

  # Compute means
  mean_low  <- mean(trait[low_idx],  na.rm = TRUE)
  mean_high <- mean(trait[high_idx], na.rm = TRUE)

  # Normalize by max absolute trait
  denom <- max(abs(trait), na.rm = TRUE)
  if (denom == 0) {
    warning("Maximum absolute trait value is zero; RTR cannot be calculated.")
    return(NA_real_)
  }

  # Calculate RTR
  (mean_high - mean_low) / denom
}




## test - passed on synthetic dataset
#calculate_RTR(df_test2,trait_col=2,env_col=1,env_low=1,env_high=2)


######################################

#' @title Phenotypic Instability Ratio (PIR)
#'
#' @description
#' Computes the Phenotypic Instability Ratio (PIR) for a trait measured across environments,
#' following Robinson (1989). PIR is defined as the difference between the maximum and minimum
#' mean trait values across environments, divided by the mean trait value in the environment
#' with the highest relative growth rate (RGR).
#'
#' @param trait_values
#'   Numeric vector of trait measurements.
#' @param env_values
#'   Optional vector (numeric, character, or factor) of the same length as `trait_values`
#'   indicating the environment for each measurement. If `NULL` (default), environments
#'   are treated as equally spaced (1, 2, …).
#' @param rgr_values
#'   Optional numeric vector of relative growth rates corresponding to `trait_values`.
#'   If `NULL` (default), RGR is estimated from successive differences in environment means.
#'
#' @details
#' 1. Convert or generate `env_values` into a factor.
#' 2. Compute the mean trait value in each environment.
#' 3. Determine `MaxMean` and `MinMean` from these environment means.
#' 4. Obtain or estimate RGR per environment; identify the environment with maximum RGR.
#' 5. Compute \eqn{\mathrm{PIR} = (MaxMean - MinMean) / Mean_{MaxRGR}}.
#'
#' @return
#'   Single numeric PIR value.
#'
#' @references
#' Robinson, D. (1989). *Plasticity in plant growth and resource use as a trait for crop breeding*.
#' Field Crops Research, 11(2-3), 153–159.
#'
#' @examples
#' trait_values <- c(
#'   10,12,11,13,14,15,13,14,12,13,
#'   20,22,21,23,24,25,23,24,22,23,
#'   30,32,31,33,34,35,33,34,32,33
#' )
#' env_values <- rep(c("Env1","Env2","Env3"), each = 10)
#' calculate_PIR(trait_values, env_values)
#'
#' # Equidistant environments
#' calculate_PIR(trait_values)
#'
#' @export
calculate_PIR <- function(trait_values, env_values = NULL, rgr_values = NULL) {
  # Validate trait_values
  if (!is.numeric(trait_values)) {
    stop("trait_values must be a numeric vector.")
  }
  n <- length(trait_values)

  # Handle env_values
  if (is.null(env_values)) {
    env <- factor(seq_len(n))
  } else {
    if (length(env_values) != n) {
      stop("env_values must have the same length as trait_values.")
    }
    env <- factor(env_values)
  }

  # Remove missing observations
  keep <- !is.na(trait_values) & !is.na(env)
  if (!any(keep)) {
    warning("All observations are missing; PIR cannot be calculated.")
    return(NA_real_)
  }
  trait <- trait_values[keep]
  env   <- droplevels(env[keep])

  # Compute environment means
  means <- tapply(trait, env, mean, na.rm = TRUE)
  if (length(means) < 2) {
    stop("At least two environments are required to calculate PIR.")
  }
  max_mean <- max(means, na.rm = TRUE)
  min_mean <- min(means, na.rm = TRUE)

  # Compute or aggregate RGR values
  if (is.null(rgr_values)) {
    mvec <- unname(means)
    rgr <- numeric(length(mvec))
    rgr[1] <- NA_real_
    denom <- mvec[-length(mvec)]
    rgr[-1] <- ifelse(denom == 0, NA_real_, diff(mvec) / denom)
  } else {
    if (!is.numeric(rgr_values) || length(rgr_values) != n) {
      stop("rgr_values must be a numeric vector the same length as trait_values.")
    }
    rgr_keep <- rgr_values[keep]
    rgr      <- tapply(rgr_keep, env, mean, na.rm = TRUE)
    if (length(rgr) != length(means)) {
      stop("Mismatch between RGR values and environment means.")
    }
  }

  # Select environment with maximum RGR
  if (all(is.na(rgr))) {
    warning("All RGR values are NA; PIR cannot be calculated.")
    return(NA_real_)
  }
  idx <- which.max(rgr)
  mean_at_max_rgr <- means[idx]
  if (is.na(mean_at_max_rgr) || mean_at_max_rgr == 0) {
    warning("Mean at maximum RGR is zero or NA; PIR cannot be calculated.")
    return(NA_real_)
  }

  # Compute PIR
  PIR <- (max_mean - min_mean) / mean_at_max_rgr
  PIR
}
