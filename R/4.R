#Calculate Developmental Plasticity Index (DPI)(calculate_DPI) - tested,
#Calculate Coefficient of Environmental Variation (CEV)(calculate_CEV) - tested,
#Calculate Plasticity Response Index (PRI)(calculate_PRI) - tested,
#Calculate Phenotypic Flexibility Index (PFI)(calculate_PFI) - tested,
#Calculate Standardized Plasticity Index (SPI)(calculate_SPI) - tested,
#Calculate Absolute Plasticity Coefficient (APC)(calculate_APC) - tested,
#Calculate Stability Index (SI)(calculate_SI) - tested,
#Calculate Relative Stability Index (RSI)(calculate_RSI),
#Calculate Environmental Variance Sensitivity (EVS)(calculate_EVS) - tested,
#Calculate Environmental Variance Sensitivity (EVS)(calculate_EVS) - tested,
#Calculate Multivariate Plasticity Index (MVPi)(calculate_MVPi) NEEDS TO BE TESTED WITH THE REQUESTED DATASET FROM PROF. BARBOSA,
#Calculate Standardized Plasticity Metric (SPM)(calculate_SPM) - tested,
#Calculate SSpop/SStotal Plasticity Ratio(calculate_Plasticity_Ratio) - tested



##################### datasets for testing

#df_test1 = data.frame(Column1 = c(rep(4, 10), rep(2, 10)),
#                      Column2 = c(rep(10, 10), rep(1, 10)))
#
#df_test2 = data.frame(Column0 = c(rep(1, 10), rep(2, 10)),
#                      Column1 = c(rep(2, 10), rep(1, 10)),
#                      Column2 = c(rep(2, 10), rep(4, 10)))
#
#df_test2.2 = data.frame(Column0 = c(rep(1, 10), rep(2, 10)),
#                        Column1 = c(rep(2, 10), rep(1, 5),rep(2,5)),
#                        Column2 = c(rep(2, 10), rep(4, 10)))
#
#
#
#df_test3=data.frame(Column0 = c(rep(3, 10), rep(2, 10),rep(1,10)),
#                    Column1 = c(rep(2, 10), rep(1, 15),rep(3,5)),
#                    Column2 = c(rep(2, 10), rep(4, 10),rep(3,10)))
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
#  Column1 = c(2, 2, 3, 3),      # Control (2) and Treatment (3)
#  Column0 = c(10, 12, 20, 22),  # Response variable (trait)
#  Column2 = c(1, 1, 2, 2)       # Covariate (e.g., biomass)
#)



#########################################

#' @title Calculate Developmental Plasticity Index (DPI)
#' @description This function calculates the Developmental Plasticity Index (DPI) by quantifying how a trait value
#' changes between two time points during development.
#'
#' The formula used is:
#' \deqn{DPI = \frac{Trait Value at Time2 - Trait Value at Time1}{Time Interval}}
#'
#' @param trait_values_time1 A numeric vector containing trait values at Time 1.
#' @param trait_values_time2 A numeric vector containing trait values at Time 2.
#' @param time_interval Optional numeric vector of the time intervals between the two time points for each sample.
#' If not provided, a vector of ones is assumed.
#'
#' @return A numeric vector of DPI values.
#'
#' @examples
#' # Example usage with synthetic data and explicit time intervals
#' trait_time1 = c(10, 12, 14)
#' trait_time2 = c(15, 18, 22)
#' time_intervals = c(5, 5, 5)
#' dpi_values = calculate_DPI(trait_time1, trait_time2, time_intervals)
#' print(dpi_values)
#'
#' # Example usage with default time interval (all ones)
#' dpi_values_default = calculate_DPI(trait_time1, trait_time2)
#' print(dpi_values_default)
#'
#' @export
calculate_DPI = function(trait_values_time1, trait_values_time2, time_interval = NULL) {
  # Ensure trait values are numeric vectors
  if (!is.numeric(trait_values_time1) || !is.numeric(trait_values_time2)) {
    stop("trait_values_time1 and trait_values_time2 must be numeric vectors.")
  }

  # If no time_interval vector is provided, assume a vector of ones
  if (is.null(time_interval)) {
    time_interval = rep(1, length(trait_values_time1))
  } else if (!is.numeric(time_interval)) {
    stop("time_interval must be a numeric vector.")
  }

  # Ensure all vectors have the same length
  if (length(trait_values_time1) != length(trait_values_time2) || length(trait_values_time1) != length(time_interval)) {
    stop("trait_values_time1, trait_values_time2, and time_interval must all have the same length.")
  }

  # Ensure all elements of time_interval are positive numeric values
  if (any(time_interval <= 0)) {
    stop("All elements of time_interval must be positive numeric values.")
  }

  # Calculate DPI for each sample
  dpi = (trait_values_time2 - trait_values_time1) / time_interval

  return(dpi)
}




## test - passed on synthetic dataset - however it is not clear how one would structure a dataset with time resolved data (will the measurements at time x for one sample be in different columns all with equally separated time intervals or in rows?)
# mathematically the function does the wanted if one knows how to structure the data

#calculate_DPI(df_test2,2,3,3)



###############################



#' @title Calculate Coefficient of Environmental Variation (CEV)
#'
#' @description
#' Computes the Coefficient of Environmental Variation (CEV) for a single trait,
#' defined as the standard deviation of the trait measurements across environments
#' divided by their mean, then multiplied by 100 to express as a percentage:
#' \deqn{\mathrm{CEV} = \frac{\mathrm{sd}(x)}{\mathrm{mean}(x)} \times 100.}
#' If fewer than two non‐missing values are provided, or if the mean is zero,
#' the function returns \code{NA}.
#'
#' @param trait_values
#'   A numeric vector of trait values measured across environments.
#'
#' @return
#' A single numeric value giving the CEV (in percent). Returns \code{NA} if:
#' \itemize{
#'   \item fewer than two non‐missing values are supplied,
#'   \item or the mean of \code{trait_values} is zero.
#' }
#'
#' @examples
#' # A simple example with equidistant environments
#' trait_values <- c(10, 15, 20, 25, 30)
#' calculate_CEV(trait_values)
#'
#' # With missing values
#' trait_values2 <- c(5, NA, 15, 20)
#' calculate_CEV(trait_values2)
#'
#' # Too few values returns NA
#' calculate_CEV(c(42))
#'
#' @export
calculate_CEV = function(trait_values) {
  # Ensure input is a numeric vector
  if (!is.numeric(trait_values)) {
    stop("trait_values must be a numeric vector.")
  }

  # Check that there are enough data points to compute a standard deviation
  if (length(trait_values) < 2) {
    return(list(CEV = NA, Mean = NA, SD = NA, Valid = FALSE))
  }

  # Calculate mean and standard deviation (ignoring NA values)
  mean_val = mean(trait_values, na.rm = TRUE)
  sd_val = sd(trait_values, na.rm = TRUE)

  # Avoid division by zero if the mean is zero
  if (mean_val == 0) {
    return(list(CEV = NA, Mean = mean_val, SD = sd_val, Valid = FALSE))
  }

  # Calculate the Coefficient of Environmental Variation (CEV)
  cev = (sd_val / mean_val) * 100

  return(cev)
}



## tested - passed on synthetic dataset

#calculate_CEV(df_test2,c(2,3))
#sd(df_test2[,3])/mean(df_test2[,3])*100

################################


#' @title Calculate Plasticity Response Index (PRI)
#'
#' @description
#' Computes the Plasticity Response Index (PRI) for a single trait, defined as the
#' difference between the mean trait value under extreme conditions and the mean under
#' control conditions, normalized by the overall mean:
#' \deqn{PRI = \frac{\overline{x}_{\rm extreme} - \overline{x}_{\rm control}}{\overline{x}_{\rm overall}}.}
#'
#' @param trait_values
#'   Numeric vector of trait measurements.
#' @param env_indicator
#'   Binary vector (same length as \code{trait_values}) indicating environment:
#'   \code{0} = control, \code{1} = extreme.
#'
#' @return
#' A single numeric value: the calculated PRI.
#' Returns \code{NA} (with a warning) if the overall mean of \code{trait_values} is zero.
#'
#' @examples
#' # Example: half control, half extreme
#' trait_values  <- c(10, 15, 20, 25, 30)
#' env_indicator  <- c(0, 0, 1, 1, 1)
#' calculate_PRI(trait_values, env_indicator)
#'
#' # Edge case: overall mean zero
#' tv <- c(-1, 1, -1, 1)
#' ei <- c(0, 0, 1, 1)
#' calculate_PRI(tv, ei)
#'
#' @export
calculate_PRI <- function(trait_values, env_indicator) {
  # Validate input
  if (!is.numeric(trait_values)) {
    stop("trait_values must be a numeric vector.")
  }
  if (length(trait_values) != length(env_indicator)) {
    stop("trait_values and env_indicator must have the same length.")
  }
  if (!all(env_indicator %in% c(0, 1))) {
    stop("env_indicator must be a binary vector with values 0 (control) and 1 (extreme).")
  }

  # Calculate overall mean
  overall_mean <- mean(trait_values, na.rm = TRUE)
  if (overall_mean == 0) {
    warning("Overall mean of trait values is 0; PRI is undefined. Returning NA.")
    return(NA)
  }

  # Calculate means for extreme and control environments
  mean_extreme <- mean(trait_values[env_indicator == 1], na.rm = TRUE)
  mean_control <- mean(trait_values[env_indicator == 0], na.rm = TRUE)

  # Calculate PRI
  pri_value <- (mean_extreme - mean_control) / overall_mean

  return(pri_value)
}


## test - passed on synthetic dataset

#calculate_PRI(df_test2,1,2,1,2)

############################

#' @title Calculate Phenotypic Flexibility Index (PFI) for a Single Trait
#'
#' @description
#' Computes the Phenotypic Flexibility Index (PFI) for a single trait. PFI is defined as the
#' maximum absolute deviation of any sample from its baseline value, divided by that baseline.
#' If no `baseline_values` are provided, the overall mean of `trait_values` is used as the baseline.
#'
#' @param trait_values
#'   Numeric vector of trait measurements.
#' @param baseline_values
#'   Optional numeric vector of baseline values (same length as `trait_values`).
#'   If \code{NULL}, uses \code{mean(trait_values)} for all samples.
#'
#' @return
#' A single numeric value: the PFI, i.e.
#' \deqn{\max_i \frac{|x_i - b_i|}{b_i},}
#' where \(x_i\) are the trait values and \(b_i\) the corresponding baselines.
#' Returns \code{NA} (with a warning) if the baseline for the sample with maximum deviation is zero.
#'
#' @examples
#' # Provided baseline:
#' trait_values    <- c(18, 22, 25, 20, 19)
#' baseline_values <- c(20, 20, 20, 20, 20)
#' calculate_PFI(trait_values, baseline_values)
#'
#' # No baseline provided (uses overall mean):
#' trait_values <- c(18, 22, 25, 20, 19)
#' calculate_PFI(trait_values)
#'
#' @export
calculate_PFI <- function(trait_values, baseline_values = NULL) {
  # Validate input for trait_values
  if (!is.numeric(trait_values)) {
    stop("trait_values must be a numeric vector.")
  }

  n <- length(trait_values)

  # If baseline_values is not provided, use the overall mean for all samples.
  if (is.null(baseline_values)) {
    baseline_values <- rep(mean(trait_values, na.rm = TRUE), n)
  } else {
    # Validate baseline_values
    if (!is.numeric(baseline_values)) {
      stop("baseline_values must be a numeric vector.")
    }
    if (length(baseline_values) != n) {
      stop("baseline_values must have the same length as trait_values.")
    }
  }

  # Calculate the absolute deviations from the baseline for each sample
  deviations <- abs(trait_values - baseline_values)

  # Identify the sample with the maximum deviation
  max_deviation_index <- which.max(deviations)
  max_deviation <- deviations[max_deviation_index]

  # Retrieve the baseline corresponding to the maximum deviation sample
  baseline_for_max <- baseline_values[max_deviation_index]

  if (baseline_for_max == 0) {
    warning("The baseline value corresponding to the maximum deviation is 0; PFI is undefined. Returning NA.")
    return(NA)
  }

  # Calculate the PFI as the ratio of the maximum deviation to its corresponding baseline
  pfi_value <- max_deviation / baseline_for_max

  return(pfi_value)
}

## test - passed on synthetic dataset - however I am unsure wether this is supposed to be one score or a separate score for each timepoint/sample

#calculate_PFI(df_test2,baseline_cols=c(1,1),trait_cols=c(2,3))


##########################

#' @title Calculate Standardized Plasticity Index (SPI)
#'
#' @description
#' Computes the Standardized Plasticity Index (SPI) for one or more traits by comparing trait
#' values between two environments and scaling the difference by the variability in a reference environment.
#'
#' The index is defined as:
#' \deqn{SPI = \frac{\bar{x}_{env2} - \bar{x}_{env1}}{\mathrm{sd}(\,x_{ref}\,)},}
#' where \(\bar{x}_{env1}\) and \(\bar{x}_{env2}\) are the mean trait values in the two focal environments,
#' and \(\mathrm{sd}(x_{ref})\) is the standard deviation of trait values in the reference environment.
#'
#' @param data
#'   A data frame containing the observations.
#' @param env_col
#'   Either a column name or index in \code{data} that holds the environment labels, or an external
#'   vector of length \code{nrow(data)} specifying the environment for each row.
#' @param trait_cols
#'   A character or numeric vector indicating which columns in \code{data} contain the trait measurements.
#' @param env1
#'   The label (value) of the first environment in \code{env_col} to compare.
#' @param env2
#'   The label (value) of the second environment in \code{env_col} to compare.
#' @param reference_env
#'   The label (value) of the reference environment in \code{env_col} whose standard deviation is used for scaling.
#'
#' @return
#' A named numeric vector of SPI values, one entry for each trait in \code{trait_cols}.
#'
#' @examples
#' df <- data.frame(
#'   Height = c(10, 12, 15, 18),
#'   Weight = c(5,  6,  7,  8),
#'   Env    = factor(c("A", "A", "B", "B"))
#' )
#' # Compute SPI comparing B vs. A, with A as the reference
#' calculate_SPI(
#'   data          = df,
#'   env_col       = "Env",
#'   trait_cols    = c("Height", "Weight"),
#'   env1          = "A",
#'   env2          = "B",
#'   reference_env = "A"
#' )
#'
#' @export
calculate_SPI = function(data, env_col, trait_cols, env1, env2, reference_env) {

  # Handle env_col
  if (is.numeric(env_col) && length(env_col) == 1) {
    env_col = data[[env_col]]
  } else if (is.vector(env_col) && length(env_col) == nrow(data)) {
    env_col = env_col
  } else {
    env_col = data[[env_col]]
  }

  # Initialize a list to store SPI values for each trait
  spi_results = numeric(length(trait_cols))
  names(spi_results) = trait_cols

  i=0
  # Loop through each trait column and calculate SPI
  for (trait_col in trait_cols) {
    i=i+1
    # Extract trait data
    trait_values = if (is.numeric(trait_col)) data[[trait_col]] else data[[trait_col]]

    # Calculate mean trait values for env1 and env2
    trait_env1_value = mean(trait_values[env_col == env1], na.rm = TRUE)
    trait_env2_value = mean(trait_values[env_col == env2], na.rm = TRUE)

    # Calculate the standard deviation of the reference environment's trait values
    sd_reference = sd(trait_values[env_col == reference_env], na.rm = TRUE)

    # Calculate the Standardized Plasticity Index (SPI) for this trait
    spi = (trait_env2_value - trait_env1_value) / sd_reference

    # Store the result
    spi_results[i] = spi
  }

  return(spi_results)
}

## test - passed on synthetic dataset

#calculate_SPI(df_test2.2,1,2,1,2,2)
#-0.5/sd(c(rep(1, 5),rep(2,5)))




########################

#' @title Absolute Plasticity Coefficient (APC) for a Single Trait
#'
#' @description
#' Computes the Absolute Plasticity Coefficient (APC) by averaging the absolute differences
#' between mean trait values in consecutive environments. If no environment labels are provided,
#' environments are assumed equidistant (1, 2, …).
#'
#' The index is defined as:
#' \deqn{APC = \frac{1}{n - 1} \sum_{i=1}^{n-1} \left| \mu_{i+1} - \mu_i \right|,}
#' where \(\mu_i\) is the mean trait value in environment \(i\) and \(n\) is the number of unique environments.
#'
#' @param trait_values
#'   A numeric vector of trait measurements.
#' @param env_labels
#'   (Optional) A vector of environment identifiers (numeric, character, or factor) of the same length
#'   as `trait_values`. If `NULL` (default), environments are taken to be `1:length(trait_values)`.
#' @param sequential_env
#'   Logical; if `TRUE` (default), environments are ordered by their natural/ascending order before
#'   computing means. If `FALSE`, environments are used in the order of first appearance.
#'
#' @return
#' A single numeric value: the APC, or `NA` (with a warning) if fewer than two unique environments are found.
#'
#' @examples
#' # With explicit environment labels:
#' trait_values <- c(10, 12, 15, 20, 22, 25)
#' env_labels   <- c(1, 1, 2, 2, 3, 3)
#' calculate_APC(trait_values, env_labels, sequential_env = TRUE)
#'
#' # Without env_labels (assumes equidistant environments):
#' calculate_APC(trait_values)
#'
#' @export
calculate_APC <- function(trait_values, env_labels = NULL, sequential_env = TRUE) {
  # Input validation for trait_values
  if (!is.numeric(trait_values)) {
    stop("trait_values must be a numeric vector.")
  }

  # If no environment labels are provided, assume equidistant environments.
  if (is.null(env_labels)) {
    env_labels <- seq_along(trait_values)
  } else {
    if (length(trait_values) != length(env_labels)) {
      stop("trait_values and env_labels must have the same length.")
    }
  }

  # Order the data based on the environment labels.
  if (sequential_env) {
    if (is.numeric(env_labels)) {
      ordering <- order(env_labels)
    } else if (is.factor(env_labels)) {
      if (is.ordered(env_labels)) {
        ordering <- order(as.numeric(env_labels))
      } else {
        ordering <- seq_along(env_labels)
      }
    } else {  # For character vectors
      ordering <- order(env_labels)
    }
  } else {
    # Use the order of first appearance
    uniq_env <- unique(env_labels)
    ordering <- order(match(env_labels, uniq_env))
  }

  trait_values <- trait_values[ordering]
  env_labels <- env_labels[ordering]

  # Calculate mean trait value for each unique environment
  env_means <- tapply(trait_values, env_labels, mean, na.rm = TRUE)
  if (length(env_means) < 2) {
    warning("Less than two unique environments found; APC is undefined. Returning NA.")
    return(NA)
  }

  # Compute absolute differences between consecutive environment means
  abs_diff <- abs(diff(env_means))

  # Calculate APC as the mean of these absolute differences
  apc_value <- mean(abs_diff, na.rm = TRUE)

  return(apc_value)
}

## test - passed on synthetic dataset
#calculate_APC(df_test2,1,2)

#############################


#' @title Stability Index (SI) for a Single Trait
#'
#' @description
#' Computes the Stability Index (SI), quantifying how stable a trait is across environments.
#' SI is defined as the ratio of the variance among environment means to their overall mean:
#' \deqn{SI = \frac{\mathrm{Var}(\bar{x}_e)}{\mathrm{Mean}(\bar{x}_e)},}
#' where \(\bar{x}_e\) are the mean trait values in each environment.
#'
#' @param trait_values
#'   A numeric vector of trait measurements.
#' @param env
#'   (Optional) A vector (numeric, character, or factor) of the same length as `trait_values`
#'   indicating the environment for each measurement. If `NULL`, environments are assumed
#'   to be equidistant (`1:length(trait_values)`).
#'
#' @return
#' A single numeric value: the Stability Index (SI). Returns `NA` (with a warning)
#' if the overall mean of environment means is zero.
#'
#' @examples
#' # Explicit environments:
#' trait_values <- c(10, 12, 15, 20)
#' env          <- c(1, 2, 3, 4)
#' calculate_SI(trait_values, env)
#'
#' # Equidistant environments assumed:
#' trait_values <- c(10, 12, 15, 20)
#' calculate_SI(trait_values)
#'
#' @export
calculate_SI <- function(trait_values, env = NULL) {
  # Validate trait_values
  if (!is.numeric(trait_values)) {
    stop("trait_values must be a numeric vector.")
  }

  # If no environment vector is provided, assume equidistant environments.
  if (is.null(env)) {
    env <- seq_along(trait_values)
  } else {
    if (length(trait_values) != length(env)) {
      stop("trait_values and env must have the same length.")
    }
  }

  # Compute the mean trait value for each unique environment.
  env_means <- tapply(trait_values, env, mean, na.rm = TRUE)

  # Calculate the variance among environment means and their overall mean.
  variance_env <- var(env_means, na.rm = TRUE)
  overall_mean <- mean(env_means, na.rm = TRUE)

  if (overall_mean == 0) {
    warning("The overall mean of environment means is 0; SI is undefined. Returning NA.")
    return(NA)
  }

  # Compute the Stability Index (SI)
  SI_value <- variance_env / overall_mean

  return(SI_value)
}


### test - passed on synthetic dataset
#var(df_test2[,2])/mean(df_test2[,2])
#calculate_SI(df_test2,1,2)
#######################


#' @title Relative Stability Index (RSI) for a Single Trait
#'
#' @description
#' Computes the Relative Stability Index (RSI), quantifying the relative stability of a trait
#' across environments. RSI is defined as:
#' \deqn{RSI = 1 - \frac{\mathrm{SD}(\bar{x}_e)}{\mathrm{Mean}(\bar{x}_e)},}
#' where \(\bar{x}_e\) are the mean trait values in each environment.
#'
#' @param trait_values
#'   A numeric vector of trait measurements.
#' @param env
#'   (Optional) A vector (numeric, character, or factor) of the same length as `trait_values`
#'   indicating the environment for each measurement. If `NULL`, environments are assumed
#'   to be equidistant (`1:length(trait_values)`).
#'
#' @return
#' A single numeric value: the Relative Stability Index (RSI). Returns `NA` (with a warning)
#' if the overall mean of environment means is zero.
#'
#' @examples
#' # Explicit environments:
#' trait_values <- c(10, 12, 15, 20)
#' env          <- c(1, 2, 3, 4)
#' calculate_RSI(trait_values, env)
#'
#' # Equidistant environments assumed:
#' trait_values <- c(10, 12, 15, 20)
#' calculate_RSI(trait_values)
#'
#' @export
calculate_RSI <- function(trait_values, env = NULL) {
  # Validate trait_values
  if (!is.numeric(trait_values)) {
    stop("trait_values must be a numeric vector.")
  }

  # If no environment vector is provided, assume equidistant environments.
  if (is.null(env)) {
    env <- seq_along(trait_values)
  } else {
    if (length(trait_values) != length(env)) {
      stop("trait_values and env must have the same length.")
    }
  }

  # Compute the mean trait value for each unique environment.
  env_means <- tapply(trait_values, env, mean, na.rm = TRUE)

  # Calculate the standard deviation and mean of the environment means.
  overall_sd <- sd(env_means, na.rm = TRUE)
  overall_mean <- mean(env_means, na.rm = TRUE)

  if (overall_mean == 0) {
    warning("The overall mean of environment means is 0; RSI is undefined. Returning NA.")
    return(NA)
  }

  # Compute the Relative Stability Index (RSI)
  RSI_value <- 1 - (overall_sd / overall_mean)

  return(RSI_value)
}


### test - passed on synthetic dataset
#1-sd(df_test2[,2])/mean(df_test2[,2])
#calculate_RSI(df_test2,1,2)
#
###test2
#1-sd(synthetic_data1[,3])/mean(synthetic_data1[,3])
#calculate_RSI(synthetic_data1,1,3)
#
###test3
#1-sd(df_test5[1:20,2])/mean(df_test5[1:20,2])
#calculate_RSI(df_test5,1,2,specific_envs=c(3,2))



##########################

#' @title Calculate Environmental Variance Sensitivity (EVS) for a Single Trait
#'
#' @description
#' Computes the Environmental Variance Sensitivity (EVS), which quantifies how sensitive
#' a trait is to changes in the environment. EVS is defined as the ratio of the variance
#' in the trait values to the variance in the environmental values:
#' \deqn{EVS = \frac{\mathrm{Var}(\text{trait values})}{\mathrm{Var}(\text{environment})}.}
#'
#' @param trait_values
#'   A numeric vector of trait measurements.
#' @param env
#'   (Optional) A vector (numeric, character, or factor) of the same length as `trait_values`
#'   indicating the environment for each measurement. If `NULL`, environments are assumed
#'   to be equidistant (`1:length(trait_values)`).
#'
#' @return
#' A single numeric value: the Environmental Variance Sensitivity (EVS). Returns `NA`
#' (with a warning) if the variance of the environmental values is zero.
#'
#' @examples
#' # Explicit environments:
#' trait_values <- c(10, 12, 15, 18)
#' env          <- c(1, 2, 3, 4)
#' calculate_EVS(trait_values, env)
#'
#' # Equidistant environments assumed:
#' trait_values <- c(10, 12, 15, 18)
#' calculate_EVS(trait_values)
#'
#' @export
calculate_EVS <- function(trait_values, env = NULL) {
  # Validate trait_values
  if (!is.numeric(trait_values)) {
    stop("trait_values must be a numeric vector.")
  }

  # If no environment vector is provided, assume equidistant environments.
  if (is.null(env)) {
    env <- seq_along(trait_values)
  } else {
    if (length(trait_values) != length(env)) {
      stop("trait_values and env must have the same length.")
    }
  }

  # Compute the variance in trait values.
  trait_variance <- var(trait_values, na.rm = TRUE)

  # If env is not numeric, convert it to numeric via as.factor.
  if (!is.numeric(env)) {
    env_numeric <- as.numeric(as.factor(env))
  } else {
    env_numeric <- env
  }

  # Compute the variance in the environment.
  env_variance <- var(env_numeric, na.rm = TRUE)

  if (env_variance == 0) {
    warning("Environmental variance is zero; EVS is undefined. Returning NA.")
    return(NA)
  }

  # Calculate EVS as the ratio of the two variances.
  EVS_value <- trait_variance / env_variance

  return(EVS_value)
}

## test - passed on synthetic dataset

#var(df_test2[,2])/var(df_test2[,1])
#calculate_EVS(df_test2,1,2)
#
###test2 - passed
#var(synthetic_data1[,3])/var(synthetic_data1[,2])
#calculate_EVS(synthetic_data1,2,3)

#####################

#' @title Calculate Multivariate Plasticity Index (MVPi)
#'
#' @description
#' Computes the Multivariate Plasticity Index (MVPi) for a set of phenotypic traits by
#' performing PCA and measuring the average Euclidean distance between environment centroids
#' in the reduced-dimensional space. This quantifies overall plasticity across environments.
#'
#' @param trait_data
#'   A data frame, matrix, or numeric vector of phenotypic trait values (rows = samples).
#' @param env
#'   (Optional) A vector (numeric, character, or factor) of length equal to the number of rows
#'   in `trait_data` indicating each sample's environment. If `NULL`, equidistant environments
#'   (`1:nrow(trait_data)`) are assumed.
#' @param n_axes
#'   Integer ≥ 1 giving the number of principal component axes to retain from the PCA.
#'
#' @return
#' A single numeric value: the mean Euclidean distance between all pairs of environment centroids
#' in the first `n_axes` PCA dimensions. Returns `NA` (with a warning) if fewer than two environments.
#'
#' @examples
#' # Example with explicit environments:
#' trait_data <- data.frame(
#'   Trait1 = c(1.2, 2.3, 1.9, 3.1, 2.0, 2.8),
#'   Trait2 = c(2.5, 3.0, 2.8, 3.6, 2.9, 3.3)
#' )
#' env <- factor(c(1, 1, 2, 2, 3, 3))
#' mvpi_value <- calculate_MVPi(trait_data, env, n_axes = 2)
#' print(mvpi_value)
#'
#' # Example with no env vector (equidistant environments assumed):
#' mvpi_value2 <- calculate_MVPi(trait_data, n_axes = 2)
#' print(mvpi_value2)
#'
#' @export
calculate_MVPi <- function(trait_data, env = NULL, n_axes = 1) {
  # Convert trait_data to a numeric matrix if not already one.
  if (!is.matrix(trait_data)) {
    trait_data <- as.matrix(trait_data)
  }
  if (!is.numeric(trait_data)) {
    stop("trait_data must be numeric.")
  }

  # Determine number of samples.
  n_samples <- nrow(trait_data)

  # If no environment vector is provided, assume equidistant environments.
  if (is.null(env)) {
    env <- seq_len(n_samples)
  } else {
    if (length(env) != n_samples) {
      stop("The length of env must equal the number of rows in trait_data.")
    }
  }

  # Convert env to a factor.
  env <- factor(env)

  # Perform PCA on the trait data.
  pca_result <- prcomp(trait_data, center = TRUE, scale. = FALSE)

  # Ensure that n_axes does not exceed the available dimensions.
  if (n_axes > ncol(pca_result$x)) {
    stop("n_axes exceeds the number of available principal components.")
  }

  # Retain the first n_axes of the PCA scores.
  scores <- pca_result$x[, 1:n_axes, drop = FALSE]

  # Compute centroids for each environment.
  env_levels <- levels(env)
  if (length(env_levels) < 2) {
    warning("Less than two environments found; MVPi is undefined. Returning NA.")
    return(NA)
  }
  centroids <- sapply(env_levels, function(x) {
    colMeans(scores[env == x, , drop = FALSE], na.rm = TRUE)
  })
  centroids <- t(centroids)  # rows = environment centroids

  # Compute pairwise Euclidean distances between centroids.
  distances <- as.vector(dist(centroids))

  # Calculate MVPi as the mean Euclidean distance.
  mvpi <- mean(distances, na.rm = TRUE)

  return(mvpi)
}



## test - to be tested with dataset from DOI: 10.1093/jxb/eraa545


#calculate_MVPi(synthetic_data1,1,c(3,4,5),3)


#######################

#' @title Calculate Standardized Plasticity Metric (SPM)
#'
#' @description
#' Computes the Standardized Plasticity Metric (SPM) by comparing mean trait values
#' between a resident and a nonresident environment. SPM quantifies plasticity as the
#' proportional difference relative to the resident environment:
#'
#' \deqn{SPM = \frac{|\,\bar{x}_{resident} - \bar{x}_{nonresident}\,|}{\bar{x}_{resident}}}
#'
#' @param data
#'   A data frame containing the trait measurements and environment labels.
#' @param env_col
#'   Either the name or index of the column in `data` that holds the environment labels,
#'   or a vector of length `nrow(data)` giving environment membership.
#' @param trait_col
#'   The name or index of the column in `data` containing the numeric trait values.
#' @param resident_env
#'   The label or value (matching `env_col`) of the resident environment.
#' @param nonresident_env
#'   The label or value (matching `env_col`) of the nonresident environment.
#'
#' @return
#' A single numeric value: the SPM, i.e. the absolute difference between the mean trait
#' in the resident vs. nonresident environment divided by the resident mean.
#'
#' @examples
#' synthetic_data <- data.frame(
#'   Trait_Values = c(10, 12, 15, 18),
#'   Environment  = factor(c(1, 1, 2, 2))
#' )
#' # Resident = 1, Nonresident = 2
#' spm_value <- calculate_SPM(
#'   data            = synthetic_data,
#'   env_col         = "Environment",
#'   trait_col       = "Trait_Values",
#'   resident_env    = 1,
#'   nonresident_env = 2
#' )
#' print(spm_value)
#'
#' @export
calculate_SPM = function(data, env_col, trait_col, resident_env, nonresident_env) {

  # Handle env_col
  if (is.numeric(env_col) && length(env_col) == 1) {
    env_col = data[[env_col]]
  } else if (is.vector(env_col) && length(env_col) == nrow(data)) {
    env_col = env_col
  } else {
    env_col = data[[env_col]]
  }

  # Extract trait data
  trait_values = data[[trait_col]]

  # Calculate mean trait values for resident and nonresident environments
  trait_resident = mean(trait_values[env_col == resident_env], na.rm = TRUE)
  trait_nonresident = mean(trait_values[env_col == nonresident_env], na.rm = TRUE)

  # Calculate the Standardized Plasticity Metric (SPM)
  spm = abs(trait_resident - trait_nonresident) / trait_resident

  return(spm)
}

#
### test - passed on synthetic dataset
#(mean(df_test2[df_test2[, 1] == 1, 2])-mean(df_test2[df_test2[, 1] == 2, 2]))/mean(df_test2[df_test2[, 1] == 1, 2])
#calculate_SPM(df_test2,1,2,1,2)
#
###test2
#(mean(synthetic_data1[synthetic_data1[, 1] == 1, 3])-mean(synthetic_data1[synthetic_data1[, 1] == 3, 3]))/mean(synthetic_data1[synthetic_data1[, 1] == 1, 3])
#calculate_SPM(synthetic_data1,1,3,1,3)




############################


#' @title Calculate SSpop/SStotal Plasticity Ratio for Multiple Traits
#'
#' @description
#' Computes the plasticity ratio for each specified trait, defined as the
#' proportion of variance attributable to population differences (SSpop)
#' relative to the total phenotypic variance (SStotal) across environments.
#' Internally, a two‐factor ANOVA (population + environment) is fitted and
#' the ratio of the population sum of squares to the total sum of squares is returned.
#'
#' @param data
#'   A data frame containing the trait measurements, environment labels, and population labels.
#' @param env_col
#'   Either the name or index of the column in `data` that holds environment labels,
#'   or an external vector of length `nrow(data)` specifying environment membership.
#' @param trait_cols
#'   A vector of column names or indices in `data` corresponding to the traits
#'   for which to compute the plasticity ratio.
#' @param pop_col
#'   Either the name or index of the column in `data` that holds population labels,
#'   or an external vector of length `nrow(data)` specifying population membership.
#'
#' @return
#' A named numeric vector, one entry per trait in `trait_cols`, giving
#' \eqn{SSpop / SStotal} for that trait.
#'
#' @examples
#' synthetic_data <- data.frame(
#'   Trait_1    = c(10, 12, 15, 18, 11, 17),
#'   Trait_2    = c(20, 22, 25, 28, 21, 27),
#'   Environment = factor(c(1, 1, 2, 2, 3, 3)),
#'   Population  = factor(c("A", "A", "B", "B", "C", "C"))
#' )
#' plasticity_ratios <- calculate_Plasticity_Ratio(
#'   data       = synthetic_data,
#'   env_col    = "Environment",
#'   trait_cols = c("Trait_1", "Trait_2"),
#'   pop_col    = "Population"
#' )
#' print(plasticity_ratios)
#'
#' @export
calculate_Plasticity_Ratio = function(data, env_col, trait_cols, pop_col) {

  # Handle env_col
  if (is.numeric(env_col) && length(env_col) == 1) {
    env_col = data[[env_col]]
  } else if (is.vector(env_col) && length(env_col) == nrow(data)) {
    env_col = env_col
  } else {
    env_col = data[[env_col]]
  }

  # Handle pop_col
  if (is.numeric(pop_col) && length(pop_col) == 1) {
    pop_col = data[[pop_col]]
  } else if (is.vector(pop_col) && length(pop_col) == nrow(data)) {
    pop_col = pop_col
  } else {
    pop_col = data[[pop_col]]
  }

  # Initialize a vector to store Plasticity Ratios for each trait
  plasticity_ratios = numeric(length(trait_cols))
  names(plasticity_ratios) = trait_cols

  # Loop over each trait column
  for (i in seq_along(trait_cols)) {
    # Extract the trait values for the current trait
    trait_values = if (is.numeric(trait_cols[i])) data[[trait_cols[i]]] else data[[trait_cols[i]]]

    # Create a temporary data frame for the ANOVA
    temp_data = data.frame(trait_values, pop_col, env_col)

    # Perform a one-way ANOVA to obtain SSpop and SStotal
    anova_result = aov(trait_values ~ pop_col + env_col, data = temp_data)

    # Extract the sum of squares between populations (SSpop) and total sum of squares (SStotal)
    ss_total = sum(anova(anova_result)[, "Sum Sq"])
    ss_pop = anova(anova_result)["pop_col", "Sum Sq"]

    # Calculate and store the Plasticity Ratio
    plasticity_ratios[i] = ss_pop / ss_total
  }

  return(plasticity_ratios)
}


#pop=c(rep(1,150),rep(2,150))
#calculate_Plasticity_Ratio(synthetic_data1,1,3,pop)
#
### test - passed on synthetic dataset
#
#overall_mean = mean(synthetic_data1[, 3])
#
#mean_pop1 = mean(synthetic_data1[1:150, 3])   # Mean for Population 1 (first 150 rows)
#mean_pop2 = mean(synthetic_data1[151:300, 3])  # Mean for Population 2 (next 150 rows)
#
#n_pop1 = 150  # Number of observations in Population 1
#n_pop2 = 150  # Number of observations in Population 2
#
#ss_pop1 = n_pop1 * (mean_pop1 - overall_mean)^2
#ss_pop2 = n_pop2 * (mean_pop2 - overall_mean)^2
#SS_pop = ss_pop1 + ss_pop2  # Total SS_pop
#
#SS_total = sum((synthetic_data1[, 3] - overall_mean)^2)
#
#plasticity_ratio = SS_pop / SS_total
#print(paste("Plasticity Ratio:", plasticity_ratio))



