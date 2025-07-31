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
calculate_rdpi = function(trait_values, env_values = NULL) {
  # Input validation
  if (!is.numeric(trait_values)) {
    stop("trait_values must be numeric")
  }

  n_envs = length(trait_values)

  # If no environment values provided, create sequential environments
  if (is.null(env_values)) {
    env_values = seq_len(n_envs)
  }

  # Ensure same length
  if (length(trait_values) != length(env_values)) {
    stop("trait_values and env_values must have the same length")
  }

  # Get all pairs of environment indices
  env_pairs = combn(n_envs, 2)
  n_pairs = ncol(env_pairs)

  # Calculate RDPIs for each pair
  rdpis = numeric(n_pairs)

  for (i in seq_len(n_pairs)) {
    idx1 = env_pairs[1, i]
    idx2 = env_pairs[2, i]

    # Calculate relative distance for this pair
    abs_diff = abs(trait_values[idx2] - trait_values[idx1])
    sum_vals = trait_values[idx1] + trait_values[idx2]

    # Handle potential division by zero
    if (sum_vals == 0) {
      rdpis[i] = 0
    } else {
      rdpis[i] = abs_diff / sum_vals
    }
  }

  # Calculate final RDPI as mean of all RDPIs
  rdpi = mean(rdpis, na.rm = TRUE)

  return(rdpi)
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





#' @title Relative Distance Plasticity Index for Mean Phenotypic Values (RDPI)
#'
#' @description
#' Computes the Relative Distance Plasticity Index (RDPI) for mean phenotypic values of one or more traits.
#' RDPI for a given trait and genotype is defined as the absolute difference between the mean trait values
#' in two environments, divided by the smaller of those two means.  When multiple environments are present,
#' all pairwise combinations are evaluated and the results collated.
#'
#' @param dataframe
#'   A data.frame containing your phenotypic data.
#' @param trait_cols
#'   A character vector or numeric indices identifying the column(s) in `dataframe` that hold trait values.
#' @param sp
#'   (Optional) A column name or index giving the genotype (or species) identifier.  If `NULL` (default),
#'   the entire data set is treated as a single group.
#' @param factors
#'   (Optional) Column name(s) or indices of internal environmental factors in `dataframe`.  These will be
#'   combined (via `interaction()`) to define unique environments.
#' @param factors_not_in_dataframe
#'   (Optional) A named list of external factor vectors (all of the same length as `nrow(dataframe)`)
#'   to be combined with the internal `factors` into a single `Combined_Factors` column.
#' @param stat_analysis
#'   (Optional) Logical.  If `TRUE`, performs ANOVA and Tukey HSD on each trait across the environmental
#'   combinations, and generates boxplots.  Default is `NULL` (no statistical tests).
#'
#' @return
#' A list with components:
#' \describe{
#'   \item{rdpi_results}{A data.frame of per‑genotype, per‑trait, per‑environment‑pair RDPI values.}
#'   \item{trait_boxplots}{A ggplot2 object of boxplots (only if `stat_analysis = TRUE`).}
#'   \item{anova_results}{A list of ANOVA summaries for each trait (only if `stat_analysis = TRUE`).}
#'   \item{tukey_results}{A list of Tukey HSD test results for each trait (only if `stat_analysis = TRUE`).}
#' }
#'
#' @details
#' 1. Internal and external factors are combined via `interaction()` to form `Combined_Factors`.
#' 2. For each genotype (or overall if `sp = NULL`), and each trait, all pairwise environment means
#'    are computed and RDPI = |mean₁ − mean₂| / min(mean₁, mean₂).
#' 3. If `stat_analysis = TRUE`, an ANOVA (`aov()`) and Tukey HSD (`agricolae::HSD.test()`) are run for
#'    each trait across `Combined_Factors`, and boxplots produced via ggplot2.
#'
#' @examples
#' df <- data.frame(
#'   Genotype   = rep(c("G1","G2"), each = 6),
#'   Environment = rep(c("E1","E2","E3"), times = 4),
#'   Height     = c(10,12,11,14,16,15, 20,22,21,23,25,24),
#'   Biomass    = rnorm(12, 5, 1)
#' )
#'
#' # Compute RDPI for Height, by Genotype and Environment:
#' res <- rdpi_mean_calculation(
#'   dataframe   = df,
#'   trait_cols  = "Height",
#'   sp          = "Genotype",
#'   factors     = "Environment"
#' )
#' head(res$rdpi_results)
#'
#' # With ANOVA/Tukey and boxplots:
#' res2 <- rdpi_mean_calculation(
#'   dataframe     = df,
#'   trait_cols    = "Height",
#'   sp            = "Genotype",
#'   factors       = "Environment",
#'   stat_analysis = TRUE
#' )
#' print(res2$trait_boxplots)
#'
#' @export
rdpi_mean_calculation = function(dataframe, trait_cols, sp = NULL, factors = NULL, factors_not_in_dataframe = NULL, stat_analysis = NULL) {

  # Convert column indices to names if necessary
  if (!is.null(sp)) {
    sp = if (is.numeric(sp)) names(dataframe)[sp] else sp
  }
  traits = if (is.numeric(trait_cols)) names(dataframe)[trait_cols] else trait_cols

  # Combine internal and external factors into a single dataframe
  if (!is.null(factors_not_in_dataframe)) {
    if (length(factors_not_in_dataframe[[1]]) != nrow(dataframe)) {
      stop("The length of external factors must match the number of rows in the dataframe.")
    }
    external_factors_df = as.data.frame(factors_not_in_dataframe)
    if (!is.null(factors)) {
      factors = if (is.numeric(factors)) names(dataframe)[factors] else factors
      combined_factors_df = cbind(dataframe[factors], external_factors_df)
    } else {
      combined_factors_df = external_factors_df
    }
    dataframe$Combined_Factors = interaction(combined_factors_df, drop = TRUE)
  } else if (!is.null(factors)) {
    factors = if (is.numeric(factors)) names(dataframe)[factors] else factors
    dataframe$Combined_Factors = interaction(dataframe[factors], drop = TRUE)
  } else {
    stop("You must provide either internal factors, external factors, or both.")
  }

  dataframe$Combined_Factors = as.factor(dataframe$Combined_Factors)

  all_rdpi_results = data.frame()  # Initialize a dataframe to store all RDPI values

  if (is.null(sp)) {
    unique_species = list("Single_Group" = dataframe)
  } else {
    unique_species = split(dataframe, dataframe[[sp]])
  }

  for (species_name in names(unique_species)) {
    species_data = unique_species[[species_name]]

    for (trait in traits) {
      RDPI_values = data.frame(sp = character(), env1 = character(), env2 = character(), rdpi = numeric())

      # Get unique environment combinations
      env_levels = levels(species_data$Combined_Factors)
      n_env_levels = length(env_levels)

      mean_values = aggregate(species_data[[trait]], list(species_data$Combined_Factors), mean)
      colnames(mean_values) = c("Env_Combination", "Mean_Trait")

      for (i in 1:(n_env_levels - 1)) {
        for (j in (i + 1):n_env_levels) {
          mean_i = mean_values$Mean_Trait[mean_values$Env_Combination == env_levels[i]]
          mean_j = mean_values$Mean_Trait[mean_values$Env_Combination == env_levels[j]]

          # Calculate RDPI for this trait between the two environments
          rdpi_value = abs(mean_i - mean_j) / min(mean_i, mean_j)

          # Append the RDPI value for this species, trait, and environment pair
          RDPI_values = rbind(RDPI_values, data.frame(sp = species_name, env1 = env_levels[i], env2 = env_levels[j], rdpi = rdpi_value))
        }
      }

      all_rdpi_results = rbind(all_rdpi_results, RDPI_values)
    }
  }

  # If statistical analysis is requested, perform ANOVA and Tukey's HSD
  if (!is.null(stat_analysis)) {
    all_trait_data = data.frame()

    for (species_name in names(unique_species)) {
      species_data = unique_species[[species_name]]

      for (trait in traits) {
        trait_data = data.frame(
          sp = species_name,
          trait = trait,
          Combined_Factors = species_data$Combined_Factors,
          Trait_Value = species_data[[trait]]
        )

        all_trait_data = rbind(all_trait_data, trait_data)
      }
    }

    # Perform ANOVA and Tukey's HSD test
    anova_results = list()
    tukey_results = list()

    for (trait in traits) {
      fit = aov(Trait_Value ~ Combined_Factors, data = subset(all_trait_data, trait == trait))
      anova_results[[trait]] = summary(fit)

      # Perform Tukey's HSD test
      Tukey = agricolae::HSD.test(fit, trt = "Combined_Factors")
      tukey_results[[trait]] = Tukey
    }

    # Create boxplots for the traits across environmental combinations
    boxplot_traits = ggplot2::ggplot(all_trait_data, ggplot2::aes(x = Combined_Factors, y = Trait_Value, fill = trait)) +
      ggplot2::geom_boxplot() +
      ggplot2::facet_wrap(~trait, scales = "free_y") +
      ggplot2::theme_bw() +
      ggplot2::xlab("Environmental Combinations") +
      ggplot2::ylab("Trait Values") +
      ggplot2::ggtitle("Boxplots of Trait Values Across Environmental Combinations") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

    return(list(
      rdpi_results = all_rdpi_results,
      trait_boxplots = boxplot_traits,
      anova_results = anova_results,
      tukey_results = tukey_results
    ))
  }

  # If stat_analysis is NULL, return only the RDPI values
  return(all_rdpi_results)
}




## Example usage with synthetic data
#external_light = rep(c(0.4, 0.6, 0.8), each = 100)
#external_water =rep(c("Low", "High"), each = 150)
#
### test - passed on synthetic dataset
#
#rdpi_mean_calculation(df_test20,sp=1, trait_cols=3,factors = 2)


#########################################


#' @title Environmental Sensitivity Performance Index (ESPI)
#'
#' @description
#' Calculates the Environmental Sensitivity Performance Index (ESPI) for one or more traits
#' across different environmental conditions. ESPI quantifies the sensitivity of a trait by
#' comparing the difference between its maximum and minimum mean values to the range of the
#' environmental gradient.
#'
#' ESPI is defined as:
#' \deqn{ESPI = \frac{\max(\bar{z}_e) - \min(\bar{z}_e)}{\lvert \max(e) - \min(e)\rvert},}
#' where \(\bar{z}_e\) is the mean trait value in environment \(e\) and \(e\) is the environmental value.
#'
#' @param trait_values
#'   A numeric vector for a single trait, or a data frame with one column per trait.
#' @param env_values
#'   (Optional) A numeric vector of environmental values corresponding to `trait_values`. If `NULL`,
#'   environments are assumed to be equidistant (1, 2, …).
#'
#' @return
#' For a single trait vector, returns a single numeric ESPI value. For a data frame of traits,
#' returns a named numeric vector of ESPI values (one per column).
#'
#' @examples
#' # Single trait with explicit environments
#' trait <- c(10, 15, 20, 25, 30)
#' env   <- c(1, 2, 3, 4, 5)
#' calculate_ESPI(trait, env)
#'
#' # Multiple traits in a data frame
#' df <- data.frame(
#'   Height = c(10, 15, 20, 25, 30),
#'   Weight = c(5,  7,  9, 11, 13)
#' )
#' calculate_ESPI(df, env)
#'
#' # Without specifying env_values (assumes equidistant)
#' calculate_ESPI(trait)
#'
#' @export
calculate_ESPI = function(trait_values, env_values = NULL) {

  if (is.null(env_values)) {
    env_values = seq_along(trait_values)

  }

  if (!is.vector(trait_values) && !is.data.frame(trait_values)) {
    stop("trait_values must be a numeric vector or a data frame where each column is a trait.")
  }

  if (length(env_values) != length(trait_values)) {
    stop("env_values must be the same length as trait_values.")
  }

  # Function to calculate ESPI for a single trait
  calculate_single_espi = function(single_trait) {
    means = tapply(single_trait, env_values, mean, na.rm = TRUE)
    max_mean = max(means, na.rm = TRUE)
    min_mean = min(means, na.rm = TRUE)

    abs_env_distance = abs(max(as.numeric(env_values), na.rm = TRUE) - min(as.numeric(env_values), na.rm = TRUE))
    if (abs_env_distance > 0) {
      return((max_mean - min_mean) / abs_env_distance)
    } else {
      return(NA)
    }
  }

  # Handle single trait
  if (is.vector(trait_values)) {
    return(calculate_single_espi(trait_values))
  } else if (is.data.frame(trait_values)) {
    # Handle multiple traits
    espi_results = sapply(trait_values, calculate_single_espi)
    names(espi_results) = colnames(trait_values)
    return(espi_results)
  }
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

  if (!is.numeric(trait_values)) {
    stop("`trait_values` must be a numeric vector.")
  }
  n <- length(trait_values)
  if (n < 2) {
    stop("At least two trait values are required.")
  }

  # Equidistant numeric environments
  env_values <- factor(seq_len(n))

  env_levels <- levels(env_values)
  m <- length(env_levels)
  if (m < 2) {
    stop("At least two environments are required.")
  }

  espiid_list <- numeric()
  names(espiid_list) <- character()

  for (i in seq_len(m - 1)) {
    for (j in (i + 1):m) {
      # values in env i and j
      zi <- trait_values[env_values == env_levels[i]]
      zj <- trait_values[env_values == env_levels[j]]

      # all pairwise absolute differences
      diffs <- abs(outer(zi, zj, "-"))
      dist_ij <- if (use_median) median(diffs, na.rm = TRUE) else mean(diffs, na.rm = TRUE)

      env_dist <- abs(as.numeric(env_levels[i]) - as.numeric(env_levels[j]))
      espiid_ij <- if (env_dist > 0) dist_ij / env_dist else NA_real_

      pair_name <- paste0(env_levels[i], "-", env_levels[j])
      espiid_list[pair_name] <- espiid_ij
    }
  }

  if (aggregate) {
    mean(espiid_list, na.rm = TRUE)
  } else {
    espiid_list
  }
}

#l=list(external_light)
#
#external_factors=list(external_light,external_water)
#
### test - passed on synthetic dataset
#
#espiid_calculation(df_test2.2, trait_cols=2,factors=1)




