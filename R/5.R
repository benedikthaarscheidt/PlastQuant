


#' @title Calculate Plasticity Score from a Polynomial Reaction Norm
#'
#' @description
#' Fits polynomial regression models of degrees 1 through \code{max_degree} to describe the relationship
#' between a phenotype and an environmental gradient. The best-fitting model is selected based on the
#' specified information criterion (\code{"AIC"} or \code{"BIC"}), with preference for lower-degree models
#' when criterion differences are small. The plasticity score is then computed as the sum of the absolute
#' values of the estimated polynomial coefficients (excluding the intercept).
#'
#' @param trait_values
#'   A numeric vector of phenotype values measured across environments.
#' @param env_values
#'   (Optional) A numeric vector of environmental values corresponding to \code{trait_values}.
#'   If \code{NULL} (default), environments are assumed equidistant and set to \code{1:length(trait_values)}.
#' @param max_degree
#'   Integer ≥ 1. The maximum polynomial degree to consider in model fitting. Default is 3.
#' @param criterion
#'   A character string, either \code{"AIC"} or \code{"BIC"}, indicating which information criterion
#'   to use for model selection. Default is \code{"BIC"}.
#'
#' @return
#' A named list with components:
#' \describe{
#'   \item{\code{best_degree}}{The degree of the selected polynomial model.}
#'   \item{\code{plasticity_score}}{The sum of absolute values of the polynomial coefficients (excluding intercept).}
#'   \item{\code{coefficients}}{A numeric vector of the estimated polynomial coefficients (excluding intercept).}
#' }
#'
#' @examples
#' trait_values <- c(100, 110, 120, 105, 115, 125)
#' env_values   <- c(1, 2, 3, 1, 2, 3)
#'
#' # Default: BIC-based selection
#' result1 <- calculate_plasticity(trait_values, env_values)
#' print(result1)
#'
#' # AIC-based selection, equidistant environments
#' result2 <- calculate_plasticity(trait_values, max_degree = 4, criterion = "AIC")
#' print(result2)
#'
#' @references
#' Scheiner, S. M. (1993). Genetics and Evolution of Phenotypic Plasticity. Annual Review of Ecology and Systematics, 24, 35–68.
#'
#' @export
calculate_plasticity = function(trait_values, env_values = NULL, max_degree = 3, criterion = "BIC") {
  # Input validation
  if (!is.numeric(trait_values)) {
    stop("trait_values must be a numeric vector")
  }

  n_values = length(trait_values)

  # If no environmental values are provided, create equidistant values
  if (is.null(env_values)) {
    env_values = seq_len(n_values)
  }

  # Ensure environmental values are numeric
  if (!is.numeric(env_values)) {
    stop("env_values must be numeric")
  }

  # Check for length mismatch
  if (length(trait_values) != length(env_values)) {
    stop("trait_values and env_values must have the same length")
  }

  # Store best model information
  best_degree = 1
  best_criterion_value = Inf
  best_model = NULL

  # Try polynomial degrees from 1 to max_degree
  for (degree in 1:max_degree) {
    formula = as.formula(paste("trait_values ~ poly(env_values, ", degree, ", raw = TRUE)", sep = ""))
    model = lm(formula)

    # Use AIC or BIC for model selection
    model_criterion = ifelse(criterion == "AIC", AIC(model), BIC(model))

    # Penalize higher-degree models if difference is small (<2)
    if ((model_criterion) < best_criterion_value) {
      best_criterion_value = model_criterion
      best_degree = degree
      best_model = model
    }
  }


  coefficients = coef(best_model)[-1]

  # Compute plasticity score as sum of absolute values of coefficients
  plasticity_score = sum(abs(coefficients))

  return(list(
    best_degree = best_degree,
    plasticity_score = plasticity_score,
    coefficients = coefficients
  ))
}

##############################################


#' @title Cross‐Environment Covariance and Correlation for Linear Reaction Norms
#'
#' @description
#' Computes the covariance of a trait measured across environments, and optionally its Pearson correlation
#' with the environmental gradient. If \code{env_values} is not supplied, environments are assumed
#' equidistant (1, 2, …, length of trait_values).
#'
#' @param trait_values
#'   A numeric vector of trait measurements taken across environments.
#' @param env_values
#'   (Optional) A numeric vector of the same length as \code{trait_values}, giving the environment
#'   associated with each measurement. If \code{NULL}, defaults to \code{seq_along(trait_values)}.
#' @param return_correlation
#'   Logical; if \code{TRUE}, returns both covariance and correlation (default is \code{FALSE}).
#'
#' @return
#' A named list containing:
#' \itemize{
#'   \item{\code{covariance}}{Numeric, the covariance between \code{trait_values} and \code{env_values}.}
#'   \item{\code{correlation}}{(Only if \code{return_correlation = TRUE}) Numeric, the Pearson correlation.}
#' }
#'
#' @examples
#' trait_values <- c(10, 12, 15, 17, 20, 22, 25, 27, 30, 32)
#'
#' # 1) Covariance with equidistant environments
#' cross_env_cov(trait_values)
#'
#' # 2) Covariance and correlation with explicit environments
#' env_values <- seq_along(trait_values)
#' cross_env_cov(trait_values, env_values, return_correlation = TRUE)
#'
#' @export
cross_env_cov = function(trait_values, env_values = NULL, return_correlation = FALSE) {
  # Input validation
  if (!is.numeric(trait_values)) {
    stop("trait_values must be a numeric vector")
  }

  n = length(trait_values)

  # If no env_values are provided, assume equidistant environments
  if (is.null(env_values)) {
    env_values = seq_len(n)  # Equidistant environments
  }

  # Ensure correct length
  if (length(env_values) != n) {
    stop("trait_values and env_values must have the same length")
  }

  # Calculate covariance
  covariance = cov(trait_values, env_values)

  # Calculate correlation if requested
  if (return_correlation) {
    correlation = cor(trait_values, env_values)
    return(list(covariance = covariance, correlation = correlation))
  } else {
    return(list(covariance = covariance))
  }
}

