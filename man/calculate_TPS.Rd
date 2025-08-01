% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/3.R
\name{calculate_TPS}
\alias{calculate_TPS}
\title{Calculate Transplant Plasticity Score (TPS)}
\usage{
calculate_TPS(trait_values, env_values, native_env, transplanted_env)
}
\arguments{
\item{trait_values}{A numeric vector containing trait values measured across different environments.}

\item{env_values}{A numeric vector specifying the environment associated with each trait measurement.
If NULL, environments are assumed to be equidistant.}

\item{native_env}{A value indicating the native environment. Must be present in `env_values`.}

\item{transplanted_env}{A value indicating the transplanted environment. Must be present in `env_values`.}
}
\value{
A numeric value representing the mean TPS for the trait.
}
\description{
This function calculates the Transplant Plasticity Score (TPS) for a given trait when an organism is transplanted to a different environment.
The TPS quantifies the relative change in a trait between the native environment and the transplanted environment.
}
\examples{
trait_values = c(50, 60, 70, 80)
env_values = c(1, 1, 2, 2)  # Native = 1, Transplanted = 2

# Calculate TPS between native (1) and transplanted (2)
tps_result = calculate_TPS(trait_values, env_values, native_env = 1, transplanted_env = 2)
print(tps_result)

}
