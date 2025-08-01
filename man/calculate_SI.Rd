% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/4.R
\name{calculate_SI}
\alias{calculate_SI}
\title{Stability Index (SI) for a Single Trait}
\usage{
calculate_SI(trait_values, env = NULL)
}
\arguments{
\item{trait_values}{A numeric vector of trait measurements.}

\item{env}{(Optional) A vector (numeric, character, or factor) of the same length as `trait_values`
indicating the environment for each measurement. If `NULL`, environments are assumed
to be equidistant (`1:length(trait_values)`).}
}
\value{
A single numeric value: the Stability Index (SI). Returns `NA` (with a warning)
if the overall mean of environment means is zero.
}
\description{
Computes the Stability Index (SI), quantifying how stable a trait is across environments.
SI is defined as the ratio of the variance among environment means to their overall mean:
\deqn{SI = \frac{\mathrm{Var}(\bar{x}_e)}{\mathrm{Mean}(\bar{x}_e)},}
where \(\bar{x}_e\) are the mean trait values in each environment.
}
\examples{
# Explicit environments:
trait_values <- c(10, 12, 15, 20)
env          <- c(1, 2, 3, 4)
calculate_SI(trait_values, env)

# Equidistant environments assumed:
trait_values <- c(10, 12, 15, 20)
calculate_SI(trait_values)

}
