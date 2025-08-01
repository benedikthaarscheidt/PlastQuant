% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/3.R
\name{calculate_PR}
\alias{calculate_PR}
\title{Calculate Phenotypic Range (PR)}
\usage{
calculate_PR(trait_values, env_values = NULL, across = TRUE)
}
\arguments{
\item{trait_values}{A numeric vector representing the measured trait values.}

\item{env_values}{(Optional) A numeric vector representing the environmental conditions.
If NULL, equidistant environments are assumed.}

\item{across}{A logical value. If TRUE, calculates the PR across all environments combined.}
}
\value{
A numeric value representing the Phenotypic Range (PR).
}
\description{
Computes the Phenotypic Range (PR) for a trait across different environments.
}
\details{
The Phenotypic Range (PR) is defined as:

\deqn{PR = \max(trait\_values) - \min(trait\_values)}

If `across = TRUE`, PR is calculated across **all** environments.

If `env_values` is **not provided**, equidistant environments are assumed.
}
\examples{
# Example 1: Using explicitly defined environmental values
trait_values = c(10, 12, 15, 17, 20)
env_values = c(5, 10, 15, 20, 25)
calculate_PR(trait_values, env_values)

# Example 2: Using equidistant environments (default behavior)
trait_values = c(10, 15, 20, 25, 30)
calculate_PR(trait_values)

# Example 3: Calculating PR across all environments
calculate_PR(trait_values, env_values, across = TRUE)

}
