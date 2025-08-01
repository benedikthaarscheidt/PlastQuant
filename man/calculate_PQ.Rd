% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/3.R
\name{calculate_PQ}
\alias{calculate_PQ}
\title{Calculate Plasticity Quotient (PQ)}
\usage{
calculate_PQ(trait_values, env_values = NULL)
}
\arguments{
\item{trait_values}{A numeric vector representing the measured trait values across environments.}

\item{env_values}{(Optional) A numeric vector representing the environmental conditions
in which the traits were measured. If NULL, equidistant environments are assumed.}
}
\value{
A numeric value representing the Plasticity Quotient (PQ).
}
\description{
Computes the Plasticity Quotient (PQ) based on the range of trait values
and the corresponding environmental factor range.
}
\details{
The Plasticity Quotient (PQ) is calculated as:

\deqn{PQ = \frac{\max(trait\_values) - \min(trait\_values)}{\max(env\_values) - \min(env\_values)}}

If `env_values` is not provided, it is assumed that the environments are equidistant (i.e.,
`env_values = 1, 2, 3, ..., n` where `n` is the number of trait values).

If the environmental range is zero (i.e., all values are identical), a warning is issued and `NA` is returned.
}
\examples{
# Example 1: Using explicitly defined environmental values
trait_values = c(10, 12, 15, 17, 20)
env_values = c(5, 10, 15, 20, 25)
calculate_PQ(trait_values, env_values)

# Example 2: Using equidistant environments (default behavior)
trait_values = c(10, 15, 20, 25, 30)
calculate_PQ(trait_values)

# Example 3: Handling zero environment range (returns NA with warning)
trait_values = c(10, 15, 20, 25, 30)
env_values = rep(5, 5)
calculate_PQ(trait_values, env_values)

}
