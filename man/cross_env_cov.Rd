% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/5.R
\name{cross_env_cov}
\alias{cross_env_cov}
\title{Cross‐Environment Covariance and Correlation for Linear Reaction Norms}
\usage{
cross_env_cov(trait_values, env_values = NULL, return_correlation = FALSE)
}
\arguments{
\item{trait_values}{A numeric vector of trait measurements taken across environments.}

\item{env_values}{(Optional) A numeric vector of the same length as \code{trait_values}, giving the environment
associated with each measurement. If \code{NULL}, defaults to \code{seq_along(trait_values)}.}

\item{return_correlation}{Logical; if \code{TRUE}, returns both covariance and correlation (default is \code{FALSE}).}
}
\value{
A named list containing:
\itemize{
  \item{\code{covariance}}{Numeric, the covariance between \code{trait_values} and \code{env_values}.}
  \item{\code{correlation}}{(Only if \code{return_correlation = TRUE}) Numeric, the Pearson correlation.}
}
}
\description{
Computes the covariance of a trait measured across environments, and optionally its Pearson correlation
with the environmental gradient. If \code{env_values} is not supplied, environments are assumed
equidistant (1, 2, …, length of trait_values).
}
\examples{
trait_values <- c(10, 12, 15, 17, 20, 22, 25, 27, 30, 32)

# 1) Covariance with equidistant environments
cross_env_cov(trait_values)

# 2) Covariance and correlation with explicit environments
env_values <- seq_along(trait_values)
cross_env_cov(trait_values, env_values, return_correlation = TRUE)

}
