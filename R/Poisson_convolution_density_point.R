#' Convolution Density for Poisson
#'
#' @param x Numeric evaluation point for the aggregated density.
#' @param lambda_point Vector of lambda parameters for the Poisson
#' distribution for each bottom level observation.
#'
#' @details
#' This function computes an approximation of the density of an
#' aggregated random variable formed from a sum of Poisson-distributed
#' components. The density is evaluated using \code{stats::dpois()}.
#'
#' @return A numeric value representing the estimated aggregated density evaluated at \code{z}.
#'
#' @examples
#' set.seed(1)
#'
#' # Simulation setup
#' n_sims <- 50
#' n_bottom <- 5
#'
#' # Simulated lambda samples
#' lambda_vector <- rnorm(n_bottom)
#'
#' Poisson_convolution_density_point(x = 5, lambda_vector = lambda_vector)
#' @import stats dpois
#' @export

Poisson_convolution_density_point <- function(x, lambda_vector) {

  # Compute convoluted sum
  convolution_lambda <- sum(lambda_vector)

  # Density of mixture
  dens <- stats::dpois(x = x, lambda = convolution_lambda)

  return(dens)
}


