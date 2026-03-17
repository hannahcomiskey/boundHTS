#' Monte Carlo Convolution Density for Poisson
#'
#' @param x Numeric evaluation point for the aggregated density.
#' @param lambda_matrix Matrix of lambda parameters for the Poisson
#' distribution for each bottom level observation (rows correspond to simulations).
#'
#' @details
#' This function computes a Monte Carlo approximation of the density of an
#' aggregated random variable formed from a sum of Poisson-distributed
#' components. For each Monte Carlo simulation and posterior draw, the density
#' is evaluated using \code{stats::dpois()}, and the result is averaged
#' across simulations.
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
#' lambda_matrix <- matrix(rnorm(n_sims * n_bottom), nrow = n_sims)
#'
#' Poisson_convolution_density(x = 5, lambda_matrix = lambda_matrix)
#' @importFrom stats dpois
#' @export

Poisson_convolution_density <- function(x, lambda_matrix) {
  n_sims <- dim(lambda_matrix)[1]

  # Compute convoluted sum
  convolution_lambda <- apply(lambda_matrix, 1, sum)

  # Monte Carlo mixture
  avg_over_draws <- mean(stats::dpois(x = x, lambda = convolution_lambda), na.rm = TRUE)

  return(avg_over_draws)
}


