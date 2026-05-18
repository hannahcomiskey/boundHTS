#' Monte Carlo Convolution Density for Gamma Distributions
#'
#' @param z Numeric evaluation point for the aggregated density.
#' @param shape_point Vector of shape parameters for the Gamma distribution for each observation.
#' @param rate_point Matrix of rate parameters for the Gamma distribution for each observation.
#' @param bottom_samps Array of bottom series samples for each element contributing to the
#'   aggregate (dimensions: n_sims x n_draws x b).
#'
#' @details
#' This function computes a Monte Carlo approximation of the density of an
#' aggregated random variable formed from a sum of Gamma distributed
#' components. For each Monte Carlo simulation and point estimate, the density
#' is evaluated using \code{stats::rgamma()}, and the result is averaged
#' across simulations.
#'
#' @return A numeric value representing the estimated aggregated density evaluated at \code{z}.
#'
#' @examples
#' set.seed(1)
#'
#' # Simulation setup
#' n_sims <- 50
#' n_draws <- 10
#' b <- 2
#'
#' # Simulated weighted samples
#' bottom_samps <- array(runif(n_sims * n_draws * b),
#'                         dim = c(n_sims, n_draws, b))
#'
#' shape_point <- runif(2, 2, 5)
#' rate_point  <- runif(2, 2, 5)
#'
#' weights <- c(1, 1)
#'
#' Gamma_convolution_density_point(
#'   z = 0.5,
#'   shape_point = shape_point,
#'   rate_point = rate_point,
#'   bottom_samps = bottom_samps
#' )
#'
#' @export

Gamma_convolution_density_point <- function(z, shape_point, rate_point, bottom_samps) {
  N <- dim(bottom_samps)[3]

  # Sum simulated values from first N-1 components
  partial_sum <- apply(bottom_samps[, , 1:(N - 1), drop = FALSE], c(1, 2), sum)

  # Remaining value for final gamma density
  x <- z - partial_sum

  # Gamma density is zero for x <= 0
  dens <- stats::dgamma(x, shape = shape_point[N], rate  = rate_point[N])

  # Monte Carlo estimate of convolution density
  Density <- mean(dens, na.rm = TRUE)

  return(Density)
}
