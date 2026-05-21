#' Monte Carlo Convolution Density for Gamma Distributions
#'
#' @param z Numeric evaluation point for the aggregated density.
#' @param shape_matrix Matrix of shape parameters for the Gamma
#'   distribution for each observation (rows correspond to draws, columns correspond to nodes).
#' @param rate_matrix Matrix of rate parameters for the Gamma
#'   distribution for each observation (rows correspond to draws, columns correspond to nodes).
#' @param bottom_samps Array of bottom series samples for each element contributing to the
#'   aggregate (dimensions: n_sims x n_draws x b).
#'
#' @details
#' This function computes a Monte Carlo approximation of the density of an
#' aggregated random variable formed from a sum of Gamma-distributed
#' components. For each Monte Carlo simulation and posterior draw, the density
#' is evaluated using \code{stats::rgamma()}, and the result is averaged
#' across simulations and posterior draws.
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
#' # Simulated samples
#' bottom_samps <- array(rgamma(n_sims * n_draws * b, shape = 2, rate= 3),
#'  dim = c(n_sims, n_draws, b))
#'
#' shape_matrix <- matrix(runif(n_draws * b, 2, 5), nrow = n_draws)
#' rate_matrix  <- matrix(runif(n_draws * b, 2, 5), nrow = n_draws)
#'
#' weights <- c(1, 1)
#'
#' Gamma_convolution_density(z = 0.5,
#'   shape_matrix = shape_matrix,
#'   rate_matrix = rate_matrix,
#'   bottom_samps = bottom_samps
#' )
#'
#' @export

Gamma_convolution_density <- function(z, shape_matrix, rate_matrix, bottom_samps) {
  N <- dim(bottom_samps)[3]
  n_draws <- dim(bottom_samps)[2]
  n_sims <- dim(bottom_samps)[1]
  conv_pdf <- matrix(0, nrow = n_sims, ncol=n_draws)

  # Generate convolution
  for(s in 1:n_sims) {
    for(m in 1:n_draws) {
      conv_pdf[s,m] <- stats::dgamma(x = z - sum(bottom_samps[s,m, 1:(N-1)]),
                                     shape = shape_matrix[m,N],
                                     rate = rate_matrix[m,N])
    }
  }
  avg_over_sims <- apply(conv_pdf, 2, mean)

  # Compute final result
  avg_over_draws <- mean(avg_over_sims, na.rm=TRUE)

  return(avg_over_draws)
}
