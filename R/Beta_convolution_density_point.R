#' Monte Carlo Convolution Density for Four-Parameter Beta Distributions
#'
#' @param z Numeric evaluation point for the aggregated density.
#' @param alpha_point Point estimates of alpha (shape1) parameters for the Beta
#'   distribution for each observation.
#' @param beta_point Point estimates of beta (shape2) parameters for the Beta
#'   distribution for each observation.
#' @param weighted_samps Array of weighted samples for each element in the
#'   aggregate (dimensions: n_sims x n_draws x b).
#' @param weights Numeric vector of weights used to combine the components
#'   into the aggregated density \eqn{Z} (length b).
#'
#' @details
#' This function computes a Monte Carlo approximation of the density of an
#' aggregated random variable formed from a weighted sum of Beta-distributed
#' components. For each Monte Carlo simulation and point estimate, the density
#' is evaluated using \code{ExtDist::dBeta_ab()}, and the result is averaged
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
#' weighted_samps <- array(runif(n_sims * n_draws * b),
#'                         dim = c(n_sims, n_draws, b))
#'
#' alpha_point <- matrix(runif(n_sims * b, 2, 5), nrow = n_sims)
#' beta_point  <- matrix(runif(n_sims * b, 2, 5), nrow = n_sims)
#'
#' weights <- c(1, 1)
#'
#' Beta_convolution_density(
#'   z = 0.5,
#'   alpha_point = alpha_point,
#'   beta_point = beta_point,
#'   weighted_samps = weighted_samps,
#'   weights = weights
#' )
#'
#' @export

Beta_convolution_density_point <- function(z, alpha_point, beta_point, weighted_samps, weights) {
  N <- dim(weighted_samps)[3]
  n_sims <- dim(weighted_samps)[1]
  n_draws <- dim(weighted_samps)[2]
  conv_pdf <- matrix(0, nrow = n_sims, ncol=n_draws)

  for(m in 1:n_draws) {
    for(s in 1:n_sims) {
      conv_pdf[s,m] <- ExtDist::dBeta_ab(x = z - sum(weighted_samps[s,m, 1:(N-1)]),
                                         alpha = alpha_point[N],
                                         beta = beta_point[N],
                                         lower=0,
                                         upper=weights[N])
    }
  }
  avg_over_sims <- apply(conv_pdf, 2, mean)

  # Compute final result
  avg_over_draws <- mean(avg_over_sims, na.rm=TRUE)

  return(avg_over_draws)
}
