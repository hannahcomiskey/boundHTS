#' Monte Carlo Convolution Density for Zero-One-Inflated Beta Distributions
#'
#' @param z Numeric evaluation point for the aggregated density.
#' @param alpha_matrix matrix of shape parameters for each element (b) in the aggregate over each of the N observations (N rows x n_nodes)
#' @param beta_matrix matrix of shape parameters for each element in the aggregate (N rows x n_nodes).
#' @param coi_matrix Numeric matrix; Monte Carlo draws of conditional one-inflation probability (n_draws x n_nodes).
#' @param zoi_matrix Numeric matrix; Monte Carlo draws of zero-inflation probability (n_draws x n_nodes).
#' @param weighted_samps matrix of weighted samples for each element in the aggregate (N rows x n_nodes)
#' @param weights vector of weights that are used to combine to form the aggregated density Z. (length n_nodes)
#'
#' @details
#' This function computes a Monte Carlo approximation of the density of an
#' aggregated Zero-One-Inflated Beta (ZOIB) distribution. Posterior draws are passed to \code{dZOIB_4p()},
#' which evaluates the ZOIB density for each Monte Carlo draw.
#' The resulting density is normalised using numerical integration.
#'
#' @return A numeric vector containing the estimated aggregate density evaluated at \code{z_values}.
#'
#' @examples
#' set.seed(1)
#'
#' # Simulation setup
#' n_draws <- 200
#' n_nodes <- 2
#'
#' # Simulated posterior draws
#' Y_mc <- matrix(runif(n_draws * n_nodes, 0.1, 0.9), nrow = n_draws)
#' phi_array <- matrix(rexp(n_draws * n_nodes, 1), nrow = n_draws)
#' zoi_array <- matrix(runif(n_draws * n_nodes, 0, 0.2), nrow = n_draws)
#' coi_array <- matrix(runif(n_draws * n_nodes, 0, 0.2), nrow = n_draws)
#'
#' weights <- c(1, 1)
#' z_values <- seq(0, 2, length.out = 50)
#'
#' dens <- ZOIB_convolution_density(
#'   Y_mc = Y_mc,
#'   phi_array = phi_array,
#'   zoi_array = zoi_array,
#'   coi_array = coi_array,
#'   weights = weights,
#'   z_values = z_values,
#'   n_mc = 100
#' )
#'
#' head(dens)
#'
#' @export

ZOIB_convolution_density <- function(z, alpha_matrix, beta_matrix, coi_matrix,
                                     zi_matrix, weighted_samps, weights) {
  n_draws <- dim(weighted_samps)[1]
  n_mc <- dim(weighted_samps)[2]
  N <- dim(weighted_samps)[3]

  conv_pdf <- matrix(0, nrow = n_mc, ncol=n_draws)

  for(m in 1:n_draws) {
    for(s in 1:n_mc) {
      conv_pdf[s,m] <- dZOIB_4p(x = z - sum(weighted_samps[s,m, 1:(N-1)]),
                               alpha_point = alpha_matrix[s,N],
                               beta_point = beta_matrix[s,N],
                               zi_point = zi_matrix[s,N],
                               coi_point = coi_matrix[s,N],
                               weight = weights[N])
    }
  }
  avg_over_sims <- apply(conv_pdf, 2, mean)

  # Compute final result
  avg_over_draws <- mean(avg_over_sims, na.rm=TRUE)

  return(avg_over_draws)

}
