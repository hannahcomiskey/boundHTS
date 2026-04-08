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
#' n_sims <- 10
#' n_draws <- 200
#' N <- 2
#'
#' # Simulated posterior draws
#'  weighted_samps <- array(runif(n_sims * n_draws * N, min = 0, max = 0.2),
#'  dim = c(n_sims, n_draws, N))
#'  alpha_matrix <- matrix(runif(n_draws * N, min = 2, max = 10),
#'  nrow = n_draws, ncol = N)
#'  beta_matrix <- matrix(runif(n_draws * N, min = 2, max = 10),
#'  nrow = n_draws, ncol = N)
#'  zoi_matrix <- matrix(runif(n_draws * N, 0.1, 0.4), nrow = n_draws, ncol = N)
#'  coi_matrix <- matrix(runif(n_draws * N, 0.2, 0.8), nrow = n_draws, ncol = N)
#'
#' weights <- c(1, 1)
#' z_values <- seq(0, 2, length.out = 50)
#'
#' dens <- ZOIB_convolution_density(z=0.5, alpha_matrix=alpha_matrix,
#' beta_matrix = beta_matrix, zoi_matrix = zoi_matrix, coi_matrix = coi_matrix,
#' weights = weights)
#'
#' head(dens)
#'
#' @export

ZOIB_convolution_density <- function(z, alpha_matrix, beta_matrix, coi_matrix,
                                     zoi_matrix, weighted_samps, weights) {
  n_sims <- dim(weighted_samps)[1]
  n_draws <- dim(weighted_samps)[2]
  N <- dim(weighted_samps)[3]

  conv_pdf <- matrix(0, nrow = n_sims, ncol=n_draws)

  for(s in 1:n_sims) {
    for(m in 1:n_draws) {
      conv_pdf[s,m] <- dZOIB_4p(x = z - sum(weighted_samps[s,m, 1:(N-1)]),
                               alpha_point = alpha_matrix[m,N],
                               beta_point = beta_matrix[m,N],
                               zi_point = zoi_matrix[m,N],
                               coi_point = coi_matrix[m,N],
                               weight = weights[N])
    }
  }
  avg_over_sims <- apply(conv_pdf, 2, mean)

  # Compute final result
  avg_over_draws <- mean(avg_over_sims, na.rm=TRUE)

  return(avg_over_draws)

}
