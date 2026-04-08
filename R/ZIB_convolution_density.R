#' Monte Carlo Estimate of Aggregated Zero-Inflated Four-Parameter Beta (ZIB) Density
#'
#' @param z Numeric evaluation point for the aggregated density.
#' @param alpha_matrix matrix of shape parameters for each element (b) in the aggregate over each of the N observations (N rows x n_nodes)
#' @param beta_matrix matrix of shape parameters for each element in the aggregate (N rows x n_nodes)
#' @param zi_matrix Numeric matrix; Monte Carlo draws of zero-inflation probability (n_draws x n_nodes).
#' @param weighted_samps matrix of weighted samples for each element in the aggregate (N rows x n_nodes)
#' @param weights vector of weights that are used to combine to form the aggregated density Z. (length n_nodes)
#'
#' @details
#' For each evaluation point `z`, a Monte Carlo estimate of the aggregated zero-inflated Beta
#' density is computed by resampling posterior draws from `alpha_matrix`, `beta_matrix,` and `zi_matrix` and
#' averaging `dZIB_4p` across the selected draws. The resulting density is normalized
#' using the trapezoidal rule (`pracma::trapz`) to ensure integration to 1.
#'
#' @return Numeric vector of the same length as `z_values` representing the normalized
#'   aggregated zero-inflated Beta density.
#'
#' @examples
#' set.seed(1)
#'
#' n_draws <- 5
#' n_nodes <- 2
#' n_mc <- 10
#' n_z <- 50
#'
#' weighted_samps <- matrix(runif(n_draws * n_nodes), nrow = n_draws)
#' phi_array <- matrix(rexp(n_draws * n_nodes, 1), nrow = n_draws)
#' zi_array <- matrix(runif(n_draws * n_nodes, 0, 0.2), nrow = n_draws)
#' weights <- rep(1, n_nodes)
#' z_values <- seq(0, sum(weights), length.out = n_z)
#'
#' density <- ZIB_convolution_density(weighted_samps, phi_array, zi_array, weights, z_values, n_mc)
#' plot(z_values, density, type = "l")
#'
#' @export

ZIB_convolution_density <- function(z, alpha_matrix, beta_matrix,
                                    zi_matrix, weighted_samps, weights) {

  n_sims <- dim(weighted_samps)[1]
  n_draws <- dim(weighted_samps)[2]
  N <- dim(weighted_samps)[3]

  conv_pdf <- matrix(0, nrow = n_sims, ncol=n_draws)

  for(s in 1:n_sims) {
    for(m in 1:n_draws) {
      conv_pdf[s,m] <- dZIB_4p(x = z - sum(weighted_samps[s,m, 1:(N-1)]),
                               alpha_point = alpha_matrix[m,N],
                               beta_point = beta_matrix[m,N],
                               zi_point = zi_matrix[m,N],
                               weight = weights[N])
    }
  }
  avg_over_sims <- apply(conv_pdf, 2, mean)

  # Compute final result
  avg_over_draws <- mean(avg_over_sims, na.rm=TRUE)

  return(avg_over_draws)

}
