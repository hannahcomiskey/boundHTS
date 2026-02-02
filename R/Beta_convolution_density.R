#' Calculate 4 parameter Beta PDF
#'
#' @param z evaluation points
#' @param alpha_matrix matrix of shape parameters for each element (b) in the aggregate over each of the N observations (N rows x b columns)
#' @param beta_matrix matrix of shape parameters for each element in the aggregate (N rows x b columns)
#' @param weighted_samps matrix of weighted samples for each element in the aggregate (N rows x b columns)
#' @param weights vector of weights that are used to combine to form the aggregated density Z. (length b)
#' @details
#' A wrapper function for calculating Monte Carlo estimates for the aggregated Beta density Z.
#'
#' @return The aggregate density Z over a grid of values.
#' @export

Beta_convolution_density <- function(z, alpha_matrix, beta_matrix, weighted_samps, weights) {
  N <- dim(weighted_samps)[3]
  n_sims <- dim(weighted_samps)[1]
  n_draws <- dim(weighted_samps)[2]
  conv_pdf <- matrix(NA, nrow = n_sims, ncol=n_draws)

  for(m in 1:n_draws) {
    for(s in 1:n_sims) {
      conv_pdf[s,m] <- dBeta_4p(x = z - sum(weighted_samps[s,m, 1:(N-1)]),
                                   upper = weights[N],
                                   alpha = alpha_matrix[s,N],
                                   beta = beta_matrix[s,N])
    }
  }
  avg_over_sims <- apply(conv_pdf, 2, mean)

  # Compute final result
  avg_over_draws <- mean(avg_over_sims, na.rm=TRUE)

  return(avg_over_draws)
}
