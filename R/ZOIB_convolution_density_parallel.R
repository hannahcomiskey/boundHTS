#' Parallelized normalised predictive zero-one inflated Beta density over vector z
#'
#' @param z_values evaluation points
#' @param alpha_matrix matrix of shape parameters for each element (b) in the aggregate over each of the N observations (N rows x b columns)
#' @param beta_matrix matrix of shape parameters for each element in the aggregate (N rows x b columns)
#' @param zoi_matrix matrix of zero-inflation probability parameters for each element in the aggregate (N rows x b columns)
#' @param coi_matrix Numeric matrix; Monte Carlo draws of conditional one-inflation probability (n_draws x n_nodes).
#' @param weighted_samps matrix of weighted samples for each element in the aggregate (N rows x b columns)
#' @param weights vector of weights that are used to combine to form the aggregated density Z. (length b)
#' @details
#' A wrapper function for calculating Monte Carlo estimates for the aggregate density Z using point estimates.
#'
#' @return The aggregate density Z over a grid of values.
#' @export

ZOIB_convolution_density_parallel <- function(z_values, alpha_matrix, beta_matrix,
                                              coi_matrix, zoi_matrix, weighted_samps, weights) {
  Density <- future.apply::future_sapply(z_values,
                                         ZOIB_convolution_density,
                                         alpha_matrix=alpha_matrix,
                                         beta_matrix=beta_matrix,
                                         zoi_matrix = zoi_matrix,
                                         coi_matrix = coi_matrix,
                                         weighted_samps=weighted_samps,
                                         weights=weights)
  return(Density / pracma::trapz(z_values, Density))
}
