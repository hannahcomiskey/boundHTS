#' Parallelized normalised predictive zero-inflated Beta density over vector z
#'
#' @param z_values evaluation points
#' @param alpha_matrix matrix of shape parameters for each element (b) in the aggregate over each of the N observations (N rows x b columns)
#' @param beta_matrix matrix of shape parameters for each element in the aggregate (N rows x b columns)
#' @param zi_matrix matrix of zero-inflation probability parameters for each element in the aggregate (N rows x b columns)
#' @param weighted_samps matrix of weighted samples for each element in the aggregate (N rows x b columns)
#' @param weights vector of weights that are used to combine to form the aggregated density Z. (length b)
#' @details
#' A wrapper function for calculating Monte Carlo estimates for the aggregate density Z using point estimates.
#'
#' @return The aggregate density Z over a grid of values.
#' @export

ZIB_convolution_density_parallel <- function(z_values, alpha_matrix, beta_matrix, zi_matrix,
                                             weighted_samps, weights) {
  Density <- future.apply::future_sapply(z_values,
                                         ZIB_convolution_density,
                                         alpha_matrix=alpha_matrix,
                                         beta_matrix=beta_matrix,
                                         zi_matrix = zi_matrix,
                                         weighted_samps=weighted_samps,
                                         weights=weights)
  return(Density / pracma::trapz(z_values, Density))
}
