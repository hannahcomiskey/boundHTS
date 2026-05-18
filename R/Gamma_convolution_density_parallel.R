#' Parallelized normalised predictive density over vector z
#'
#' @param z_values evaluation points
#' @param shape_matrix matrix of shape parameters for each element (b) in the aggregate over each of the N observations (N rows x b columns)
#' @param rate_matrix matrix of rate parameters for each element in the aggregate (N rows x b columns)
#' @param bottom_samps matrix of samples for each element in the aggregate (N rows x b columns)
#' @details
#' A wrapper function for calculating Monte Carlo estimates for the aggregate density Z.
#'
#' @return The aggregate density Z over a grid of values using a convolution of Gamma distributions.
#' @export

Gamma_convolution_density_parallel <- function(z_values, shape_matrix, rate_matrix, bottom_samps) {
  Density <- future.apply::future_sapply(z_values,
                                         Gamma_convolution_density,
                                         shape_matrix=shape_matrix,
                                         rate_matrix=rate_matrix,
                                         bottom_samps=bottom_samps)

  if(sum(Density)>0) {
    norm_dens <- Density / pracma::trapz(z_values, Density) # normalise density
  } else {
    norm_dens <- Density # leave if Density is 0 everywhere (out of bounds)
  }
  return(norm_dens)
}
