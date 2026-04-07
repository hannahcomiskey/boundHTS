#' Parallelized normalised predictive density over vector z
#'
#' @param z_values evaluation points
#' @param alpha_input Point estimates (point=TRUE) or matrix (point=FALSE) of
#' shape parameters for each element (b) in the aggregate over each of the N observations.
#' @param beta_input Point estimates (point=TRUE) or matrix (point=FALSE) of shape parameters for each element in the aggregate.
#' @param zi_input Point estimates (point=TRUE) or matrix (point=FALSE) of zero-inflation probability parameters for each element in the aggregate.
#' @param weighted_samps Array of weighted samples for each element in the
#'   aggregate (dimensions: n_sims x n_draws x b).
#' @param weights vector of weights that are used to combine to form the aggregated density Z. (length b)
#' @param point A true/false indicator to denote whether you are using
#' point estimates (point=TRUE) or posterior samples (point=FALSE) of the Beta parameters.
#' @details
#' A wrapper function for calculating Monte Carlo estimates for the aggregate density Z using the zero-inflated Beta distribution.
#' @return The aggregate density Z over a grid of values using a convolution of zero-inflated Beta distributions.
#' @export

ZIB_convolution <- function(z_values, alpha_input, beta_input, zi_input, weighted_samps, weights, point) {
  if(point==TRUE & is.vector(alpha_input)==TRUE & is.vector(beta_input)==TRUE & is.vector(zi_input)==TRUE) {
    dens <- ZIB_convolution_density_point_parallel(z_values = z_values,
                                                   alpha_point = alpha_input,
                                                   beta_point = beta_input,
                                                   zi_point = zi_input,
                                                   weighted_samps = weighted_samps,
                                                   weights = weights)
  }
  if(point==FALSE & is.matrix(alpha_input)==TRUE & is.matrix(beta_input)==TRUE & is.matrix(zi_input)==TRUE) {
    dens <- ZIB_convolution_density_parallel(z_values = z_values,
                                             alpha_matrix = alpha_input,
                                             beta_matrix = beta_input,
                                             zi_matrix = zi_input,
                                             weighted_samps = weighted_samps,
                                             weights = weights)
  }
  return(dens)
}
