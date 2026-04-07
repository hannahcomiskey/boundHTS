#' Parallelized normalised predictive density over vector z
#'
#' @param z_values evaluation points
#' @param alpha_input Point estimates (point=TRUE) or matrix (point=FALSE) of
#' shape parameters for each element (b) in the aggregate over each of the N observations (N rows x b columns)
#' @param beta_input Point estimates (point=TRUE) or matrix (point=FALSE) of shape parameters for each element in the aggregate (N rows x b columns)
#' @param weighted_samps Array of weighted samples for each element in the
#'   aggregate (dimensions: n_sims x n_draws x b).
#' @param weights vector of weights that are used to combine to form the aggregated density Z. (length b)
#' @param point A true/false indicator to denote whether you are using
#' point estimates (point=TRUE) or posterior samples (point=FALSE) of the Beta parameters.
#' @details
#' A wrapper function for calculating Monte Carlo estimates for the aggregate density Z using the Beta distribution.
#' @return The aggregate density Z over a grid of values using a convolution of Beta distributions.
#' @export

Beta_convolution <- function(z_values, alpha_input, beta_input, weighted_samps, weights) {
  if(point==TRUE & is.vector(alpha_input)==TRUE & is.vector(beta_input)==TRUE) {
    dens <- Beta_convolution_density_point_parallel(z_values = z_values,
                                                    alpha_point = alpha_input,
                                                    beta_point = beta_input,
                                                    weighted_samps = weighted_samps,
                                                    weights = weights)
  }
  if(point==FALSE & is.matrix(alpha_input)==TRUE & is.matrix(beta_input)==TRUE) {
    dens <- Beta_convolution_density_parallel(z_values = z_values,
                                              alpha_matrix = alpha_input,
                                              beta_matrix = beta_input,
                                              weighted_samps = weighted_samps,
                                              weights = weights)
  }
  return(dens)
}
