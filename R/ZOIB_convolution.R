#' Parallelized normalised predictive density over vector z
#'
#' @param z_values evaluation points
#' @param alpha_input Point estimates (point=TRUE) or matrix (point=FALSE) of
#' shape parameters for each element (b) in the aggregate over each of the N observations.
#' @param beta_input Point estimates (point=TRUE) or matrix (point=FALSE) of shape parameters for each element in the aggregate.
#' @param zoi_input Point estimates (point=TRUE) or matrix (point=FALSE) of zero-inflation probability parameters for each element in the aggregate.
#' @param coi_input Point estimates (point=TRUE) or matrix (point=FALSE) of conditional probability of getting exactly 1 for each element in the aggregate.
#' @param weighted_samps Array of weighted samples for each element in the
#'   aggregate (dimensions: n_sims x n_draws x b).
#' @param weights vector of weights that are used to combine to form the aggregated density Z. (length b)
#' @param point A true/false indicator to denote whether you are using
#' point estimates (point=TRUE) or posterior samples (point=FALSE) of the Beta parameters.
#' @details
#' A wrapper function for calculating Monte Carlo estimates for the aggregate density Z using the zero-one inflated Beta (ZOIB) distribution.
#' @return The aggregate density Z over a grid of values using a convolution of zero-one inflated Beta (ZOIB) distributions.
#' @export

ZOIB_convolution <- function(z_values, alpha_input, beta_input, zi_input, weighted_samps, weights, point) {
  if(point==TRUE & is.vector(alpha_input)==TRUE & is.vector(beta_input)==TRUE & is.vector(zoi_input)==TRUE & is.vector(coi_input)==TRUE) {
    dens <- ZOIB_convolution_density_point_parallel(z_values = z_values,
                                                   alpha_point = alpha_input,
                                                   beta_point = beta_input,
                                                   zoi_point = zi_input,
                                                   coi_point = coi_input,
                                                   weighted_samps = weighted_samps,
                                                   weights = weights)
  }
  if(point==FALSE & is.matrix(alpha_input)==TRUE & is.matrix(beta_input)==TRUE & is.matrix(zoi_input)==TRUE  & is.matrix(coi_input)==TRUE) {
    dens <- ZOIB_convolution_density_parallel(z_values = z_values,
                                             alpha_matrix = alpha_input,
                                             beta_matrix = beta_input,
                                             zoi_matrix = zoi_input,
                                             coi_matrix = coi_input,
                                             weighted_samps = weighted_samps,
                                             weights = weights)
  }
  return(dens)
}
