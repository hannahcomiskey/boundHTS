#' Parallelized normalised predictive zero-inflated Beta density over vector z
#'
#' @param z_values evaluation points
#' @param alpha_point Point estimates of alpha (shape1) parameters for the Beta
#'   distribution for each observation.
#' @param beta_point Point estimates of beta (shape2) parameters for the Beta
#'   distribution for each observation.
#' @param zi_point Numeric vector; zero-inflation probability (length = n_nodes).
#' @param weighted_samps Array of weighted samples for each element in the
#'   aggregate (dimensions: n_sims x n_draws x b).
#' @param weights Numeric vector of weights used to combine the components
#'   into the aggregated density \eqn{Z} (length b).
#' @details
#' A wrapper function for calculating Monte Carlo estimates for the aggregate density Z using point estimates.
#'
#' @return The aggregate density Z over a grid of values.
#' @export

ZIB_convolution_density_point_parallel <- function(z_values, alpha_point,
                                                    beta_point, zi_point,
                                                    weighted_samps, weights) {
  Density <- future.apply::future_sapply(z_values,
                                         ZIB_convolution_density_point,
                                         alpha_point = alpha_point,
                                         beta_point = beta_point,
                                         zi_point = zi_point,
                                         weighted_samps = weighted_samps,
                                         weights = weights)
  return(Density / pracma::trapz(z_values, Density))
}
