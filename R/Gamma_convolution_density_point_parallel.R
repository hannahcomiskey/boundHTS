#' Parallelized normalised predictive density over vector z
#'
#' @param z_values evaluation points
#' @param shape_point Point estimates of shape parameters for the Gamma distribution for each observation.
#' @param rate_point Point estimates of rate parameters for the Gamma distribution for each observation.
#' @param bottom_samps Array of samples for each element in the aggregate (dimensions: n_sims x n_draws x b).
#' @details
#' A wrapper function for calculating Monte Carlo estimates for the aggregate density Z using point estimates.
#'
#' @return The aggregate density Z over a grid of values using a convolution of Gamma distributions.
#' @export

Gamma_convolution_density_point_parallel <- function(z_values, shape_point, rate_point, bottom_samps) {
  Density <- future.apply::future_sapply(z_values,
                                         Gamma_convolution_density_point,
                                         shape_point = shape_point,
                                         rate_point = rate_point,
                                         bottom_samps = bottom_samps)
  return(Density / pracma::trapz(z_values, Density))
}
