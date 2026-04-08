#' Parallelized normalised predictive Poisson density over vector z
#'
#' @param z_values evaluation points
#' @param lambda_vector Point estimates of lambda parameters for the Poisson
#'   distribution for each bottom series observation.
#' @details
#' A wrapper function for calculating Monte Carlo estimates for the aggregate density Z using point estimates.
#'
#' @return The aggregate density Z over a grid of values.
#' @export

Poisson_convolution_density_point_parallel <- function(z_values, lambda_vector) {
  Density <- future.apply::future_sapply(z_values,
                                         Poisson_convolution_density_point,
                                         lambda_vector = lambda_vector)
  return(Density / sum(Density)) # normalised
}
