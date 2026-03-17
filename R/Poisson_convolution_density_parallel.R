#' Parallelized normalised predictive Poisson density over vector z
#'
#' @param z_values evaluation points
#' @param lambda_matrix Matrix of lambda parameters for the Poisson
#' distribution for each bottom level observation (rows correspond to simulations).
#' @details
#' A wrapper function for calculating Monte Carlo estimates for the aggregate density Z.
#'
#' @return The aggregate density Z over a grid of values.
#' @export

Poisson_convolution_density_parallel <- function(z_values, lambda_matrix) {
  Density <- future.apply::future_sapply(z_values,
                                         Poisson_convolution_density,
                                         lambda_matrix=lambda_matrix)
  return(Density / sum(Density)) # normalise
}
