#' Monte Carlo Convolution Density for Zero-One-Inflated Beta Distributions
#'
#'
#' @param weighted_samps Monte Carlo draws of weighted bottom-series samples
#'   (matrix: n_draws x n_nodes).
#' @param phi_point Beta precision parameter (vector: length n_nodes).
#' @param zoi_point Zero-one inflation probability (vector: length n_nodes).
#' @param coi_point Conditional one inflation probability (vector: length n_nodes).
#' @param weights Node weights defining the aggregation structure (vector: length n_nodes).
#' @param z_values Numeric vector of evaluation points for the aggregate density.
#' @param n_mc Number of Monte Carlo samples used in the convolution estimate.
#'
#' @details
#' This function computes a Monte Carlo approximation of the density of an
#' aggregated Zero-One-Inflated Beta (ZOIB) distribution. Posterior draws are passed to \code{dZOIB_4p()},
#' which evaluates the ZOIB density for each Monte Carlo draw.
#' The resulting density is normalised using numerical integration.
#'
#' @return A numeric vector containing the estimated aggregate density evaluated at \code{z_values}.
#'
#' @examples
#' set.seed(1)
#'
#' # Simulation setup
#' n_draws <- 200
#' n_nodes <- 2
#'
#' # Simulated data
#' weighted_samps <- matrix(runif(n_draws * n_nodes, 0.1, 0.9), nrow = n_draws)
#' phi_point <- rexp(n_nodes, 1)
#' zoi_point <- runif(n_nodes, 0, 0.2)
#' coi_point <- runif(n_nodes, 0, 0.2)
#'
#' weights <- c(1, 1)
#' z_values <- seq(0, 2, length.out = 50)
#'
#' dens <- ZOIB_convolution_density_point(
#'   weighted_samps = weighted_samps,
#'   phi_point = phi_point,
#'   zoi_point = zoi_point,
#'   coi_point = coi_point,
#'   weights = weights,
#'   z_values = z_values,
#'   n_mc = 100
#' )
#'
#' head(dens)
#'
#' @export

ZOIB_convolution_density_point <- function(weighted_samps, phi_point, zoi_point,
                                           coi_point, weights, z_values, n_mc) {

  n_draws <- dim(weighted_samps)[1]

  future::plan(future::sequential)

  Density <- future.apply::future_sapply(
    z_values,
    function(z) {
      library(boundHTS)
      mean(
        dZOIB_4p(z = z,
                 Y_mc = weighted_samps,
                 phi_mc = phi_point,
                 zoi_mc = zoi_point,
                 coi_mc = coi_point,
                 upper = weights), na.rm = TRUE)
    },
    future.seed = TRUE
  )

  norm_dens <- Density / pracma::trapz(z_values, Density)

  return(norm_dens)

}
