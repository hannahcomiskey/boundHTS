#' Monte Carlo Convolution Density for Zero-One-Inflated Beta Distributions
#'
#'
#' @param Y_mc Monte Carlo draws of weighted bottom-series samples
#'   (matrix: n_draws x n_nodes).
#' @param phi_array Monte Carlo draws of the Beta precision parameter
#'   (matrix: n_draws x n_nodes).
#' @param zoi_array Monte Carlo draws of the zero-one inflation probability
#'   (matrix: n_draws x n_nodes).
#' @param coi_array Monte Carlo draws of the conditional one inflation probability
#'   (matrix: n_draws x n_nodes).
#' @param weights Node weights (vector of length n_nodes) defining the
#'   aggregation structure.
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
#' # Simulated posterior draws
#' Y_mc <- matrix(runif(n_draws * n_nodes, 0.1, 0.9), nrow = n_draws)
#' phi_array <- matrix(rexp(n_draws * n_nodes, 1), nrow = n_draws)
#' zoi_array <- matrix(runif(n_draws * n_nodes, 0, 0.2), nrow = n_draws)
#' coi_array <- matrix(runif(n_draws * n_nodes, 0, 0.2), nrow = n_draws)
#'
#' weights <- c(1, 1)
#' z_values <- seq(0, 2, length.out = 50)
#'
#' dens <- ZOIB_convolution_density(
#'   Y_mc = Y_mc,
#'   phi_array = phi_array,
#'   zoi_array = zoi_array,
#'   coi_array = coi_array,
#'   weights = weights,
#'   z_values = z_values,
#'   n_mc = 100
#' )
#'
#' head(dens)
#'
#' @export

ZOIB_convolution_density <- function(Y_mc,
                                phi_array,
                                zoi_array,
                                coi_array,
                                weights,
                                z_values,
                                n_mc) {

  n_draws <- dim(Y_mc)[1]

  # Resample posterior draws ONCE
  draw_id <- sample(seq_len(n_draws), n_mc, replace = TRUE)

  phi_mc <- phi_array[draw_id, , drop=FALSE]
  zoi_mc <- zoi_array[draw_id, , drop=FALSE]
  coi_mc <- coi_array[draw_id, , drop=FALSE]

  future::plan(future::sequential)

  Density <- future.apply::future_sapply(
    z_values,
    function(z) {
      library(boundHTS)
      mean(
        dZOIB_4p(z = z,
                 Y_mc = Y_mc,
                 phi_mc = phi_mc,
                 zoi_mc = zoi_mc,
                 coi_mc = coi_mc,
                 upper = weights), na.rm = TRUE)
    },
    future.seed = TRUE
  )

  norm_dens <- Density / pracma::trapz(z_values, Density)

  return(norm_dens)

}
