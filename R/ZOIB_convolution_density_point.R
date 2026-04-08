#' Monte Carlo Convolution Density for Zero-One-Inflated Beta Distributions
#'
#'
#' @param weighted_samps Monte Carlo draws of weighted bottom-series samples
#'   (matrix: n_draws x n_nodes).
#' @param alpha_point Point estimates of alpha (shape1) parameters for the Beta
#'   distribution for each observation.
#' @param beta_point Point estimates of beta (shape2) parameters for the Beta
#'   distribution for each observation.
#' @param zoi_point Zero-one inflation probability (vector: length n_nodes).
#' @param coi_point Conditional one inflation probability (vector: length n_nodes).
#' @param weights Node weights defining the aggregation structure (vector: length n_nodes).
#' @param z_values Numeric vector of evaluation points for the aggregate density.
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
#' alpha_point <- runif(n_nodes, 2, 5)
#' beta_point  <- runif(n_nodes, 2, 5)
#' zoi_point <- runif(n_nodes, 0, 0.2)
#' coi_point <- runif(n_nodes, 0, 0.2)
#'
#' weights <- c(1, 1)
#' z_values <- seq(0, 2, length.out = 50)
#'
#' dens <- ZOIB_convolution_density_point(
#'   weighted_samps = weighted_samps,
#'   alpha_point = alpha_point,
#'   beta_point = beta_point,
#'   zoi_point = zoi_point,
#'   coi_point = coi_point,
#'   weights = weights,
#'   z_values = z_values
#' )
#'
#' head(dens)
#'
#' @export

ZOIB_convolution_density_point <- function(weighted_samps, alpha_point, beta_point,
                                           zoi_point, coi_point, weights, z_values) {


  # Dims
  N <- dim(weighted_samps)[3]

  partial_sum <- rowSums(weighted_samps[, , 1:(N-1), drop = FALSE], dims = 2)

  # Remaining x value for parent
  x <- z_values - partial_sum

  Density <- mean(dZOIB_4p(x = x,
                          alpha_point = alpha_point[N],
                          beta_point = beta_point[N],
                          zi_point = zoi_point[N],
                          coi_point = coi_point[N],
                          weight = weights[N]), na.rm = TRUE)

  return(Density)

}
