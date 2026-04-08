#' Vectorized Zero-Inflated Four-Parameter Beta Density
#'
#' @param x evaluation points
#' @param alpha_point Point estimates of alpha (shape1) parameters for the Beta
#'   distribution for each observation.
#' @param beta_point Point estimates of beta (shape2) parameters for the Beta
#'   distribution for each observation.
#' @param zi_point Point estimates of zero inflation
#' @param weight Numeric vector of weights used to combine the components
#'   into the aggregated density \eqn{Z} (length b).
#' @details
#' Computes the density of a zero-inflated four-parameter Beta distribution using Monte Carlo draws.
#' The function evaluates the density at points `z` while accounting for hierarchical aggregation across nodes.
#'
#' @return A numeric vector containing the estimated density values for each Monte Carlo draw.
#'
#' @examples
#' set.seed(1)
#' n_draws <- 5
#' n_nodes <- 2
#' n_mc <- 10
#' weighted_samps <- array(runif(n_mc*n_draws* n_nodes), dim = c(n_mc, n_draws, n_nodes))
#' zi_point <- runif(n_nodes, 0, 0.2)
#' alpha_point <- runif(n_nodes, 2, 5)
#' beta_point  <- runif(n_nodes, 2, 5)
#' weights <- rep(1, n_nodes)
#' dZIB_4p(0.5, alpha_point[n_nodes], beta_point[n_nodes], zi_point[n_nodes],
#'  weights[n_nodes])
#'
#' @export
#' @importFrom ExtDist dBeta_ab
#' @export

dZIB_4p <- function(x, alpha_point, beta_point, zi_point, weight) {

  # Scale to [0,1] interval
  x_scaled <- x / weight

  # Initialize density to 0
  dens <- x
  dens[] <- 0

  # Boundary handling
  at0 <- x_scaled == 0 # handles values == 0
  inside <- x_scaled > 0 & x_scaled < 1 # handles values inside range
  outside <- x_scaled < 0 | x_scaled > 1 # handles values outside range

  dens[outside] <- 0
  dens[at0] <- zi_point
  dens[inside] <- (1 - zi_point) * ExtDist::dBeta_ab(x[inside],
                                                     alpha_point,
                                                     beta_point,
                                                     0, weight)

  return(dens)

}
