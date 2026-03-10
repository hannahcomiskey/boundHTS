#' Vectorized Zero-Inflated Four-Parameter Beta Density
#'
#' @param z evaluation points
#' @param Y_mc Monte Carlo draws of Y
#' @param phi_mc Monte Carlo draws of Beta precision (n_draws x n_nodes x n_years)
#' @param zi_mc Monte Carlo draws of zero inflation (n_draws x n_nodes x n_years)
#' @param lower lower bound
#' @param upper upper bound
#' @details
#' Computes the density of a zero-inflated four-parameter Beta distribution using Monte Carlo draws.
#' The function evaluates the density at points `z` while accounting for hierarchical aggregation across nodes.
#'
#' @return A numeric vector containing the estimated density values for each Monte Carlo draw.
#'
#' @examples
#' set.seed(1)
#'
#' # Monte Carlo size
#' n_mc <- 100
#' n_nodes <- 1
#'
#' # Simulated Monte Carlo draws
#' Y_mc  <- matrix(runif(n_mc * n_nodes, 0.2, 0.8), nrow = n_mc)
#' phi_mc <- matrix(rexp(n_mc * n_nodes, 1), nrow = n_mc)
#' zi_mc  <- matrix(runif(n_mc * n_nodes, 0, 0.2), nrow = n_mc)
#'
#' # Evaluation points
#' z <- seq(0, 1, length.out = 50)
#'
#' # Bounds
#' lower <- 0
#' upper <- 1
#'
#' dens <- dZIB_4p(z[25], Y_mc, phi_mc, zi_mc, upper, lower)
#' head(dens)
#'
#' @importFrom ExtDist dBeta_ab
#' @export

dZIB_4p <- function(z, Y_mc, phi_mc, zi_mc, upper, lower = 0) {

  n_nodes <- ncol(Y_mc)
  n_mc <- nrow(phi_mc)

  # Determine parent node and child sum
  if (n_nodes > 1) {
    parent   <- n_nodes
    child_sum <- rowSums(Y_mc[, -parent, drop = FALSE]) # all nodes except the last
    mu <- Y_mc[, parent]
    zi_vec <- zi_mc[, parent]
    upper_vec <- upper[parent]
  } else {
    child_sum <- 0 # no other nodes
    parent <- 1
    mu <- Y_mc[, 1]
    zi_vec <- zi_mc[, 1]
    upper_vec <- upper
  }

  # Remaining x value for parent
  x_vals <- z - child_sum

  # Beta parameters
  alpha <- pmax(mu * phi_mc[, parent], 1e-4)
  beta  <- pmax((1 - mu) * phi_mc[, parent], 1e-4)

  # Scale to [0,1] interval
  x_scaled <- x_vals / upper_vec

  # Initialize density to 0
  dens <- numeric(n_mc)

  # Boundary handling
  at0 <- x_scaled == 0 # handles values == 0
  inside <- x_scaled > 0 & x_scaled < 1 # handles values inside range
  outside <- x_scaled < 0 | x_scaled > 1 # handles values outside range

  dens[outside] <- 0
  dens[at0] <- zi_vec[at0]
  dens[inside] <- (1 - zi_vec[inside]) * ExtDist::dBeta_ab(x_vals[inside],
                                                           alpha[inside],
                                                           beta[inside],
                                                           lower,
                                                           upper_vec)

  return(dens)
}
