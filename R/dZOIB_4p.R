#' Vectorized Zero-One-Inflated Four-Parameter Beta Density
#'
#' @param z Numeric evaluation point for the density.
#' @param Y_mc Matrix of Monte Carlo draws of the mean parameter \eqn{\mu}
#'   (n_mc x n_nodes).
#' @param phi_mc Matrix of Monte Carlo draws of the precision parameter
#'   \eqn{\phi} (n_mc x n_nodes).
#' @param zoi_mc Matrix of Monte Carlo draws of the zero-one inflation
#'   probability (n_mc x n_nodes).
#' @param coi_mc Matrix of Monte Carlo draws of the conditional one inflation
#'   probability (n_mc x n_nodes).
#' @param upper Numeric vector giving the upper bound for each node.
#' @param lower Numeric scalar giving the lower bound (default is 0).
#'
#' @details
#' The Zero-One-Inflated Beta (ZOIB) distribution places probability mass at both
#' boundaries (0 and 1) and follows a Beta distribution on the interior of the
#' interval. This function evaluates the density conditional on Monte Carlo
#' draws of the model parameters and supports hierarchical aggregation by
#' subtracting the contribution of child nodes before evaluating the density
#' of the parent node.
#'
#' The Beta density is evaluated using \code{ExtDist::dBeta_ab()} for the
#' interior of the support.
#'
#' @return A numeric vector of length \code{n_mc} containing the density
#' evaluated at \code{z} for each Monte Carlo draw.
#'
#' @examples
#' set.seed(1)
#'
#' n_mc <- 100
#' n_nodes <- 2
#'
#' Y_mc <- matrix(runif(n_mc * n_nodes, 0.2, 0.8), nrow = n_mc)
#' phi_mc <- matrix(rexp(n_mc * n_nodes, 1), nrow = n_mc)
#' zoi_mc <- matrix(runif(n_mc * n_nodes, 0, 0.2), nrow = n_mc)
#' coi_mc <- matrix(runif(n_mc * n_nodes, 0, 0.2), nrow = n_mc)
#'
#' upper <- c(1, 1)
#'
#' dZOIB_4p(
#'   z = 0.5,
#'   Y_mc = Y_mc,
#'   phi_mc = phi_mc,
#'   zoi_mc = zoi_mc,
#'   coi_mc = coi_mc,
#'   upper = upper
#' )
#'
#' @export

dZOIB_4p <- function(z, Y_mc, phi_mc, zoi_mc, coi_mc, upper, lower = 0) {

  n_nodes <- ncol(Y_mc)
  n_mc <- nrow(phi_mc)

  # Determine parent node and child sum
  if (n_nodes > 1) {
    parent   <- n_nodes
    child_sum <- rowSums(Y_mc[, -parent, drop = FALSE]) # all nodes except the last
    mu <- Y_mc[, parent]
    zoi_vec <- zoi_mc[, parent]
    coi_vec <- coi_mc[, parent]
    upper_vec <- upper[parent]
  } else {
    child_sum <- 0 # no other nodes
    parent <- 1
    mu <- Y_mc[, 1]
    zoi_vec <- zoi_mc[, 1]
    coi_vec <- coi_mc[, 1]
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
  at1 <- x_scaled == 1 # handles values == 1
  inside <- x_scaled > 0 & x_scaled < 1 # handles values inside range
  outside <- x_scaled < 0 | x_scaled > 1 # handles values outside range

  dens[outside] <- 0
  dens[at0] <- zoi_vec[at0] * (1 - coi_vec[at0])
  dens[at1] <- zoi_vec[at1] * coi_vec[at1]
  dens[inside] <- (1 - zoi_vec[inside]) *
    ExtDist::dBeta_ab(x_vals[inside],
                      alpha[inside],
                      beta[inside],
                      lower,
                      upper_vec)

  return(dens)
}
