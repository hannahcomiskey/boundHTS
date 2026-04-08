#' Vectorized Zero-One-Inflated Four-Parameter Beta Density
#'
#' @param x evaluation points
#' @param alpha_point Point estimates of alpha (shape1) parameters for the Beta
#'   distribution for each observation.
#' @param beta_point Point estimates of beta (shape2) parameters for the Beta
#'   distribution for each observation.
#' @param zi_point Point estimates of zero inflation (n_nodes).
#' @param coi_point Point estimate of conditional one inflation probability (n_nodes).
#' @param weight Weight used to combine the components into the aggregated density \eqn{Z}.
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
#' weighted_samps  <- array(runif(n_mc * n_draws * n_nodes, 0.2, 0.8), dim=c(n_mc, n_draws, n_nodes))
#' alpha_point <- runif(n_nodes, 2, 5)
#' beta_point  <- runif(n_nodes, 2, 5)
#' zi_point <- runif(n_nodes, 0, 0.05)
#' coi_point <- runif(n_nodes, 0, 0.02)
#' # Evaluation points
#' z <- seq(0, 1, length.out = 50)
#' weights <- 0.5
#' dens <- dZOIB_4p(z[25], alpha_point[2], beta_point[2], zi_point[2], coi_point[2], weights)
#' dens
#' @importFrom ExtDist dBeta_ab
#' @export

dZOIB_4p <- function(x, alpha_point, beta_point, zi_point, coi_point, weight) {
  # Scale to [0,1] interval
  x_scaled <- x / weight

  if(length(x) > 1) {
    # Initialize density to 0
    dens <- vector()

    # Boundary handling
    at0 <- x_scaled == 0 # handles values == 0
    at1 <- x_scaled == 1 # handles values == 1
    inside <- x_scaled > 0 & x_scaled < 1 # handles values inside range
    outside <- x_scaled < 0 | x_scaled > 1 # handles values outside range

    dens[outside] <- 0
    dens[at0] <- zi_point * (1-coi_point)
    dens[at1] <- zi_point * coi_point
    dens[inside] <- (1 - zi_point) * ExtDist::dBeta_ab(x[inside],
                                                       alpha_point,
                                                       beta_point,
                                                       0, weight)
    } else {
        # Outside support
        if (x_scaled < 0 || x_scaled > 1) {
          return(0)
        }

        # At boundaries
        if (x_scaled == 0) {
          return(zi_point * (1 - coi_point))
        }

        if (x_scaled == 1) {
          return(zi_point * coi_point)
        }

        # Inside (0,1)
        dens <- (1 - zi_point) * ExtDist::dBeta_ab(x, alpha_point, beta_point, 0, weight)
      }
  return(dens)
}
