#' Monte Carlo Estimate of Aggregated Zero-Inflated Four-Parameter Beta (ZIB) Density
#'
#' @param z Numeric evaluation point for the aggregated density.
#' @param weighted_samps Numeric matrix; Monte Carlo draws of weighted
#' bottom-series samples (n_draws x n_nodes).
#' @param alpha_point Point estimates of alpha (shape1) parameters for the Beta
#'   distribution for each observation.
#' @param beta_point Point estimates of beta (shape2) parameters for the Beta
#'   distribution for each observation.
#' @param zi_point Numeric vector; zero-inflation probability (length = n_nodes).
#' @param weights Numeric vector; node-specific upper bounds for the
#' Beta distribution (length = n_nodes).
#'
#' @details
#' For each evaluation point `z`, a Monte Carlo estimate of the aggregated
#' zero-inflated Beta density is computed from `alpha_point`, `beta_point`,
#' and `zi_point` and averaging `dZIB_4p` across the selected draws.
#' The resulting density is normalized using the trapezoidal rule
#' (`pracma::trapz`) to ensure integration to 1.
#'
#' @return Numeric vector of the same length as `z_values` representing
#' the normalized aggregated zero-inflated Beta density.
#'
#' @examples
#' set.seed(1)
#' n_draws <- 5
#' n_nodes <- 2
#' n_mc <- 10
#' weighted_samps <- array(runif(n_mc*n_draws* n_nodes), dim = c(n_mc, n_draws, n_nodes))
#' z_values <- seq(0, 1, out.length=100)
#' zi_point <- runif(n_nodes, 0, 0.2)
#' alpha_point <- runif(n_nodes, 2, 5)
#' beta_point  <- runif(n_nodes, 2, 5)
#' weights <- rep(1, n_nodes)
#' for(z in 1:length(z_values)) {
#'   denisty <- vector()
#'   density[z] <- ZIB_convolution_density_point(z_values[z],
#'   weighted_samps = weighted_samps, alpha_point[n_nodes], beta_point[n_nodes],
#'    zi_point[n_nodes], weights[n_nodes])
#'   }
#'   plot(z_values, density, type = "l")
#' @export

ZIB_convolution_density_point <- function(z, weighted_samps,
                                          alpha_point, beta_point, zi_point,
                                          weights) {


  # Dims
  N <- dim(weighted_samps)[3]

  partial_sum <- rowSums(weighted_samps[, , 1:(N-1), drop = FALSE], dims = 2)

  # Remaining x value for parent
  x <- z - partial_sum

  Density <- mean(dZIB_4p(x = x,
                          alpha_point = alpha_point,
                          beta_point = beta_point,
                          zi_point = zi_point,
                          weight = weights), na.rm = TRUE)

  return(Density)

}
