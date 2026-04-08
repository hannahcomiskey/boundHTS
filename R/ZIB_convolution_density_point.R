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
#' # Simulation setup
#' n_sims <- 50
#' n_draws <- 10
#' b <- 2
#'
#' # Simulated weighted samples
#' weighted_samps <- array(runif(n_sims * n_draws * b),
#'                         dim = c(n_sims, n_draws, b))
#'
#' alpha_point <- runif(b, 2, 5)
#' beta_point  <- runif(b, 2, 5)
#' zi_point  <- runif(b, 0, 0.1)
#' weights <- c(1, 1)
#'
#' density <- ZIB_convolution_density_point(z = 0.5, alpha_point = alpha_point,
#' beta_point = beta_point, zi_point = zi_point,
#' weighted_samps = weighted_samps, weights = weights)
#'
#' density
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
