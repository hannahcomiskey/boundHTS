#' Monte Carlo Estimate of Aggregated Zero-Inflated Four-Parameter Beta (ZIB) Density
#'
#' @param weighted_samps Numeric matrix; Monte Carlo draws of weighted
#' bottom-series samples (n_draws x n_nodes).
#' @param phi_point Numeric vector; Beta precision parameter (length = n_nodes).
#' @param zi_point Numeric vector; zero-inflation probability (length = n_nodes).
#' @param weights Numeric vector; node-specific upper bounds for the
#' Beta distribution (length = n_nodes).
#' @param z_values Numeric vector; evaluation points over which the
#' density is calculated.
#' @param n_mc Integer; number of Monte Carlo samples to use for aggregation.
#'
#' @details
#' For each evaluation point `z`, a Monte Carlo estimate of the aggregated
#' zero-inflated Beta density is computed by resampling posterior draws
#' from `phi_point` and `zi_point` and averaging `dZIB_4p` across the
#' selected draws. The resulting density is normalized using the
#' trapezoidal rule (`pracma::trapz`) to ensure integration to 1.
#'
#' @return Numeric vector of the same length as `z_values` representing
#' the normalized aggregated zero-inflated Beta density.
#'
#' @examples
#' set.seed(1)
#'
#' n_draws <- 5
#' n_nodes <- 2
#' n_mc <- 10
#' n_z <- 50
#'
#' weighted_samps <- matrix(runif(n_draws * n_nodes), nrow = n_draws)
#' phi_point <- rexp(n_nodes, 1)
#' zi_point <- runif(n_nodes, 0, 0.2)
#' weights <- rep(1, n_nodes)
#' z_values <- seq(0, sum(weights), length.out = n_z)
#'
#' density <- ZIB_convolution_density_point(weighted_samps, phi_point, zi_point, weights, z_values, n_mc)
#' plot(z_values, density, type = "l")
#'
#' @export

ZIB_convolution_density_point <- function(weighted_samps,
                                phi_point,
                                zi_point,
                                weights,
                                z_values,
                                n_mc) {

  n_draws <- dim(weighted_samps)[1]

  future::plan(future::sequential)

  Density <- future.apply::future_sapply(
    z_values,
    function(z) {
      library(boundHTS)
      mean(
        dZIB_4p(z = z,
                weighted_samps = weighted_samps,
                phi_mc = phi_point,
                zi_mc = zi_point,
                upper = weights), na.rm = TRUE)
    },
    future.seed = TRUE
  )

  norm_dens <- Density / pracma::trapz(z_values, Density)

  return(norm_dens)

}
