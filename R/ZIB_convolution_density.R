#' Monte Carlo Estimate of Aggregated Zero-Inflated Four-Parameter Beta (ZIB) Density
#'
#' @param Y_mc Numeric matrix; Monte Carlo draws of weighted bottom-series samples (n_draws x n_nodes).
#' @param phi_array Numeric matrix; Monte Carlo draws of the Beta precision parameter (n_draws x n_nodes).
#' @param zi_array Numeric matrix; Monte Carlo draws of zero-inflation probability (n_draws x n_nodes).
#' @param weights Numeric vector; node-specific upper bounds for the Beta distribution (length = n_nodes).
#' @param z_values Numeric vector; evaluation points over which the density is calculated.
#' @param n_mc Integer; number of Monte Carlo samples to use for aggregation.
#'
#' @details
#' For each evaluation point `z`, a Monte Carlo estimate of the aggregated zero-inflated Beta
#' density is computed by resampling posterior draws from `phi_array` and `zi_array` and
#' averaging `dZIB_4p` across the selected draws. The resulting density is normalized
#' using the trapezoidal rule (`pracma::trapz`) to ensure integration to 1.
#'
#' @return Numeric vector of the same length as `z_values` representing the normalized
#'   aggregated zero-inflated Beta density.
#'
#' @examples
#' set.seed(1)
#'
#' n_draws <- 5
#' n_nodes <- 2
#' n_mc <- 10
#' n_z <- 50
#'
#' Y_mc <- matrix(runif(n_draws * n_nodes), nrow = n_draws)
#' phi_array <- matrix(rexp(n_draws * n_nodes, 1), nrow = n_draws)
#' zi_array <- matrix(runif(n_draws * n_nodes, 0, 0.2), nrow = n_draws)
#' weights <- rep(1, n_nodes)
#' z_values <- seq(0, sum(weights), length.out = n_z)
#'
#' density <- ZIB_convolution_density(Y_mc, phi_array, zi_array, weights, z_values, n_mc)
#' plot(z_values, density, type = "l")
#'
#' @export

ZIB_convolution_density <- function(Y_mc,
                                phi_array,
                                zi_array,
                                weights,
                                z_values,
                                n_mc) {

  n_draws <- dim(Y_mc)[1]

  # Resample posterior draws ONCE
  draw_id <- sample(seq_len(n_draws), n_mc, replace = TRUE)

  phi_mc <- phi_array[draw_id, , drop=FALSE]
  zi_mc <- zi_array[draw_id, , drop=FALSE]

  future::plan(future::sequential)

  Density <- future.apply::future_sapply(
    z_values,
    function(z) {
      mean(dZIB_4p(z = z,
                   Y_mc = Y_mc,
                   phi_mc = phi_mc,
                   zi_mc = zi_mc,
                   upper = weights), na.rm = TRUE)
    },
    future.seed = TRUE
  )

  norm_dens <- Density / pracma::trapz(z_values, Density)

  return(norm_dens)

}
