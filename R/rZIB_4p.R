#' Vectorized Zero-Inflated Four-Parameter Beta Density Sampler
#'
#' @param n_mc Integer; number of Monte Carlo samples to generate.
#' @param sub_obs_data Numeric matrix of observed bottom-level series to use as
#'   the mean parameter in the sampler (dimensions: n_years x n_nodes).
#' @param phi_array Numeric array of posterior draws for the Beta precision
#'   parameter (dimensions: n_draws x n_nodes x n_years).
#' @param zi_array Numeric array of posterior draws for the zero-inflation
#'   probability (dimensions: n_draws x n_nodes x n_years).
#' @param weights Numeric vector of node-specific upper bounds for the
#'   Beta distribution (length = n_nodes).
#'
#' @details
#' For each Monte Carlo sample, a posterior draw is selected from `phi_array` and `zi_array`.
#' Zero inflation is applied according to the sampled `zi_r` probabilities.
#' For nodes not set to zero by the inflation, samples are drawn from a Beta distribution
#' scaled to `[0, weight]` using `ExtDist::rBeta_ab`.
#'
#' @return A numeric array of dimension \code{n_mc x n_nodes x n_years} containing
#'   the Monte Carlo samples of the zero-inflated Beta distribution.
#'
#' @examples
#' set.seed(1)
#'
#' n_mc <- 10
#' n_years <- 3
#' n_nodes <- 2
#' n_draws <- 5
#'
#' sub_obs_data <- matrix(runif(n_years * n_nodes, 0.2, 0.8),
#'                        nrow = n_years, ncol = n_nodes)
#' phi_array <- array(rexp(n_draws * n_nodes * n_years, 1),
#'                    dim = c(n_draws, n_nodes, n_years))
#' zi_array <- array(runif(n_draws * n_nodes * n_years, 0, 0.2),
#'                   dim = c(n_draws, n_nodes, n_years))
#' weights <- rep(1, n_nodes)
#'
#' samples <- rZIB_4p(n_mc, sub_obs_data, phi_array, zi_array, weights)
#' dim(samples)
#'
#' @export

rZIB_4p <- function(n_mc, sub_obs_data, phi_array, zi_array, weights) {

  sub_obs_data <- as.matrix(sub_obs_data)

  n_years <- nrow(sub_obs_data)
  n_nodes <- ncol(sub_obs_data)
  n_draws <- dim(phi_array)[1]

  # Monte Carlo particles
  Y <- array(NA_real_, dim = c(n_mc, n_nodes, n_years))

  # sample posterior draw index for each particle
  draw_id <- sample(seq_len(n_draws), n_mc, replace = TRUE)

  for (m in seq_len(n_mc)) {

    r <- draw_id[m]

    phi_r <- phi_array[r, , ]   # [nodes × years]
    zi_r <- zi_array[r, , ]

    for (t in seq_len(n_years)) {

      mu <- sub_obs_data[t, ]

      inflate <- stats::runif(n_nodes) < zi_r[,t]

      # zero inflation
      if (any(inflate)) {
        Y[m, inflate, t] <- 0
      }

      # beta part
      if (any(!inflate)) {
        idx <- which(!inflate)

        alpha <- mu[idx] * phi_r[idx, t]
        beta  <- (1 - mu[idx]) * phi_r[idx, t]

        Y[m, idx, t] <-
          ExtDist::rBeta_ab(
            length(idx),
            alpha,
            beta,
            0,
            weights[idx]
          )
      }
    }
  }

  return(Y)
}
