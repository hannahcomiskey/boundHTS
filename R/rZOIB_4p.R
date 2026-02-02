#' Vectorized Zero-One Inflated 4 parameter Beta sampler
#'
#' @param n_mc number of Monte Carlo samples
#' @param sub_obs_data matrix (n_years x n_nodes)
#' @param phi_array array (n_draws x n_nodes x n_years)
#' @param zoi_array array (n_draws x n_nodes x n_years)
#' @param coi_array array (n_draws x n_nodes x n_years)
#' @param weights node weights (vector of length n_nodes)
#' @return 3D array (n_draws x n_nodes x n_years)

rZOIB_4p <- function(n_mc, sub_obs_data, phi_array, zoi_array, coi_array, weights) {

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
    zoi_r <- zoi_array[r, , ]
    coi_r <- coi_array[r, , ]

    for (t in seq_len(n_years)) {

      mu <- sub_obs_data[t, ]

      inflate <- stats::runif(n_nodes) < zoi_r[t]

      # zero/one inflation
      if (any(inflate)) {
        Y[m, inflate, t] <-
          stats::rbinom(sum(inflate), 1, prob = coi_r[t]) *
          weights[inflate]
      }

      # beta part
      if (any(!inflate)) {
        idx <- which(!inflate)

        alpha <- mu[idx] * phi_r[t]
        beta  <- (1 - mu[idx]) * phi_r[t]

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
