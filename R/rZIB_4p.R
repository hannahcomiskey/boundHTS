#' Vectorized Zero Inflated 4 parameter Beta sampler
#'
#' @param n_mc number of Monte Carlo samples
#' @param sub_obs_data matrix of observed data (n_years x n_nodes)
#' @param phi_array Monte Carlo draws of Beta precision (n_draws x n_nodes x n_years)
#' @param zi_array Monte Carlo draws of zero inflation (n_draws x n_nodes x n_years)
#' @param weights node weights for 4 parameter Beta (vector of length n_nodes)
#' @return 3D array (n_draws x n_nodes x n_years)

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
