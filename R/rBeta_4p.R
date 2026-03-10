#' Vectorized 4 parameter Beta sampler using `ExtDist::rBeta_ab`.
#'
#' @param n_mc number of Monte Carlo samples
#' @param weights node weights (vector of length n_nodes)
#' @param sub_obs_data matrix of bottom-level series to use as expectation in sampler (n_years x n_nodes)
#' @param phi_array array (n_draws x n_nodes x n_years)
#' @return 3D array (n_draws x n_nodes x n_years)

rBeta_4p <- function(n_mc, sub_obs_data, phi_array, weights) {

  sub_obs_data <- as.matrix(sub_obs_data)

  n_years <- nrow(sub_obs_data)
  n_nodes <- ncol(sub_obs_data)
  n_draws <- dim(phi_array)[1]

  if(n_nodes != dim(phi_array)[2] | n_nodes != length(weights)) {
    message("Error: Dimension mismatch across inputs. Check data, phi and weights. ")
    stop()
  }

  # Monte Carlo particles
  Y <- array(NA_real_, dim = c(n_mc, n_nodes, n_years))

  # sample posterior draw index for each sample
  draw_id <- sample(seq_len(n_draws), n_mc, replace = TRUE)

  for (m in seq_len(n_mc)) {

    r <- draw_id[m]

    phi_r <- phi_array[r, , ]   # matrix [n_nodes × n_years]

    for (t in seq_len(n_years)) {

      mu <- sub_obs_data[t, ] # vector [n_nodes]

      # beta part
      alpha <- mu * phi_r[,t] # vector [n_nodes]
      beta  <- (1 - mu) * phi_r[,t] # vector [n_nodes]

      for(j in 1:n_nodes){
        Y[m, j, t] <- ExtDist::rBeta_ab(1, alpha[j], beta[j], 0, weights[j])
      }
    }
  }

  return(Y)
}
