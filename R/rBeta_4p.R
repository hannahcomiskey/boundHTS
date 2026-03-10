#' Vectorized Four-Parameter Beta Density Sampler
#'
#' @param n_mc Integer; number of Monte Carlo samples to generate.
#' @param sub_obs_data Numeric matrix of bottom-level series to use as the
#'   mean parameter in the sampler (dimensions: n_years x n_nodes).
#' @param phi_array Numeric array of posterior draws for the precision
#'   parameter (dimensions: n_draws x n_nodes x n_years).
#' @param weights Numeric vector of node-specific upper bounds (length = n_nodes).
#'
#' @details
#' For each Monte Carlo sample, a posterior draw is selected from `phi_array`.
#' The Beta shape parameters are computed as \code{alpha = mu * phi} and
#' \code{beta = (1 - mu) * phi}, where `mu` comes from `sub_obs_data`.
#' The resulting samples are drawn using \code{ExtDist::rBeta_ab} and returned
#' as a 3D array (n_mc x n_nodes x n_years).
#'
#' @return A numeric array of dimension \code{n_mc x n_nodes x n_years} containing the Monte Carlo samples.
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
#' weights <- rep(1, n_nodes)
#'
#' samples <- rBeta_4p(n_mc, sub_obs_data, phi_array, weights)
#' dim(samples)
#'
#' @export

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
