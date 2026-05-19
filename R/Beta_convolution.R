#' Parallelized normalised predictive density over vector z
#'
#' @param S Binary summation matrix with rows corresponding to all
#' hierarchy nodes and columns corresponding to bottom-level nodes.
#'
#' @param alpha_mat Numeric vector (point = TRUE) or matrix
#' (point = FALSE) of Beta shape1/alpha parameters for all hierarchy nodes.
#'
#' @param beta_mat Numeric vector (point = TRUE) or matrix
#' (point = FALSE) of Beta shape2/beta parameters for all hierarchy nodes.
#'
#' @param weights_list A named list of local aggregation weights used
#' to construct the aggregation weight matrix.
#'
#' @param point Logical indicating whether point estimates
#' (`point = TRUE`) or posterior samples (`point = FALSE`)
#' are supplied.
#'
#' @param n_draws Number of posterior draws.
#'
#' @param n_sims Number of Monte Carlo simulations per draw.
#'
#' @param z_values Density evaluation grid.
#'
#' @details
#' A wrapper function for calculating Monte Carlo estimates for the
#' aggregate density Z using four-parameter Beta distributions.
#'
#' Aggregation is performed coherently using:
#'
#' \deqn{
#' y = Wb
#' }
#'
#' where:
#'
#' \itemize{
#' \item \eqn{b} are bottom-level random variables,
#' \item \eqn{W} is the aggregation weight matrix.
#' }
#'
#' Density estimation is only performed for aggregated hierarchy nodes.
#'
#' @return
#' A list of aggregate densities over a grid of values for each
#' aggregated hierarchy node.
#'
#' @examples
#' ## ---------------------------------------------------------
#' ## Example: Point-estimate Beta convolution
#' ## ---------------------------------------------------------
#'
#' set.seed(123)
#'
#' S <- rbind(
#'   Total = c(1, 1, 1, 1),
#'   A     = c(1, 1, 0, 0),
#'   B     = c(0, 0, 1, 1),
#'   A1    = c(1, 0, 0, 0),
#'   A2    = c(0, 1, 0, 0),
#'   B1    = c(0, 0, 1, 0),
#'   B2    = c(0, 0, 0, 1)
#' )
#'
#' local_weights <- list(
#'   Total = c(A1 = 0.10, A2 = 0.15, B1 = 0.30, B2 = 0.45),
#'   A     = c(A1 = 0.40, A2 = 0.60),
#'   B     = c(B1 = 0.30, B2 = 0.70)
#' )
#'
#' alpha_mat <- c(
#'   Total = 8,
#'   A = 5,
#'   B = 6,
#'   A1 = 2,
#'   A2 = 4,
#'   B1 = 3,
#'   B2 = 5
#' )
#'
#' beta_mat <- c(
#'   Total = 4,
#'   A = 3,
#'   B = 2,
#'   A1 = 6,
#'   A2 = 5,
#'   B1 = 7,
#'   B2 = 4
#' )
#'
#' dens <- Beta_convolution(
#'   S = S,
#'   alpha_mat = alpha_mat,
#'   beta_mat = beta_mat,
#'   weights_list = local_weights,
#'   point = TRUE,
#'   n_draws = 500,
#'   n_sims = 50
#' )
#'
#' head(dens[[1]])
#'
#' ## ---------------------------------------------------------
#' ## Example: Posterior-sample Beta convolution
#' ## ---------------------------------------------------------
#'
#' J <- 100
#'
#' alpha_post <- cbind(
#'   Total = rgamma(J, 8, 1),
#'   A = rgamma(J, 5, 1),
#'   B = rgamma(J, 6, 1),
#'   A1 = rgamma(J, 2, 1),
#'   A2 = rgamma(J, 4, 1),
#'   B1 = rgamma(J, 3, 1),
#'   B2 = rgamma(J, 5, 1)
#' )
#'
#' beta_post <- cbind(
#'   Total = rgamma(J, 4, 1),
#'   A = rgamma(J, 3, 1),
#'   B = rgamma(J, 2, 1),
#'   A1 = rgamma(J, 6, 1),
#'   A2 = rgamma(J, 5, 1),
#'   B1 = rgamma(J, 7, 1),
#'   B2 = rgamma(J, 4, 1)
#' )
#'
#' dens_post <- Beta_convolution(
#'   S = S,
#'   alpha_mat = alpha_post,
#'   beta_mat = beta_post,
#'   weights_list = local_weights,
#'   point = FALSE,
#'   n_draws = J,
#'   n_sims = 50
#' )
#'
#' head(dens_post[[1]])
#'
#' @export

Beta_convolution <- function(S, alpha_mat, beta_mat, weights_list, point,
                             n_draws = 2000, n_sims = 100, z_values = NULL) {

  ## --------------------------------------------------------
  ## Density grid
  ## --------------------------------------------------------

  if (is.null(z_values)) {
    z_values <- seq(0, 1, length.out = 1000)
  }

  ## --------------------------------------------------------
  ## Build aggregation matrix
  ## --------------------------------------------------------

  W <- build_weight_matrix(
    S = S,
    weights_list = weights_list
  )

  ## --------------------------------------------------------
  ## Input checks
  ## --------------------------------------------------------

  check_beta_convolution_inputs(
    S = S,
    W = W,
    alpha_mat = alpha_mat,
    beta_mat = beta_mat,
    point = point,
    n_draws = n_draws,
    n_sims = n_sims,
    z_values = z_values
  )

  ## --------------------------------------------------------
  ## Bottom and aggregated nodes
  ## --------------------------------------------------------

  bottom_nodes <- colnames(S)

  agg_nodes <- setdiff(rownames(W), bottom_nodes)

  W_agg <- W[agg_nodes, , drop = FALSE]

  n_bottom <- length(bottom_nodes)
  n_agg <- nrow(W_agg)

  ## --------------------------------------------------------
  ## Density estimation
  ## --------------------------------------------------------

  dens_list <- vector("list", n_agg)

  for (i in seq_len(n_agg)) {

    ## ------------------------------------------------------
    ## Current aggregated node
    ## ------------------------------------------------------

    parent_node <- agg_nodes[i]

    ## ------------------------------------------------------
    ## Bottom nodes contributing to aggregation
    ## ------------------------------------------------------

    weights <- W_agg[i, ]

    valid_bottom <- which(weights > 0)

    node_names <- bottom_nodes[valid_bottom]

    weights <- weights[valid_bottom]

    ## ------------------------------------------------------
    ## Extract Beta parameters
    ## ------------------------------------------------------

    if (isTRUE(point)) {

      alpha_input <- alpha_mat[match(node_names, names(alpha_mat))]

      beta_input <- beta_mat[match(node_names, names(beta_mat))]

    } else {

      alpha_input <- alpha_mat[ , match(node_names, colnames(alpha_mat)),drop = FALSE]

      beta_input <- beta_mat[ ,match(node_names, colnames(beta_mat)), drop = FALSE]
    }

    ## ------------------------------------------------------
    ## Weighted samples
    ## ------------------------------------------------------

    weighted_samps <- array(NA, dim = c(n_sims, n_draws, length(node_names)))

    ## ------------------------------------------------------
    ## Point-estimate mode
    ## ------------------------------------------------------

    if (isTRUE(point)) {

      for (m in seq_len(length(node_names))) {

        weighted_samps[, , m] <- matrix(ExtDist::rBeta_ab(n_sims * n_draws,
                                                          alpha_input[m],
                                                          beta_input[m],
                                                          0,
                                                          weights[m]),
                                        nrow = n_sims, ncol = n_draws)
      }

      dens <- Beta_convolution_density_point_parallel(
        z_values = z_values,
        alpha_point = alpha_input,
        beta_point = beta_input,
        weighted_samps = weighted_samps,
        weights = weights
      )

    } else {

      ## ----------------------------------------------------
      ## Posterior-sample mode
      ## ----------------------------------------------------

      for (m in seq_len(length(node_names))) {

        for (d in seq_len(n_draws)) {

          weighted_samps[, d, m] <-
            ExtDist::rBeta_ab(
              n_sims,
              alpha_input[d, m],
              beta_input[d, m],
              0,
              weights[m]
            )
        }
      }

      dens <- Beta_convolution_density_parallel(
        z_values = z_values,
        alpha_matrix = alpha_input,
        beta_matrix = beta_input,
        weighted_samps = weighted_samps,
        weights = weights
      )
    }

    ## ------------------------------------------------------
    ## Store densities
    ## ------------------------------------------------------

    dens_df <- tibble::tibble(
      Node = parent_node,
      Z = z_values,
      Density = dens
    )

    dens_list[[i]] <- dens_df
  }

  names(dens_list) <- agg_nodes

  return(dens_list)
}
