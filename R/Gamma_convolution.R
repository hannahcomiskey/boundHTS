#' Parallelized normalised predictive density over vector z
#'
#' @param S Binary summation matrix with rows corresponding to all hierarchy nodes and columns corresponding to bottom-level nodes.
#' @param shape_mat Numeric vector (point = TRUE) or matrix (point = FALSE) of Gamma shape parameters for all hierarchy nodes.
#' @param rate_mat Numeric vector (point = TRUE) or matrix (point = FALSE) of Gamma rate parameters for all hierarchy nodes.
#' @param point Logical. If TRUE uses point estimates, otherwise posterior draws.
#' @param n_draws Number of posterior draws (if applicable).
#' @param n_sims Number of Monte Carlo simulations per draw.
#' @param z_values Evaluation grid.
#'
#' @return List of aggregated density estimates.
#' @examples
#' ## ---------------------------------------------------------
#' ## Example: Point-estimate Gamma convolution
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
#' colnames(S) <- c("A1", "A2", "B1", "B2")
#' shape_mat <- c(
#'   Total = 8,
#'   A = 5,
#'   B = 6,
#'   A1 = 2,
#'   A2 = 4,
#'   B1 = 3,
#'   B2 = 5
#' )
#'
#' rate_mat <- c(
#'   Total = 4,
#'   A = 3,
#'   B = 2,
#'   A1 = 6,
#'   A2 = 5,
#'   B1 = 7,
#'   B2 = 4
#' )
#'
#' z_values <- seq(0, 50, length.out = 500)
#'
#' dens <- Gamma_convolution(
#'   S = S,
#'   shape_mat = shape_mat,
#'   rate_mat = rate_mat,
#'   z_values = z_values,
#'   point = TRUE,
#'   n_sims = 50
#' )
#'
#' head(dens[[1]])
#'
#'
#' ## ---------------------------------------------------------
#' ## Example: Posterior-sample Gamma convolution
#' ## ---------------------------------------------------------
#'
#' J <- 100
#'
#' shape_post <- cbind(
#'   Total = rgamma(J, 8, 1),
#'   A = rgamma(J, 5, 1),
#'   B = rgamma(J, 6, 1),
#'   A1 = rgamma(J, 2, 1),
#'   A2 = rgamma(J, 4, 1),
#'   B1 = rgamma(J, 3, 1),
#'   B2 = rgamma(J, 5, 1)
#' )
#'
#' rate_post <- cbind(
#'   Total = rgamma(J, 4, 1),
#'   A = rgamma(J, 3, 1),
#'   B = rgamma(J, 2, 1),
#'   A1 = rgamma(J, 6, 1),
#'   A2 = rgamma(J, 5, 1),
#'   B1 = rgamma(J, 7, 1),
#'   B2 = rgamma(J, 4, 1)
#' )
#'
#' dens_post <- Gamma_convolution(
#'   S = S,
#'   shape_mat = shape_post,
#'   rate_mat = rate_post,
#'   z_values = z_values,
#'   point = FALSE,
#'   n_draws = J,
#'   n_sims = 50
#' )
#'
#' head(dens_post[[1]])
#'
#' @export

Gamma_convolution <- function(S,shape_mat,rate_mat, point,
                              n_draws = NULL, n_sims = 100, z_values) {

  ## --------------------------------------------------------
  ## Input checks
  ## --------------------------------------------------------

  check_Gamma_convolution_inputs(
    S = S,
    shape_mat = shape_mat,
    rate_mat = rate_mat,
    point = point,
    n_draws = n_draws,
    n_sims = n_sims,
    z_values = z_values
  )

  ## --------------------------------------------------------
  ## Hierarchy structure
  ## --------------------------------------------------------

  bottom_nodes <- colnames(S)

  agg_nodes <- rownames(S)[
    rowSums(S) > 1
  ]

  dens_list <- vector("list", length(agg_nodes))

  ## --------------------------------------------------------
  ## Loop over aggregated nodes
  ## --------------------------------------------------------

  for (i in seq_along(agg_nodes)) {

    parent_node <- agg_nodes[i]

    node_idx <- which(S[parent_node, ] == 1)

    node_names <- bottom_nodes[node_idx]

    ## ------------------------------------------------------
    ## Extract parameters
    ## ------------------------------------------------------

    if (isTRUE(point)) {

      shape_input <- shape_mat[node_names]
      rate_input  <- rate_mat[node_names]

      n_draws <- ifelse(is.null(n_draws), 2000, n_draws)

    } else {

      shape_input <- shape_mat[, node_names, drop = FALSE]
      rate_input  <- rate_mat[, node_names, drop = FALSE]

      n_draws <- nrow(shape_mat)
    }

    ## ------------------------------------------------------
    ## Simulate bottom-level Gamma samples
    ## ------------------------------------------------------

    bottom_samps <- array(NA, dim = c(n_sims, n_draws, length(node_names)))

    if (isTRUE(point)) {

      for (m in seq_along(node_names)) {

        bottom_samps[, , m] <- matrix( stats::rgamma(n_sims * n_draws,
                                                     shape = shape_input[m],
                                                     rate = rate_input[m]),
                                       nrow = n_sims, ncol = n_draws)
      }
    } else {

      for (m in seq_along(node_names)) {
        for (d in seq_len(n_draws)) {

          bottom_samps[, d, m] <- stats::rgamma(n_sims,
                                                shape = shape_input[d, m],
                                                rate = rate_input[d, m])
        }
      }
    }

    ## ------------------------------------------------------
    ## Density estimation
    ## ------------------------------------------------------

    if (isTRUE(point)) {

      dens <- Gamma_convolution_density_point_parallel(
        z_values = z_values,
        shape_point = shape_input,
        rate_point = rate_input,
        bottom_samps = bottom_samps
      )

    } else {

      dens <- Gamma_convolution_density_parallel(
        z_values = z_values,
        shape_matrix = shape_input,
        rate_matrix = rate_input,
        bottom_samps = bottom_samps
      )
    }

    dens_list[[i]] <- tibble::tibble(
      Node = parent_node,
      Z = z_values,
      Density = dens
    )
  }

  names(dens_list) <- agg_nodes

  return(dens_list)
}
