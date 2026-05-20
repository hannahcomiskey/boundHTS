#' Parallelized normalised predictive density over vector z
#'
#' @param S Binary summation matrix with rows corresponding to all hierarchy nodes and columns corresponding to bottom-level nodes.
#' @param lambda_mat Numeric vector (point = TRUE) or matrix (point = FALSE) of Poisson shape parameters for all hierarchy nodes.
#' @param point Logical. If TRUE uses point estimates, otherwise posterior draws.
#' @param n_draws Number of posterior draws (if applicable).
#' @param n_sims Number of Monte Carlo simulations per draw.
#' @param z_values Evaluation grid.
#'
#' @return List of aggregated density estimates.
#' @examples
#' ## ---------------------------------------------------------
#' ## Example: Point-estimate Poisson convolution
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
#'
#' lambda_mat <- c(
#'   Total = 8,
#'   A = 5,
#'   B = 6,
#'   A1 = 2,
#'   A2 = 4,
#'   B1 = 3,
#'   B2 = 5
#' )
#'
#' z_values <- seq(0, 50, length.out = 500)
#'
#' dens <- Poisson_convolution(
#'   S = S,
#'   lambda_mat = lambda_mat,
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
#' ## Example: Posterior-sample Poisson convolution
#' ## ---------------------------------------------------------
#'
#' J <- 100
#'
#' lambda_post <- cbind(Total = rPoisson(J, 8, 1),
#'   A = rPoisson(J, 5, 1),
#'   B = rPoisson(J, 6, 1),
#'   A1 = rPoisson(J, 2, 1),
#'   A2 = rPoisson(J, 4, 1),
#'   B1 = rPoisson(J, 3, 1),
#'   B2 = rPoisson(J, 5, 1)
#' )
#'
#'
#' dens_post <- Poisson_convolution(
#'   S = S,
#'   lambda_mat = lambda_post,
#'   z_values = z_values,
#'   point = FALSE,
#'   n_draws = J,
#'   n_sims = 50
#' )
#'
#' head(dens_post[[1]])
#'
#' @export

Poisson_convolution <- function(S, lambda_mat, point,
                              n_draws = NULL, n_sims = 100, z_values) {

  ## --------------------------------------------------------
  ## Input checks
  ## --------------------------------------------------------

  check_Poisson_convolution_inputs(S = S, lambda_mat = lambda_mat, point = point,
                                   n_draws = n_draws, n_sims = n_sims, z_values = z_values)

  ## --------------------------------------------------------
  ## Hierarchy structure
  ## --------------------------------------------------------

  bottom_nodes <- colnames(S)

  agg_nodes <- rownames(S)[rowSums(S) > 1]

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
      lambda_input <- lambda_mat[node_names]
      n_draws <- ifelse(is.null(n_draws), 2000, n_draws)
    } else {
      lambda_input <- lambda_mat[, node_names, drop = FALSE]
      n_draws <- nrow(lambda_mat)
    }

    ## ------------------------------------------------------
    ## Density estimation
    ## ------------------------------------------------------

    if (isTRUE(point)) {

      dens <- Poisson_convolution_density_point_parallel(
        z_values = z_values,
        lambda_vector = lambda_input
      )

    } else {
      dens <- Poisson_convolution_density_parallel(
        z_values = z_values,
        lambda_matrix = lambda_input
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
