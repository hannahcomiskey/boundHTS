#' Parallelized normalised predictive density over vector z
#'
#' @param S Binary summation matrix with rows corresponding to all hierarchy nodes and columns corresponding to bottom-level nodes.
#' @param alpha_mat Numeric vector (point = TRUE) or matrix (point = FALSE) of Beta shape1/alpha parameters for all hierarchy nodes.
#'  Columns correspond to hierarchy nodes and must match
#' `rownames(W)`. When `point = TRUE`, `alpha_mat` should be a numeric vector of length equal to the number of hierarchy nodes.
#'  When `point = FALSE`, `alpha_mat` should be a matrix with:
#'\itemize{
#' \item rows corresponding to posterior draws,
#' \item columns corresponding to hierarchy nodes.
#' }
#'@param beta_mat Numeric vector (point = TRUE) or matrix (point = FALSE) of Beta shape2/beta parameters for all hierarchy nodes. Columns correspond to hierarchy nodes and must match `rownames(W)`.
#'  When `point = TRUE`, `beta_mat` should be a numeric vector of length equal to the number of hierarchy nodes.
#'  When `point = FALSE`, `beta_mat` should be a matrix with:
#'  \itemize{
#'  \item rows corresponding to posterior draws,
#'  \item columns corresponding to hierarchy nodes.
#' }
#' @param weights_list A list of weights that are used to combine to form the aggregated density Z. The length of the list is the number of tiers (excluding the top tier) in the hierarchy.
#' @param point A true/false indicator to denote whether you are using
#' point estimates (point=TRUE) or posterior samples (point=FALSE) of the Beta parameters.
#' @param n_draws The number of draws to extract from the four-parameter Beta distribution. Default to 2000.
#' @param n_sims The number of simulations to extract from the four-parameter Beta distribution per draw. Default to 100.
#' @param z_values Evaluation points. Defaults to be NULL and in-built grid is (0,1) over 1000 evenly spaced points.
#' @details
#' A wrapper function for calculating Monte Carlo estimates for the aggregate density Z using the Beta distribution.
#' Convolution is only performed for hierarchy levels containing more than one node. Levels with a single node are skipped, since no aggregation is required.
#' @return A list of aggregate densities Z over a grid of values using a convolution of Beta distributions for each of the aggregated nodes in the hierarchy.
#' @examples
#' ## ---------------------------------------------------------------
#' ## Example: Point-estimate Beta convolution
#' ## ---------------------------------------------------------------
#'
#' set.seed(123)
#'
#' ## Define a simple two-level hierarchy
#'  groups <- list(
#'   bottom = c("A", "B"),   # bottom level
#'   top = c("Total")     # aggregated level
#' )
#'
#' ## Point estimates for Beta parameters
#' alpha_list <- list(
#'   c(2, 5),
#'   c(7)
#' )
#'
#' beta_list <- list(
#'   c(6, 3),
#'   c(4)
#' )
#'
#' ## Aggregation weights
#' weights_list <- list(
#'   c(0.4, 0.6),
#'   1
#' )
#'
#' ## Estimate densities
#' dens <- Beta_convolution(
#'   groups = groups,
#'   alpha_list = alpha_list,
#'   beta_list = beta_list,
#'   weights_list = weights_list,
#'   point = TRUE,
#'   n_draws = 500,
#'   n_sims = 50
#' )
#'
#' ## Inspect first density estimate
#' head(dens[[1]])
#'
#' ## ---------------------------------------------------------
#' ## Example hierarchy
#' ## ---------------------------------------------------------
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
#'   Total = c(A = 0.2, B = 0.8),
#'   A     = c(A1 = 0.4, A2 = 0.6),
#'   B     = c(B1 = 0.3, B2 = 0.7)
#' )
#'
#' J <- 100
#'
#' alpha_list_post <- list(
#' cbind(rgamma(J, shape = 2, rate = 1), rgamma(J, shape = 5, rate = 1)),
#' as.matrix(rgamma(J, shape = 2, rate = 1))
#' )
#'
#' beta_list_post <- list(
#' cbind(rgamma(J, shape = 6, rate = 1), rgamma(J, shape = 3, rate = 1)),
#' as.matrix(rgamma(J, shape = 2, rate = 1))
#' )
#'
#'
#' weights_post <- list(c(0.4, 0.6), 1)
#'
#' dens_post <- Beta_convolution(
#'   groups = groups,
#'   alpha_list = alpha_list_post,
#'   beta_list = beta_list_post,
#'   weights_list = weights_post,
#'   point = FALSE,
#'   n_draws = J,
#'   n_sims = 50
#' )
#'
#' head(dens_post[[1]])
#' @export

old_Beta_convolution <- function(S,
                             alpha_list,
                             beta_list,
                             weights_list,
                             point,
                             n_draws = 2000,
                             n_sims = 100,
                             z_values = NULL) {

  ## --------------------------------------------------------
  ## Check inputs
  ## --------------------------------------------------------

  if(is.null(z_values)) {
    z_values <- seq(0, 1, length.out = 1000)
  }

  W <- build_weight_matrix(S = S, weights_list = weights_list)

  ## Input checks
  check_beta_convolution_inputs(
    S = S,
    alpha_list = alpha_list,
    beta_list = beta_list,
    weights_list = weights_list,
    point = point,
    n_draws = n_draws,
    n_sims = n_sims,
    z_values = z_values
  )

  ## ---------------------------------------------------------
  ## Bottom and aggregated nodes
  ## ---------------------------------------------------------

  bottom_nodes <- colnames(S)
  agg_nodes <- setdiff(rownames(W), bottom_nodes)
  W_agg <- W[agg_nodes, ,drop = FALSE]
  n_bottom <- length(bottom_nodes)
  n_agg <- nrow(W_agg)
  n_nodes  <- nrow(S)

  ## --------------------------------------------------------
  ## Density estimation
  ## --------------------------------------------------------

  dens_list <- vector("list", n_agg)

  for(i in seq_len(n_agg)) {
    x <- valid_idx[i]  # actual level index
    node_names  <- groups[[x]]

    ## ---------------------------------------------------------
    ## Safe parent naming
    ## ---------------------------------------------------------
    parent_node <- if (!is.null(names(groups))) {
      names(groups)[x+1]
    } else {
      paste0("Level_", x)
    }

    ## ---------------------------------------------------------
    ## Beta distribution parameters
    ## ---------------------------------------------------------
    alpha_input <- alpha_list[[x]]
    beta_input <- beta_list[[x]]
    weights <- weights_list[[x]]

    ## Weighted samples
    weighted_samps <- array(NA, dim = c(n_sims, n_draws, length(node_names)))

    for (m in seq_along(node_names)) {

      weighted_samps[, , m] <- matrix(
        ExtDist::rBeta_ab(n_sims * n_draws,
                          alpha_input[m], beta_input[m],
                          0, weights[m]), nrow = n_sims, ncol = n_draws
        )
    }

    if (isTRUE(point)) {

      dens <- Beta_convolution_density_point_parallel(
        z_values = z_values,
        alpha_point = alpha_input,
        beta_point = beta_input,
        weighted_samps = weighted_samps,
        weights = weights
      )

    } else {

      dens <- Beta_convolution_density_parallel(
        z_values = z_values,
        alpha_matrix = alpha_input,
        beta_matrix = beta_input,
        weighted_samps = weighted_samps,
        weights = weights
      )
    }
    dens_df <- tibble::tibble(Node = parent_node, Z = z_values, Density = dens)
    dens_list[[x]] <- dens_df
  }

  return(dens_list)
}
