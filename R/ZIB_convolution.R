#' Parallelized normalised predictive density over vector z
#'
#' @param groups List of nodes at each level in series. The length of the list is the number of tiers in the hierarchy.
#' @param alpha_list A list of point estimates (point=TRUE) or matrix (point=FALSE) of
#' Beta distribution shape1/alpha parameters for each element per level over each of the J observations (J rows x  X columns). The length of the list is the number of levels in the hierarchy.
#' @param beta_list A list of point estimates (point=TRUE) or matrix (point=FALSE) of
#' Beta distribution shape2/beta parameters for each element per level over each of the J observations (J rows x  X columns). The length of the list is the number of levels in the hierarchy.
#' @param zi_list A list of point estimates (point=TRUE) or matrix (point=FALSE) of zero-inflation probability parameters for each element per level over each of the J observations (J rows x  X columns). The length of the list is the number of levels in the hierarchy.
#' @param weights_list A list of weights that are used to combine to form the aggregated density Z. The length of the list is the number of levels (excluding the top level) in the hierarchy.
#' @param point A true/false indicator to denote whether you are using
#' point estimates (point=TRUE) or posterior samples (point=FALSE) of the zero-inflated Beta parameters.
#' @param n_draws The number of draws to extract from the four-parameter zero-inflated Beta distribution. Default to 2000.
#' @param n_sims The number of simulations to extract from the four-parameter zero-inflated Beta distribution per draw. Default to 100.
#' @param z_values Evaluation points. Defaults to be NULL and in-built grid is (0,1) over 1000 evenly spaced points.
#' @details
#' A wrapper function for calculating Monte Carlo estimates for the aggregate density Z using the zero-inflated Beta distribution.
#' @return A list of aggregate densities Z over a grid of values using a convolution of zero-inflated Beta distributions for each of the aggregated nodes in the hierarchy.
#' @examples
#' ## ---------------------------------------------------------------
#' ## Example: Point-estimate zero-inflated Beta convolution
#' ## ---------------------------------------------------------------
#'
#' set.seed(123)
#'
#' ## Define a simple two-level hierarchy
#' groups <- list(
#'   c("A", "B"),   # bottom level
#'   c("Total")     # aggregated level
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
#' ## Zero-inflation probabilities
#' zi_list <- list(
#'   c(0.10, 0.20),
#'   c(0.05)
#' )
#'
#' ## Aggregation weights
#' weights_list <- list(
#'   c(0.4, 0.6),
#'   1
#' )
#'
#' ## Estimate densities
#' dens <- ziB_convolution(
#'   groups = groups,
#'   alpha_list = alpha_list,
#'   beta_list = beta_list,
#'   zi_list = zi_list,
#'   weights_list = weights_list,
#'   point = TRUE,
#'   n_draws = 500,
#'   n_sims = 50
#' )
#'
#' ## Inspect first density estimate
#' head(dens[[1]])
#'
#' ## ---------------------------------------------------------------
#' ## Example: Posterior-sample zero-inflated Beta convolution
#' ## ---------------------------------------------------------------
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
#' zi_list_post <- list(
#' cbind(rbeta(J, 2, 20), rbeta(J, 3, 15)),
#' as.matrix(rbeta(J, 2, 20)))
#'
#' weights_list <- list(
#'   c(0.4, 0.6),
#'   1
#' )
#'
#' dens_post <- ZIB_convolution(
#'   groups = list(c("A", "B")),
#'   alpha_list = alpha_list_post,
#'   beta_list = beta_list_post,
#'   zi_list = zi_list_post,
#'   weights_list = weights_list,
#'   point = FALSE,
#'   n_draws = 500,
#'   n_sims = 50
#' )
#'
#' head(dens_post[[1]])
#' @export
#'

ZIB_convolution <- function(groups, alpha_list, beta_list, zi_list, weights_list, point, n_draws=2000, n_sims = 100, z_values=NULL) {
  if(is.null(z_values)==TRUE) {
    z_values <- seq(0, 1, length.out = 1000) # density grid
  }

  # Check data inputs for validity
  check_beta_convolution_inputs(groups, alpha_list, beta_list, weights_list, point, n_draws, n_sims, z_values)

  if (!is.list(zi_list)) {
    stop("'zi_list' must be a list.")
  }

  if (length(groups) != length(zi_list)) {
    stop("'groups' and 'zi_list' must have the same length.")
  }

  valid_idx <- which(vapply(groups, length, integer(1)) > 1L)
  dens_list <- list() # set up list to save results in

  for (i in seq_along(valid_idx)) {

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
    ## ZIB distribution parameters
    ## ---------------------------------------------------------
    alpha_input <- alpha_list[[x]]
    beta_input <- beta_list[[x]]
    zi_input <- zi_list[[x]]
    weights <- weights_list[[x]]

    ## Get weighted samples
    weighted_samps <- array(NA, dim = c(n_sims, n_draws, length(node_names)))

    for (m in seq_along(node_names)) {
      weighted_samps[,,m] <- matrix(ExtDist::rBeta_ab(n_sims * n_draws,
                                                      alpha_input[m],
                                                      beta_input[m],
                                                      0, weights[m]),
                                    nrow=n_sims, ncol= n_draws)
    }
    if(point==TRUE & is.vector(alpha_input)==TRUE & is.vector(beta_input)==TRUE & is.vector(zi_input)==TRUE) {
      dens <- ZIB_convolution_density_point_parallel(z_values = z_values,
                                                     alpha_point = alpha_input,
                                                     beta_point = beta_input,
                                                     zi_point = zi_input,
                                                     weighted_samps = weighted_samps,
                                                     weights = weights)
    }
    if(point==FALSE & is.matrix(alpha_input)==TRUE & is.matrix(beta_input)==TRUE & is.matrix(zi_input)==TRUE) {
      dens <- ZIB_convolution_density_parallel(z_values = z_values,
                                               alpha_matrix = alpha_input,
                                               beta_matrix = beta_input,
                                               zi_matrix = zi_input,
                                               weighted_samps = weighted_samps,
                                               weights = weights)
    }
    dens_df <- tibble::tibble(Node = parent_node, Z = z_values, Density = dens)
    dens_list[[x]] <- dens_df
  }
  return(dens_list)
}
