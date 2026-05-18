#' Parallelized normalised predictive density over vector z
#'
#' @param groups List of nodes at each level in series. The list must go from bottom nodes to top nodes, and have each element of the list named. The length of the list is the number of tiers in the hierarchy.
#' @param shape_list A list of point estimates (point=TRUE) or matrix (point=FALSE) of
#' Beta distribution shape1/alpha parameters for each element at the level over each of the J observations (J rows x  X columns). The length of the list is the number of tiers in the hierarchy.
#' @param rate_list A list of point estimates (point=TRUE) or matrix (point=FALSE) of
#' Beta distribution shape2/beta parameters for each element at the level over each of the J observations (J rows x  X columns). The length of the list is the number of tiers in the hierarchy.
#' @param point A true/false indicator to denote whether you are using
#' point estimates (point=TRUE) or posterior samples (point=FALSE) of the Beta parameters.
#' @param n_draws The number of draws to extract from the four-parameter Beta distribution. Default to 2000.
#' @param n_sims The number of simulations to extract from the four-parameter Beta distribution per draw. Default to 100.
#' @param z_values Evaluation points.
#' @details
#' A wrapper function for calculating Monte Carlo estimates for the aggregate density Z using the Beta distribution.
#' Convolution is only performed for hierarchy levels containing more than one node. Levels with a single node are skipped, since no aggregation is required.
#' @return A list of aggregate densities Z over a grid of values using a convolution of Beta distributions for each of the aggregated nodes in the hierarchy.
#' @examples
#' ## ---------------------------------------------------------------
#' ## Example: Point-estimate Gamma convolution
#' ## ---------------------------------------------------------------
#' set.seed(123)
#' z_values <- seq(0, 50, length.out = 1000)
#'
#' ## Define a simple two-level hierarchy
#'  groups <- list(
#'   bottom = c("AA", "AB"),   # bottom level
#'   top = c("Total")     # aggregated level
#' )
#'
#' ## Point estimates for Gamma parameters
#' shape_list <- list(
#'   c(2, 5),
#'   c(7)
#' )
#'
#' rate_list <- list(
#'   c(6, 3),
#'   c(4)
#' )
#'
#'
#' ## Estimate densities
#' dens <- Gamma_convolution(
#'   groups = groups,
#'   shape_list = shape_list,
#'   rate_list = rate_list,
#'   z_values = z_values,
#'   point = TRUE,
#'   n_draws = 500,
#'   n_sims = 50
#' )
#'
#' ## Inspect first density estimate
#' head(dens[[1]])
#'
#' ## ---------------------------------------------------------------
#' ## Example: Posterior-sample Gamma convolution
#' ## ---------------------------------------------------------------
#' z_values <- seq(0, 50, length.out = 1000)
#'
#' groups <- list(
#'   bottom = c("A", "B"),   # bottom level
#'   top = c("Total")     # aggregated level
#' )
#'
#' J <- 100
#'
#' shape_list_post <- list(
#' cbind(rgamma(J, shape = 2, rate = 1), rgamma(J, shape = 5, rate = 1)),
#' as.matrix(rgamma(J, shape = 2, rate = 1))
#' )
#'
#' rate_list_post <- list(
#' cbind(rgamma(J, shape = 6, rate = 1), rgamma(J, shape = 3, rate = 1)),
#' as.matrix(rgamma(J, shape = 2, rate = 1))
#' )
#'
#'
#' dens_post <- Gamma_convolution(
#'   groups = groups,
#'   shape_list = shape_list_post,
#'   rate_list = rate_list_post,
#'   z_values = z_values,
#'   point = FALSE,
#'   n_draws = J,
#'   n_sims = 50
#' )
#'
#' head(dens_post[[1]])
#' @export

Gamma_convolution <- function(groups,
                             shape_list,
                             rate_list,
                             z_values,
                             point,
                             n_draws = NULL,
                             n_sims = 100) {


  ## Input checks
  check_Gamma_convolution_inputs(
    groups = groups,
    shape_list = shape_list,
    rate_list = rate_list,
    point = point,
    n_draws = n_draws,
    n_sims = n_sims,
    z_values = z_values
  )

  valid_idx <- which(vapply(groups, length, integer(1)) > 1L)

  dens_list <- list()

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
    ## Gamma distribution parameters
    ## ---------------------------------------------------------
    shape_input <- shape_list[[x]]
    rate_input  <- rate_list[[x]]

    if (isTRUE(point)) {

      n_draws <- ifelse(is.null(n_draws), 2000, n_draws)

      ## Bottom samples
      bottom_samps <- array(
        NA,
        dim = c(n_sims, n_draws, length(node_names))
      )

      for (m in seq_along(node_names)) {

        bottom_samps[, , m] <- matrix(
          stats::rgamma(
            n = n_sims * n_draws,
            shape = shape_input[,m],
            rate = rate_input[,m]),
          nrow = n_sims,
          ncol = n_draws
        )
      }

      dens <- Gamma_convolution_density_point_parallel(
        z_values = z_values,
        shape_point = shape_input,
        rate_point = rate_input,
        bottom_samps = bottom_samps)

    } else {

      n_draws <- nrow(shape_input)

      ## Bottom samples
      bottom_samps <- array(NA, dim = c(n_sims, n_draws, length(node_names)))

      for (m in 1:length(node_names)) {
        for(s in 1:n_draws) {
          bottom_samps[, s, m] <- stats::rgamma(n = n_sims,
                                                shape = shape_input[s,m],
                                                rate = rate_input[s,m])
        }
      }

      dens <- Gamma_convolution_density_parallel(
        z_values = z_values,
        shape_matrix = shape_input,
        rate_matrix = rate_input,
        bottom_samps = bottom_samps)
    }
    dens_df <- tibble::tibble(Node = parent_node, Z = z_values, Density = dens)
    dens_list[[x]] <- dens_df
  }

  return(dens_list)
}
