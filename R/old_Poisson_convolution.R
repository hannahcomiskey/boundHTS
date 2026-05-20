#' Parallelized normalised predictive Poisson density over vector z
#'
#' @param groups List of nodes at each level in series. The list must go from bottom nodes to top nodes, and have each element of the list named. The length of the list is the number of tiers in the hierarchy.
#' @param z_values Evaluation points.
#' @param lambda_list A list of point estimates (point=TRUE) or matrix (point=FALSE) of lambda parameters for the Poisson distribution for each bottom series observation.
#' @param point A true/false indicator to denote whether you are using point estimates (point=TRUE) or posterior samples (point=FALSE) of the Poisson parameters.
#' @details
#' A wrapper function for calculating Monte Carlo estimates for the aggregate
#' density Z using point estimates or posterior samples of lambda.
#'
#' @return The aggregate density Z over a grid of values.
#' @examples
#' ## ---------------------------------------------------------------
#' ## Example: Point-estimate Poisson convolution
#' ## ---------------------------------------------------------------
#'
#' set.seed(123)
#'
#' ## Define hierarchy from bottom to top
#'  groups <- list(NSW = c("NSW_Male", "NSW_Female"),
#'  VIC = c("VIC_Male", "VIC_Female"),
#'  State = c("NSW", "VIC"),
#'  National = c("AUS"))
#'
#' ## Evaluation grid
#' z_values <- 0:40
#'
#' ## Point estimates for lambda
#'  lambda_list <- list(c(14, 12), # State-by-sex NSW
#'  c(13, 13), # State-by-sex VIC
#'  c(22, 28), # State
#'  52 # National
#'  )
#'
#' ## Estimate aggregate densities
#' dens <- Poisson_convolution(
#'   groups = groups,
#'   z_values = z_values,
#'   lambda_list = lambda_list,
#'   point = TRUE
#' )
#'
#' ## Inspect first level density
#' head(dens[[1]])
#'
#' ## ---------------------------------------------------------------
#' ## Example: Posterior-sample Poisson convolution
#' ## ---------------------------------------------------------------
#'
#' set.seed(123)
#'
#' ## Number of posterior draws
#' J <- 500
#' ## Posterior samples for State-by-sex lambdas
#' lambda_state_sex_NSW <- cbind(rgamma(J, shape = 12, rate = 1), rgamma(J, shape = 10, rate = 1))
#' lambda_state_sex_VIC <- cbind(rgamma(J, shape = 15, rate = 1), rgamma(J, shape = 13, rate = 1))
#' ## Posterior samples for State lambdas
#' lambda_state <- cbind(rgamma(J, shape = 22, rate = 1), rgamma(J, shape = 28, rate = 1))
#'
#' ## Posterior samples for National lambda
#' lambda_national <- matrix(rgamma(J, shape = 50, rate = 1), ncol = 1)
#'
#' ## Store posterior samples in a list
#' lambda_list_post <- list(lambda_state_sex_NSW, lambda_state_sex_VIC, lambda_state, lambda_national)
#'
#' ## Estimate aggregate densities
#' dens_post <- Poisson_convolution(
#'   groups = groups,
#'   z_values = z_values,
#'   lambda_list = lambda_list_post,
#'   point = FALSE
#' )
#'
#' ## Inspect density estimate
#' head(dens_post[[1]])
#' @export

old_Poisson_convolution <- function(groups, z_values, lambda_list, point) {

  dens_list <- list()

  check_poisson_convolution_inputs(groups, z_values, lambda_list, point)

  ## ---------------------------------------------------------
  ## Identify valid levels (skip top if singleton)
  ## ---------------------------------------------------------
  valid_idx <- which(vapply(groups, length, integer(1)) > 1L)

  for (i in seq_along(valid_idx)) {

    x <- valid_idx[i]  # actual level index

    node_names  <- groups[[x]]
    lambda_input <- lambda_list[[x]]

    ## ---------------------------------------------------------
    ## Safe parent naming
    ## ---------------------------------------------------------
    parent_node <- if (!is.null(names(groups))) {
      names(groups)[x+1]
    } else {
      paste0("Level_", x)
    }

    ## ---------------------------------------------------------
    ## Convolution
    ## ---------------------------------------------------------
    if (point && is.vector(lambda_input)) {

      dens <- Poisson_convolution_density_point_parallel(
        z_values = z_values,
        lambda_vector = lambda_input
      )

    } else if (!point && is.matrix(lambda_input)) {

      dens <- Poisson_convolution_density_parallel(
        z_values = z_values,
        lambda_matrix = lambda_input
      )

    } else {
      stop("lambda input type does not match 'point' argument")
    }

    dens_list[[i]] <- tibble::tibble(
      Node = parent_node,
      Z = z_values,
      Density = dens
    )
  }

  return(dens_list)
}
