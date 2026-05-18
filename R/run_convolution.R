#' Wrapper for generating convolution
#'
#' @param family The name of the probability distribution used in the convolution. Must be a string of one of these options: c("Poisson", "NB", "Beta", "ZIB", "ZOIB").
#' @param groups List of nodes at each level in series. The length of the list is the number of tiers in the hierarchy.
#' @param point A true/false indicator to denote whether you are using point estimates (point=TRUE) or posterior samples (point=FALSE) of the Poisson parameters.
#' @param z_values A grid of evaluation points. If the family is in the Beta-family, this can be left NULL and it defaults to a vector of length 1000 over the unit interval.
#' @param ... Parameters specific to the parameter family. See "Poisson_convolution", "Beta_convolution", "ZIB_convolution", and "ZOIB_convolution" for family-specific details.
#' @details
#' A wrapper function for calculating Monte Carlo estimates for the aggregate
#' density Z using point estimates or posterior samples of lambda.
#' @examples
#' ## ---------------------------------------------------------------
#' ## Example 1: Poisson convolution using point estimates
#' ## ---------------------------------------------------------------
#'
#' set.seed(123)
#'
#' ## Hierarchy:
#' ## State-by-sex -> State -> National
#'
#' groups <- list(
#'   NSW = c("NSW_Male", "NSW_Female"),
#'   VIC = c("VIC_Male", "VIC_Female"),
#'   State = c("NSW", "VIC"),
#'   National = c("AUS")
#' )
#'
#' ## Evaluation grid
#' z_values <- 0:60
#'
#' ## Point estimates for lambda
#' lambda_list <- list(
#'   c(14, 12),
#'   c(13, 13),
#'   c(26, 26),
#'   52
#' )
#'
#' ## Run convolution
#' dens_pois <- run_convolution(
#'   family = "Poisson",
#'   groups = groups,
#'   point = TRUE,
#'   z_values = z_values,
#'   lambda_list = lambda_list
#' )
#'
#' ## Inspect results
#' head(dens_pois[[1]])
#'
#'
#' ## ---------------------------------------------------------------
#' ## Example 2: Beta convolution using point estimates
#' ## ---------------------------------------------------------------
#'
#' set.seed(123)
#'
#' groups_beta <- list(
#'   State = c("NSW", "VIC"),
#'   National = c("AUS")
#' )
#'
#' alpha_list <- list(
#'   c(2, 5),
#'   7
#' )
#'
#' beta_list <- list(
#'   c(6, 3),
#'   4
#' )
#'
#' weights_list <- list(
#'   c(0.4, 0.6),
#'   1
#' )
#'
#' dens_beta <- run_convolution(
#'   family = "Beta",
#'   groups = groups_beta,
#'   point = TRUE,
#'   z_values = NULL,
#'   alpha_list = alpha_list,
#'   beta_list = beta_list,
#'   weights_list = weights_list,
#'   n_draws = 500,
#'   n_sims = 50
#' )
#'
#' ## Inspect results
#' head(dens_beta[[1]])
#'
#'
#' ## ---------------------------------------------------------------
#' ## Example 3: Zero-one-inflated Beta convolution
#' ## ---------------------------------------------------------------
#'
#' set.seed(123)
#'
#' zoi_list <- list(
#'   c(0.10, 0.15),
#'   0.05
#' )
#'
#' coi_list <- list(
#'   c(0.05, 0.10),
#'   0.02
#' )
#'
#' dens_zoib <- run_convolution(
#'   family = "ZOIB",
#'   groups = groups_beta,
#'   point = TRUE,
#'   z_values = NULL,
#'   alpha_list = alpha_list,
#'   beta_list = beta_list,
#'   zoi_list = zoi_list,
#'   coi_list = coi_list,
#'   weights_list = weights_list,
#'   n_draws = 500,
#'   n_sims = 50
#' )
#'
#' ## Inspect results
#' head(dens_zoib[[1]])
#' @return The aggregate density Z over a grid of values.
#' @export
#'

run_convolution <- function(family,
                            groups,
                            point,
                            z_values,
                            ...) {

  args <- list(...)

  if (!family %in% c("Poisson", "Beta", "ZIB", "ZOIB")) {
    stop(
      "Family name not recognised. Please refer to help file for naming.",
      call. = FALSE
    )
  }

  if (family == "Poisson") {

    dens_list <- Poisson_convolution(
      groups = groups,
      z_values = z_values,
      lambda_list = args$lambda_list,
      point = point
    )

  } else if (family == "Beta") {

    dens_list <- Beta_convolution(
      groups = groups,
      alpha_list = args$alpha_list,
      beta_list = args$beta_list,
      weights_list = args$weights_list,
      point = point,
      n_draws = default_settings(args$n_draws, 2000),
      n_sims = default_settings(args$n_sims, 100),
      z_values = z_values
    )

  } else if (family == "ZIB") {

    dens_list <- ZIB_convolution(
      groups = groups,
      alpha_list = args$alpha_list,
      beta_list = args$beta_list,
      zi_list = args$zi_list,
      weights_list = args$weights_list,
      point = point,
      n_draws = default_settings(args$n_draws, 2000),
      n_sims = default_settings(args$n_sims, 100),
      z_values = z_values
    )

  } else if (family == "ZOIB") {

    dens_list <- ZOIB_convolution(
      groups = groups,
      alpha_list = args$alpha_list,
      beta_list = args$beta_list,
      zoi_list = args$zoi_list,
      coi_list = args$coi_list,
      weights_list = args$weights_list,
      point = point,
      n_draws = default_settings(args$n_draws, 2000),
      n_sims = default_settings(args$n_sims, 100),
      z_values = z_values
    )
  }

  return(dens_list)
}
