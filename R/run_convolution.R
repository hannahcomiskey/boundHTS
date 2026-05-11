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
#'
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
      n_draws = args$n_draws %||% 2000,
      n_sims = args$n_sims %||% 100,
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
      n_draws = args$n_draws %||% 2000,
      n_sims = args$n_sims %||% 100,
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
      n_draws = args$n_draws %||% 2000,
      n_sims = args$n_sims %||% 100,
      z_values = z_values
    )
  }

  return(dens_list)
}
