#' Wrapper for generating convolution
#'
#' @param family The name of the probability distribution used in the convolution. Must be a string of one of these options: c("Poisson", "NB", "Beta", "ZIB", "ZOIB").
#' @param groups List of nodes at each level in series. The length of the list is the number of tiers in the hierarchy.
#' @param point A true/false indicator to denote whether you are using point estimates (point=TRUE) or posterior samples (point=FALSE) of the Poisson parameters.
#' @param z_values A grid of evaluation points. If the family is in the Beta-family, this can be left NULL and it defaults to a vector of length 1000 over the unit interval.
#' @details
#' A wrapper function for calculating Monte Carlo estimates for the aggregate
#' density Z using point estimates or posterior samples of lambda.
#'
#' @return The aggregate density Z over a grid of values.
#' @export
#'

run_convolution <- function(family, groups, point, z_values, ...) {
  if(family %in% c("Poisson", "Beta", "ZIB", "ZOIB")) { # "NB" to be done
    if(family=="Poisson") {
      dens_list <- Poisson_convolution(groups, z_values, lambda_list, point)
    } else{
      # if(family=="NB") {
      #   dens_list <- NB_convolution(groups, alpha_list, beta_list, weights_list, point, n_draws=2000, n_sims = 100, z_values=NULL)
      # } else{
        if(family=="Beta") {
          dens_list <- Beta_convolution(groups, alpha_list, beta_list, weights_list, point, n_draws=2000, n_sims = 100, z_values=NULL)
        } else{
          if(family=="ZIB") {
            dens_list <- ZIB_convolution(groups, alpha_list, beta_list, zi_list, weights_list, point, n_draws=2000, n_sims = 100, z_values=NULL)
          } else {
            dens_list <- ZOIB_convolution(groups, alpha_list, beta_list, zoi_list, coi_list, weights_list, point, n_draws=2000, n_sims = 100, z_values=NULL)
          }
        }
      # }
    }
  } else {
    stop("Family name not recongised. Please refer to Help file for naming.")
  }
}
