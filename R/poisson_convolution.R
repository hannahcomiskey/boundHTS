#' Parallelized normalised predictive Poisson density over vector z
#'
#' @param z_values evaluation points
#' @param lambda_input Point estimates (point=TRUE) or matrix (point=FALSE) of
#' lambda parameters for the Poisson distribution for each bottom series observation.
#' @param point A true/false indicator to denote whether you are using
#' point estimates (point=TRUE) or posterior samples (point=FALSE) of the Poisson parameters.
#' @details
#' A wrapper function for calculating Monte Carlo estimates for the aggregate
#' density Z using point estimates or posterior samples of lambda.
#'
#' @return The aggregate density Z over a grid of values.
#' @export

Poisson_convolution <- function(groups, z_values, lambda_input, point) {
  dens_list <- list() # set up list to save results in

  for(x in 1:length(groups)) { # convolution per level
    node_names <- groups[[x]]
    lambda_input <- lambda_list[[x]]
    if(point==TRUE & is.vector(lambda_input)==TRUE) {
      dens <- Poisson_convolution_density_point_parallel(z_values, lambda_input)
    }
    if(point==FALSE & is.matrix(lambda_input)==TRUE) {
      dens <- Poisson_convolution_density_parallel(z_values, lambda_input)
    }
    dens_list[[x]] <- dens
  }
  return(dens_list)
}
