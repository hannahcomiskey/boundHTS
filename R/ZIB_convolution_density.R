#' Vectorized Zero Inflated 4 parameter Beta density
#'
#' @param Y_mc Monte Carlo draws of weighted bottom-series samples
#' @param phi_array Monte Carlo draws of phi
#' @param zi_array Monte Carlo draws of zero inflation
#' @param weights node weights (vector of length n_nodes)
#' @param z_values evaluation points
#' @param n_mc number of Monte Carlo samples
#'
#'

ZIB_convolution_density <- function(Y_mc,
                                phi_array,
                                zi_array,
                                weights,
                                z_values,
                                n_mc) {

  n_draws <- dim(Y_mc)[1]

  # Resample posterior draws ONCE
  draw_id <- sample(seq_len(n_draws), n_mc, replace = TRUE)

  phi_mc <- phi_array[draw_id, , drop=FALSE]
  zi_mc <- zi_array[draw_id, , drop=FALSE]

  future::plan(future::sequential)

  Density <- future.apply::future_sapply(
    z_values,
    function(z) {
      library(boundHTS)
      mean(
        dZIB_4p(z = z,
                Y_mc = Y_mc,
                phi_mc = phi_mc,
                zi_mc = zi_mc,
                upper = weights), na.rm = TRUE)
    },
    future.seed = TRUE
  )

  norm_dens <- Density / pracma::trapz(z_values, Density)

  return(norm_dens)

}
