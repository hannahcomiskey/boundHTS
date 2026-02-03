#' Vectorized ZOIB density
#'
#' @param z_values evaluation points
#' @param n_mc number of Monte Carlo samples
#' @param Y_mc Monte Carlo draws of Y
#' @param phi_array Monte Carlo draws of phi
#' @param zoi_array Monte Carlo draws of zero-one inflation
#' @param coi_array Monte Carlo draws of conditional one inflation
#' @param lower lower bound
#' @param upper upper bound

ZOIB_convolution_density <- function(Y_mc,
                                phi_array,
                                zoi_array,
                                coi_array,
                                weights,
                                z_values,
                                n_mc) {

  n_draws <- dim(phi_array)[1]

  # Resample posterior draws ONCE
  draw_id <- sample(seq_len(n_draws), n_mc, replace = TRUE)

  phi_mc <- phi_array[draw_id, , drop = FALSE]
  zoi_mc <- zoi_array[draw_id, , drop = FALSE]
  coi_mc <- coi_array[draw_id, , drop = FALSE]

  future::plan(future::multisession)

  Density <- future.apply::future_sapply(
    z_values,
    function(z) {
      mean(
        dZOIB_4p(
          z       = z,
          Y_mc    = Y_mc,
          phi_mc  = phi_mc,
          zoi_mc  = zoi_mc,
          coi_mc  = coi_mc,
          upper   = weights
        ),
        na.rm = TRUE
      )
    },
    future.seed = TRUE
  )

  Density / pracma::trapz(z_values, Density)
}
