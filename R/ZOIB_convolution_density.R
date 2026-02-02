convolution_density <- function(Y_t, # [n_mc × n_nodes]
                                phi_array, # [n_draws × n_nodes]
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
          Y_mc    = Y_t,
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
