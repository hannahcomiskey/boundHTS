test_that("ZIB_convolution_density returns a finite scalar", {

  set.seed(123)

  S <- rbind(
    Total = c(1, 1, 1, 1),
    A     = c(1, 1, 0, 0),
    B     = c(0, 0, 1, 1),
    A1    = c(1, 0, 0, 0),
    A2    = c(0, 1, 0, 0),
    B1    = c(0, 0, 1, 0),
    B2    = c(0, 0, 0, 1)
  )
  colnames(S) <- c("A1", "A2", "B1", "B2")

  local_weights <- list(
    Total = c(A1 = 0.10, A2 = 0.15, B1 = 0.30, B2 = 0.45),
    A     = c(A1 = 0.40, A2 = 0.60),
    B     = c(B1 = 0.30, B2 = 0.70)
  )

  alpha_mat <- c(
    Total = 8,
    A = 5,
    B = 6,
    A1 = 2,
    A2 = 4,
    B1 = 3,
    B2 = 5
  )

  beta_mat <- c(
    Total = 4,
    A = 3,
    B = 2,
    A1 = 6,
    A2 = 5,
    B1 = 7,
    B2 = 4
  )

  zi_mat <- c(
    Total = 0.01,
    A = 0.05,
    B = 0.02,
    A1 = 0.06,
    A2 = 0.05,
    B1 = 0.01,
    B2 = 0.02
  )

  dens <- ZIB_convolution(S = S,
                          alpha_mat = alpha_mat,
                          beta_mat = beta_mat,
                          zi_mat = zi_mat,
                          weights_list = local_weights,
                          point = TRUE,
                          n_draws = 500,
                          n_sims = 50)

  dens

  # ---- expectations ----
  expect_type(dens, "list")
  expect_length(dens[[1]]$Density, 1000)
  expect_true(all(is.finite(dens[[1]]$Density)))
  expect_true(all(dens[[1]]$Density >= 0))
})

test_that("ZIB_convolution_density is stable when z is out of support", {

  set.seed(123)

  S <- rbind(
    Total = c(1, 1, 1, 1),
    A     = c(1, 1, 0, 0),
    B     = c(0, 0, 1, 1),
    A1    = c(1, 0, 0, 0),
    A2    = c(0, 1, 0, 0),
    B1    = c(0, 0, 1, 0),
    B2    = c(0, 0, 0, 1)
  )
  colnames(S) <- c("A1", "A2", "B1", "B2")

  local_weights <- list(
    Total = c(A1 = 0.10, A2 = 0.15, B1 = 0.30, B2 = 0.45),
    A     = c(A1 = 0.40, A2 = 0.60),
    B     = c(B1 = 0.30, B2 = 0.70)
  )

  alpha_mat <- c(
    Total = 8,
    A = 5,
    B = 6,
    A1 = 2,
    A2 = 4,
    B1 = 3,
    B2 = 5
  )

  beta_mat <- c(
    Total = 4,
    A = 3,
    B = 2,
    A1 = 6,
    A2 = 5,
    B1 = 7,
    B2 = 4
  )

  zi_mat <- c(
    Total = 0.01,
    A = 0.05,
    B = 0.02,
    A1 = 0.06,
    A2 = 0.05,
    B1 = 0.01,
    B2 = 0.02
  )


  ## Estimate densities
  expect_error(ZIB_convolution(S = S,
                               alpha_mat = alpha_mat,
                               beta_mat = beta_mat,
                               zi_mat = zi_mat,
                               weights_list = local_weights,
                               z_values = seq(-5, 5, length.out = 200),
                               point = TRUE,
                               n_draws = 500,
                               n_sims = 50))

})

test_that("ZIB_convolution_density errors on incompatible dimensions", {

  set.seed(123)

  S <- rbind(
    Total = c(1, 1, 1, 1),
    A     = c(1, 1, 0, 0),
    B     = c(0, 0, 1, 1),
    A1    = c(1, 0, 0, 0),
    A2    = c(0, 1, 0, 0),
    B1    = c(0, 0, 1, 0),
    B2    = c(0, 0, 0, 1)
  )
  colnames(S) <- c("A1", "A2", "B1", "B2")

  local_weights <- list(
    Total = c(A1 = 0.10, A2 = 0.15, B1 = 0.30, B2 = 0.45),
    A     = c(A1 = 0.40, A2 = 0.60),
    B     = c(B1 = 0.30, B2 = 0.70)
  )

  alpha_mat <- c(
    Total = 8,
    A = 5,
    B = 6,
    C = 9,
    A1 = 2,
    A2 = 4,
    B1 = 3,
    B2 = 5
  )

  beta_mat <- c(
    Total = 4,
    A = 3,
    B = 2,
    A1 = 6,
    A2 = 5,
    B1 = 7,
    B2 = 4
  )

  zi_mat <- c(
    Total = 0.01,
    A = 0.05,
    B = 0.02,
    A1 = 0.06,
    A2 = 0.05,
    B1 = 0.01,
    B2 = 0.02
  )


  ## Estimate densities
  expect_error(ZIB_convolution(S = S,
                               alpha_mat = alpha_mat,
                               beta_mat = beta_mat,
                               zi_mat = zi_mat,
                               weights_list = local_weights,
                               z_values = seq(-5, 5, length.out = 200),
                               point = TRUE,
                               n_draws = 500,
                               n_sims = 50))

})
