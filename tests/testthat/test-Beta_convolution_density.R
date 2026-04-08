test_that("Beta_convolution_density returns a finite scalar", {

  set.seed(123)

  n_sims  <- 30
  n_draws <- 10
  N       <- 3

  # inputs
  alpha_matrix <- matrix(
    runif(n_sims * N, min = 2, max = 10),
    nrow = n_sims,
    ncol = N
  )

  beta_matrix <- matrix(
    runif(n_sims * N, min = 2, max = 10),
    nrow = n_sims,
    ncol = N
  )

  weighted_samps <- array(
    runif(n_sims * n_draws * N, min = 0, max = 0.2),
    dim = c(n_sims, n_draws, N)
  )

  weights <- runif(N, min = 0.5, max = 2)

  # run function
  dens <- Beta_convolution(z_values =  seq(0,1, length.out = 10),
                           alpha_input = alpha_matrix,
                           beta_input = beta_matrix,
                           weighted_samps = weighted_samps,
                           weights = weights,
                           point=FALSE)
  dens

  # ---- expectations ----
  expect_type(dens, "double")
  expect_length(dens, 10)
  expect_true(all(is.finite(dens)))
  expect_true(all(dens >= 0))
})

test_that("Beta_convolution_density is stable when z is out of support", {

  set.seed(1)

  n_sims  <- 10
  n_draws <- 5
  N       <- 2

  z_values <- seq(-10, -1, length.out = 10) # well outside support

  alpha_matrix <- matrix(5, n_sims, N)
  beta_matrix  <- matrix(5, n_sims, N)

  weighted_samps <- array(
    runif(n_sims * n_draws * N, 0, 0.1),
    dim = c(n_sims, n_draws, N)
  )

  weights <- rep(1, N)

  dens <- Beta_convolution(z_values =  z_values,
                           alpha_input = alpha_matrix,
                           beta_input = beta_matrix,
                           weighted_samps = weighted_samps,
                           weights = weights,
                           point=FALSE)

  expect_true(all(is.finite(dens)))
  expect_equal(dens[1], 0, tolerance = 1e-12)
})

test_that("Beta_convolution_density errors on incompatible dimensions", {

  alpha_matrix <- matrix(2, nrow = 5, ncol = 3)
  beta_matrix  <- matrix(2, nrow = 5, ncol = 3)

  # wrong third dimension (should be N = 3)
  weighted_samps <- array(runif(5 * 4 * 2), dim = c(5, 4, 2))

  weights <- c(1, 1, 1)

  expect_error(
    Beta_convolution(z_values = seq(0,1, length.out = 10),
                     alpha_input = alpha_matrix,
                     beta_input = beta_matrix,
                     weighted_samps = weighted_samps,
                     weights = weights,
                     point=FALSE
    ),
    regexp = NA
  )
})
