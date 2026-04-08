test_that("ZIB convolution_density returns a finite scalar", {

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

  zi_matrix <- matrix(
    runif(n_draws * N, 0.1, 0.4),
    nrow = n_draws,
    ncol = N
  )

  weighted_samps <- array(
    runif(n_sims * n_draws * N, min = 0, max = 0.2),
    dim = c(n_sims, n_draws, N)
  )

  weights <- runif(N, min = 0.5, max = 2)

  # run function
  dens <- ZIB_convolution(z_values =  seq(0,1, length.out = 10),
                          alpha_input = alpha_matrix,
                          beta_input = beta_matrix,
                          zi_input = zi_matrix,
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

test_that("ZIB convolution returns a valid density over z_values", {

  set.seed(123)

  n_sims  <- 30
  n_draws <- 10
  N <- 3
  n_years <- 1

  z_values <- seq(0, 1, length.out = 101)

  # Monte Carlo draws of Y
  weighted_samps <- array(
    runif(n_sims * n_draws * N, min = 0, max = 0.2),
    dim = c(n_sims, n_draws, N)
  )

  alpha_matrix <- matrix(
    runif(n_draws * N, min = 2, max = 10),
    nrow = n_draws,
    ncol = N
  )

  beta_matrix <- matrix(
    runif(n_draws * N, min = 2, max = 10),
    nrow = n_draws,
    ncol = N
  )

  zi_matrix <- matrix(
    runif(n_draws * N, 0.1, 0.4),
    nrow = n_draws,
    ncol = N
  )

  weights <- rep(1, N)

  dens <- ZIB_convolution(z_values = z_values,
                           alpha_input = alpha_matrix,
                           beta_input = beta_matrix,
                           zi_input = zoi_matrix,
                           weighted_samps = weighted_samps,
                           weights = weights,
                           point=FALSE)

  expect_type(dens, "double")
  expect_length(dens, length(z_values))
  expect_true(all(is.finite(dens)))
  expect_true(all(dens >= 0))
})

test_that("ZIB convolution integrates to one", {

  set.seed(1)

  n_sims  <- 30
  n_draws <- 10
  N <- 3

  z_values <- seq(0, 1, length.out = 101)

  # Monte Carlo draws of Y
  weighted_samps <- array(
    runif(n_sims * n_draws * N, min = 0, max = 0.2),
    dim = c(n_sims, n_draws, N)
  )

  alpha_matrix <- matrix(
    runif(n_draws * N, min = 2, max = 10),
    nrow = n_draws,
    ncol = N
  )

  beta_matrix <- matrix(
    runif(n_draws * N, min = 2, max = 10),
    nrow = n_draws,
    ncol = N
  )

  zi_matrix <- matrix(
    runif(n_draws * N, 0.1, 0.4),
    nrow = n_draws,
    ncol = N
  )

  weights <- rep(1, N)

  dens <- ZIB_convolution(z_values = z_values,
                           alpha_input = alpha_matrix,
                           beta_input = beta_matrix,
                           zi_input = zi_matrix,
                           weights = weights,
                           weighted_samps = weighted_samps,
                           point=FALSE)

  integral <- pracma::trapz(z_values, dens)

  expect_equal(integral, 1, tolerance = 1e-6)
})

test_that("ZIB convolution is deterministic given fixed seed", {

  set.seed(42)

  n_sims  <- 30
  n_draws <- 10
  N <- 3

  z_values <- seq(0, 1, length.out = 101)

  # Monte Carlo draws of Y
  weighted_samps <- array(
    runif(n_sims * n_draws * N, min = 0, max = 0.2),
    dim = c(n_sims, n_draws, N)
  )

  alpha_matrix <- matrix(
    runif(n_draws * N, min = 2, max = 10),
    nrow = n_draws,
    ncol = N
  )

  beta_matrix <- matrix(
    runif(n_draws * N, min = 2, max = 10),
    nrow = n_draws,
    ncol = N
  )

  zi_matrix <- matrix(
    runif(n_draws * N, 0.1, 0.4),
    nrow = n_draws,
    ncol = N
  )

  weights <- rep(1, N)

  set.seed(999)

  dens1 <- ZIB_convolution(z_values = z_values,
                            alpha_input = alpha_matrix,
                            beta_input = beta_matrix,
                            zi_input = zi_matrix,
                            weights = weights,
                            weighted_samps = weighted_samps,
                            point=FALSE)

  set.seed(999)

  dens2 <- ZIB_convolution(z_values = z_values,
                          alpha_input = alpha_matrix,
                          beta_input = beta_matrix,
                          zi_input = zi_matrix,
                          weights = weights,
                          weighted_samps = weighted_samps,
                          point=FALSE)
  expect_equal(dens1, dens2)
})

test_that("ZIB convolution returns zero density outside support", {

  set.seed(7)

  n_sims  <- 20
  n_draws <- 5
  N <- 3

  z_values <- c(-1, 0, 0.5, 1, 2)

  # Monte Carlo draws of Y
  weighted_samps <- array(
    runif(n_sims * n_draws * N, min = 0, max = 0.2),
    dim = c(n_sims, n_draws, N)
  )

  alpha_matrix <- matrix(
    runif(n_draws * N, min = 2, max = 10),
    nrow = n_draws,
    ncol = N
  )

  beta_matrix <- matrix(
    runif(n_draws * N, min = 2, max = 10),
    nrow = n_draws,
    ncol = N
  )

  zi_matrix <- matrix(
    runif(n_draws * N, 0.1, 0.4),
    nrow = n_draws,
    ncol = N
  )

  weights <- rep(1, N)

  dens <- ZIB_convolution(z_values =  z_values,
                           alpha_input = alpha_matrix,
                           beta_input = beta_matrix,
                           zi_input = zi_matrix,
                           weighted_samps = weighted_samps,
                           weights = weights,
                           point=FALSE)

  expect_true(all(dens[z_values < 0] == 0))
  expect_true(all(dens[z_values > sum(weights)] == 0))
})

