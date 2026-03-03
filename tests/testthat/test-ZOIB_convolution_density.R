test_that("ZOIB_convolution_density returns a valid density over z_values", {

  set.seed(123)

  n_mc    <- 30
  n_draws <- 10
  n_nodes <- 2
  n_years <- 1

  z_values <- seq(0, 1, length.out = 101)

  # Monte Carlo draws of Y
  Y_mc <- matrix(
    runif(n_draws * n_nodes, 0.05, 0.3),
    nrow = n_draws,
    ncol = n_nodes
  )

  phi_array <- array(
    runif(n_draws * n_nodes, 5, 15),
    dim = c(n_draws, n_nodes)
  )

  zoi_array <- array(
    runif(n_draws * n_nodes, 0.1, 0.4),
    dim = c(n_draws, n_nodes)
  )

  coi_array <- array(
    runif(n_draws * n_nodes, 0.2, 0.8),
    dim = c(n_draws, n_nodes)
  )

  weights <- rep(1, n_nodes)

  dens <- ZOIB_convolution_density(
    Y_mc = Y_mc,
    phi_array = phi_array,
    zoi_array = zoi_array,
    coi_array = coi_array,
    weights = weights,
    z_values = z_values,
    n_mc = n_mc
  )

  expect_type(dens, "double")
  expect_length(dens, length(z_values))
  expect_true(all(is.finite(dens)))
  expect_true(all(dens >= 0))
})

test_that("ZOIB_convolution_density integrates to one", {

  set.seed(1)

  n_mc    <- 40
  n_draws <- 8
  n_nodes <- 2

  z_values <- seq(0, 1, length.out = 201)

  Y_mc <- matrix(
    runif(n_mc * n_nodes, 0.05, 0.25),
    n_mc, n_nodes
  )

  phi_array <- array(10, dim = c(n_mc, n_nodes))
  zoi_array <- array(0.2, dim = c(n_mc, n_nodes))
  coi_array <- array(0.5, dim = c(n_mc, n_nodes))

  weights <- rep(1, n_nodes)

  dens <- ZOIB_convolution_density(
    Y_mc = Y_mc,
    phi_array = phi_array,
    zoi_array = zoi_array,
    coi_array = coi_array,
    weights = weights,
    z_values = z_values,
    n_mc = n_mc
  )

  integral <- pracma::trapz(z_values, dens)

  expect_equal(integral, 1, tolerance = 1e-6)
})

test_that("ZOIB_convolution_density is deterministic given fixed seed", {

  set.seed(42)

  n_mc    <- 25
  n_draws <- 6
  n_nodes <- 2

  z_values <- seq(0, 1, length.out = 51)

  Y_mc <- matrix(
    runif(n_mc * n_nodes, 0.05, 0.3),
    n_draws, n_nodes
  )

  phi_array <- array(
    runif(n_draws * n_nodes, 5, 10),
    dim = c(n_draws, n_nodes)
  )

  zoi_array <- array(
    runif(n_draws * n_nodes, 0.1, 0.3),
    dim = c(n_draws, n_nodes)
  )

  coi_array <- array(
    runif(n_draws * n_nodes, 0.3, 0.7),
    dim = c(n_draws, n_nodes)
  )

  weights <- rep(1, n_nodes)

  set.seed(999)
  dens1 <- ZOIB_convolution_density(
    Y_mc = Y_mc,
    phi_array = phi_array,
    zoi_array = zoi_array,
    coi_array = coi_array,
    weights = weights,
    z_values = z_values,
    n_mc = n_mc
  )

  set.seed(999)
  dens2 <- ZOIB_convolution_density(
    Y_mc = Y_mc,
    phi_array = phi_array,
    zoi_array = zoi_array,
    coi_array = coi_array,
    weights = weights,
    z_values = z_values,
    n_mc = n_mc
  )

  expect_equal(dens1, dens2)
})

test_that("ZOIB_convolution_density returns zero density outside support", {

  set.seed(7)

  n_mc    <- 20
  n_draws <- 5
  n_nodes <- 2

  z_values <- c(-1, 0, 0.5, 1, 2)

  Y_mc <- matrix(
    runif(n_mc * n_nodes, 0.1, 0.3),
    n_draws, n_nodes
  )

  phi_array <- array(8, dim = c(n_draws, n_nodes))
  zoi_array <- array(0.2, dim = c(n_draws, n_nodes))
  coi_array <- array(0.5, dim = c(n_draws, n_nodes))

  weights <- rep(1, n_nodes)

  dens <- ZOIB_convolution_density(
    Y_mc = Y_mc,
    phi_array = phi_array,
    zoi_array = zoi_array,
    coi_array = coi_array,
    weights = weights,
    z_values = z_values,
    n_mc = n_mc
  )

  expect_true(all(dens[z_values < 0] == 0))
  expect_true(all(dens[z_values > sum(weights)] == 0))
})

