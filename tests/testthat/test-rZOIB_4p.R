test_that("rZOIB_4p returns array with correct dimensions and finite values", {

  set.seed(123)

  n_mc    <- 30
  n_years <- 4
  n_nodes <- 3
  n_draws <- 10

  sub_obs_data <- matrix(
    runif(n_years * n_nodes, 0.1, 0.8),
    nrow = n_years,
    ncol = n_nodes
  )

  phi_array <- array(
    runif(n_draws * n_nodes * n_years, 5, 20),
    dim = c(n_draws, n_nodes, n_years)
  )

  zoi_array <- array(
    runif(n_draws * n_nodes * n_years, 0.1, 0.4),
    dim = c(n_draws, n_nodes, n_years)
  )

  coi_array <- array(
    runif(n_draws * n_nodes * n_years, 0.2, 0.8),
    dim = c(n_draws, n_nodes, n_years)
  )

  weights <- runif(n_nodes, 0.5, 2)

  Y <- rZOIB_4p(
    n_mc = n_mc,
    sub_obs_data = sub_obs_data,
    phi_array = phi_array,
    zoi_array = zoi_array,
    coi_array = coi_array,
    weights = weights
  )

  expect_type(Y, "double")
  expect_true(is.array(Y))
  expect_equal(dim(Y), c(n_mc, n_nodes, n_years))
  expect_false(anyNA(Y))
  expect_true(all(is.finite(Y)))
  expect_true(all(Y >= 0))
  expect_true(all(Y <= max(weights)))
})

test_that("rZOIB_4p produces only zero/one inflation when zoi = 1", {

  set.seed(1)

  n_mc    <- 20
  n_years <- 3
  n_nodes <- 2
  n_draws <- 5

  sub_obs_data <- matrix(0.5, n_years, n_nodes)

  phi_array <- array(10, dim = c(n_draws, n_nodes, n_years))

  zoi_array <- array(1, dim = c(n_draws, n_nodes, n_years))  # always inflate
  coi_array <- array(0.7, dim = c(n_draws, n_nodes, n_years))

  weights <- c(1, 2)

  Y <- rZOIB_4p(
    n_mc = n_mc,
    sub_obs_data = sub_obs_data,
    phi_array = phi_array,
    zoi_array = zoi_array,
    coi_array = coi_array,
    weights = weights
  )

  # values must be exactly 0 or weights
  allowed <- c(0, weights)

  expect_true(all(Y %in% allowed))
})

test_that("rZOIB_4p reduces to Beta sampler when zoi = 0", {

  set.seed(42)

  n_mc    <- 25
  n_years <- 2
  n_nodes <- 3
  n_draws <- 6

  sub_obs_data <- matrix(
    runif(n_years * n_nodes, 0.2, 0.7),
    n_years, n_nodes
  )

  phi_array <- array(
    runif(n_draws * n_nodes * n_years, 8, 15),
    dim = c(n_draws, n_nodes, n_years)
  )

  zoi_array <- array(0, dim = c(n_draws, n_nodes, n_years))
  coi_array <- array(0.5, dim = c(n_draws, n_nodes, n_years)) # irrelevant

  weights <- rep(1, n_nodes)

  Y <- rZOIB_4p(
    n_mc = n_mc,
    sub_obs_data = sub_obs_data,
    phi_array = phi_array,
    zoi_array = zoi_array,
    coi_array = coi_array,
    weights = weights
  )

  # no point masses at 0 or 1 expected
  expect_false(any(Y == 0))
  expect_false(any(Y == 1))
})

test_that("rZOIB_4p works with single-node input", {

  set.seed(7)

  n_mc    <- 15
  n_years <- 4
  n_nodes <- 1
  n_draws <- 5

  sub_obs_data <- matrix(runif(n_years, 0.1, 0.9), ncol = 1)

  phi_array <- array(
    runif(n_draws * n_years, 5, 10),
    dim = c(n_draws, n_nodes, n_years)
  )

  zoi_array <- array(
    runif(n_draws * n_years, 0.2, 0.5),
    dim = c(n_draws, n_nodes, n_years)
  )

  coi_array <- array(
    runif(n_draws * n_years, 0.3, 0.6),
    dim = c(n_draws, n_nodes, n_years)
  )

  weights <- 1

  Y <- rZOIB_4p(
    n_mc = n_mc,
    sub_obs_data = sub_obs_data,
    phi_array = phi_array,
    zoi_array = zoi_array,
    coi_array = coi_array,
    weights = weights
  )

  expect_equal(dim(Y), c(n_mc, 1, n_years))
  expect_true(all(is.finite(Y)))
})
