test_that("rBeta_4p returns array with correct dimensions and support", {

  set.seed(123)

  # dimensions
  n_mc    <- 20
  n_years <- 5
  n_nodes <- 3
  n_draws <- 10

  # inputs
  sub_obs_data <- matrix(
    runif(n_years * n_nodes, min = 0.1, max = 0.9),
    nrow = n_years,
    ncol = n_nodes
  )

  phi_array <- array(
    runif(n_draws * n_nodes * n_years, min = 5, max = 20),
    dim = c(n_draws, n_nodes, n_years)
  )

  weights <- runif(n_nodes, min = 0.5, max = 2)

  # run sampler
  Y <- rBeta_4p(
    n_mc = n_mc,
    sub_obs_data = sub_obs_data,
    phi_array = phi_array,
    weights = weights
  )

  # ---- structure tests ----
  expect_type(Y, "double")
  expect_true(is.array(Y))
  expect_equal(dim(Y), c(n_mc, n_nodes, n_years))

  # ---- value tests ----
  expect_false(anyNA(Y))
  expect_true(all(Y >= 0))
  expect_true(all(Y <= max(weights)))

})

test_that("rBeta_4p errors on inconsistent dimensions", {

  n_mc <- 5

  sub_obs_data <- matrix(runif(6), nrow = 3, ncol = 2)
  phi_array    <- array(runif(4 * 3 * 3), dim = c(4, 3, 3))  # wrong n_nodes
  weights      <- c(1, 1)

  expect_error(rBeta_4p(n_mc, sub_obs_data, phi_array, weights))

})
