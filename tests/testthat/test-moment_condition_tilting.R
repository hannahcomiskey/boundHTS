test_that("moment_condition_tilting returns a scalar numeric residual", {

  y_vals <- seq(-5, 5, length.out = 1001)
  f_y <- dnorm(y_vals, mean = 0, sd = 1)

  nu <- 0.3
  mu_theory <- 1

  res <- moment_condition_tilting(
    nu = nu,
    f_y = f_y,
    y_vals = y_vals,
    mu_theory = mu_theory
  )

  expect_type(res, "double")
  expect_length(res, 1L)
  expect_true(is.finite(res))
})

test_that("moment_condition_tilting is zero at nu = 0 when target mean matches density mean", {

  y_vals <- seq(-5, 5, length.out = 2001)
  f_y <- dnorm(y_vals, mean = 0, sd = 1)

  mu_theory <- 0

  res <- moment_condition_tilting(
    nu = 0,
    f_y = f_y,
    y_vals = y_vals,
    mu_theory = mu_theory
  )

  expect_equal(res, 0, tolerance = 1e-6)
})

test_that("moment_condition_tilting recovers known mean shift under exponential tilting", {

  # For a standard normal, exponential tilting by nu shifts the mean to nu
  y_vals <- seq(-6, 6, length.out = 3001)
  f_y <- dnorm(y_vals, mean = 0, sd = 1)

  nu <- 0.5
  mu_theory <- nu

  res <- moment_condition_tilting(
    nu = nu,
    f_y = f_y,
    y_vals = y_vals,
    mu_theory = mu_theory
  )

  expect_equal(res, 0, tolerance = 5e-4)
})

test_that("moment_condition_tilting returns positive residual when tilted mean exceeds target", {

  y_vals <- seq(-5, 5, length.out = 1501)
  f_y <- dnorm(y_vals)

  nu <- 0.8
  mu_theory <- 0

  res <- moment_condition_tilting(
    nu = nu,
    f_y = f_y,
    y_vals = y_vals,
    mu_theory = mu_theory
  )

  expect_true(res > 0)
})

test_that("moment_condition_tilting is stable for small densities", {

  y_vals <- seq(-10, 10, length.out = 4001)
  f_y <- dnorm(y_vals, mean = 0, sd = 2)

  nu <- -0.7
  mu_theory <- -0.5

  res <- moment_condition_tilting(
    nu = nu,
    f_y = f_y,
    y_vals = y_vals,
    mu_theory = mu_theory
  )

  expect_true(is.finite(res))
})
