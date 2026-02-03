test_that("tilted_density_cont returns a numeric vector of correct length", {

  y_vals <- seq(-5, 5, length.out = 1001)
  f_y <- dnorm(y_vals)

  nu <- 0.3

  f_tilt <- tilted_density_cont(
    nu = nu,
    f_y = f_y,
    y_vals = y_vals
  )

  expect_type(f_tilt, "double")
  expect_length(f_tilt, length(y_vals))
  expect_true(all(is.finite(f_tilt)))
  expect_true(all(f_tilt >= 0))
})

test_that("tilted_density_cont integrates to one", {

  y_vals <- seq(-6, 6, length.out = 2001)
  f_y <- dnorm(y_vals)

  nu <- -0.4

  f_tilt <- tilted_density_cont(
    nu = nu,
    f_y = f_y,
    y_vals = y_vals
  )

  integral <- pracma::trapz(y_vals, f_tilt)

  expect_equal(integral, 1, tolerance = 1e-6)
})

test_that("tilted_density_cont is identity when nu = 0", {

  y_vals <- seq(-5, 5, length.out = 1501)
  f_y <- dnorm(y_vals)

  f_tilt <- tilted_density_cont(
    nu = 0,
    f_y = f_y,
    y_vals = y_vals
  )

  # should recover original density
  expect_equal(
    f_tilt,
    f_y / pracma::trapz(y_vals, f_y),
    tolerance = 1e-10
  )
})

test_that("tilted_density_cont shifts the mean correctly for Normal density", {

  # For N(0,1), exponential tilting by nu shifts mean to nu
  y_vals <- seq(-7, 7, length.out = 4001)
  f_y <- dnorm(y_vals)

  nu <- 0.6

  f_tilt <- tilted_density_cont(
    nu = nu,
    f_y = f_y,
    y_vals = y_vals
  )

  mean_tilt <- pracma::trapz(y_vals, y_vals * f_tilt)

  expect_equal(mean_tilt, nu, tolerance = 5e-3)
})

test_that("tilted_density_cont remains finite for small base densities", {

  y_vals <- seq(-12, 12, length.out = 5001)
  f_y <- dnorm(y_vals, sd = 2)

  nu <- -0.8

  f_tilt <- tilted_density_cont(
    nu = nu,
    f_y = f_y,
    y_vals = y_vals
  )

  expect_true(all(is.finite(f_tilt)))
  expect_true(all(f_tilt >= 0))
})
