test_that("tilted_density_discrete returns a numeric probability vector", {

  y_vals <- -3:3
  f_y <- dpois(y_vals + 3, lambda = 3)  # simple discrete pmf

  nu <- 0.4

  f_tilt <- tilted_density_discrete(
    nu = nu,
    f_y = f_y,
    y_vals = y_vals
  )

  expect_type(f_tilt, "double")
  expect_length(f_tilt, length(y_vals))
  expect_true(all(is.finite(f_tilt)))
  expect_true(all(f_tilt >= 0))
})

test_that("tilted_density_discrete sums to one", {

  y_vals <- -5:5
  f_y <- dpois(y_vals + 5, lambda = 5)

  nu <- -0.6

  f_tilt <- tilted_density_discrete(
    nu = nu,
    f_y = f_y,
    y_vals = y_vals
  )

  expect_equal(sum(f_tilt), 1, tolerance = 1e-12)
})

test_that("tilted_density_discrete is identity when nu = 0", {

  y_vals <- 0:6
  f_y <- dpois(y_vals, lambda = 2)

  f_tilt <- tilted_density_discrete(
    nu = 0,
    f_y = f_y,
    y_vals = y_vals
  )

  expect_equal(
    f_tilt,
    f_y / sum(f_y),
    tolerance = 1e-12
  )
})

test_that("tilted_density_discrete shifts the mean in the correct direction", {

  y_vals <- 0:10
  f_y <- dpois(y_vals, lambda = 3)

  nu <- 0.5

  f_tilt <- tilted_density_discrete(
    nu = nu,
    f_y = f_y,
    y_vals = y_vals
  )

  mean_orig <- sum(y_vals * f_y) / sum(f_y)
  mean_tilt <- sum(y_vals * f_tilt)

  expect_true(mean_tilt > mean_orig)
})

test_that("tilted_density_discrete remains finite for small probabilities", {

  y_vals <- -10:10
  f_y <- dnorm(y_vals, mean = 0, sd = 3)
  f_y <- f_y / sum(f_y)  # discretised pmf

  nu <- -0.9

  f_tilt <- tilted_density_discrete(
    nu = nu,
    f_y = f_y,
    y_vals = y_vals
  )

  expect_true(all(is.finite(f_tilt)))
  expect_true(all(f_tilt >= 0))
})
