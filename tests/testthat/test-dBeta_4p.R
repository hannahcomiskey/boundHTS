test_that("dBeta_4p returns valid density within support", {

  x     <- 0.25
  lower <- 0
  upper <- 0.5
  alpha <- 5
  beta  <- 2

  dens <- dBeta_4p(
    x = x,
    lower = lower,
    upper = upper,
    alpha = alpha,
    beta = beta
  )

  expect_type(dens, "double")
  expect_length(dens, 1L)
  expect_true(is.finite(dens))
  expect_true(dens > 0)
})

test_that("dBeta_4p returns zero outside support", {

  lower <- 0
  upper <- 0.5
  alpha <- 3
  beta  <- 3

  expect_equal(dBeta_4p(x = -0.1, lower = lower, upper = upper, alpha = alpha, beta = beta), 0)

  expect_equal(dBeta_4p(x = 0.6, lower = lower, upper = upper, alpha = alpha, beta = beta), 0)
})

test_that("dBeta_4p returns zero for invalid shape parameters", {

  x <- 0.2
  lower <- 0
  upper <- 1

  expect_equal(dBeta_4p(x = x, lower = lower, upper = upper, alpha = -1, beta = 2), 0)

  expect_equal(dBeta_4p(x = x, lower = lower, upper = upper, alpha = 2, beta = -1), 0)
})

test_that("dBeta_4p matches ExtDist::dBeta_ab for valid inputs", {

  x     <- 0.3
  lower <- 0
  upper <- 1
  alpha <- 4
  beta  <- 6

  expect_equal(dBeta_4p(x, lower, upper, alpha, beta), ExtDist::dBeta_ab(x, alpha, beta, lower, upper), tolerance = 1e-12)
})
