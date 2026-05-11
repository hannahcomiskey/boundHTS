# tests/testthat/test-run_convolution.R
library(testthat)
devtools::load_all()

test_that("run_convolution errors for unsupported family", {

  expect_error(
    run_convolution(
      family = "Gaussian",
      groups = list(A = c("a1", "a2")),
      point = TRUE,
      z_values = 0:10
    ),
    "Family name not recognised. Please refer to help file for naming."
  )

})

test_that("Poisson convolution runs with point estimates", {

  groups <- list(
    NSW = c("NSW_Male", "NSW_Female"),
    National = c("AUS")
  )

  lambda_list <- list(c(12, 10), 22)

  z_values <- 0:40

  result <- run_convolution(
    family = "Poisson",
    groups = groups,
    point = TRUE,
    z_values = z_values,
    lambda_list = lambda_list
  )

  expect_type(result, "list")

  expect_true(all(
    c("Node", "Z", "Density") %in%
      names(result[[1]])
  ))

  expect_equal(
    nrow(result[[1]]),
    length(z_values)
  )

})

test_that("Poisson convolution runs with posterior samples", {

  set.seed(123)

  J <- 100

  groups <- list(
    NSW = c("NSW_Male", "NSW_Female"),
    National = c("AUS")
  )

  lambda_list <- list(
    cbind(
      rgamma(J, shape = 12, rate = 1),
      rgamma(J, shape = 10, rate = 1)
    ),
    matrix(
      rgamma(J, shape = 22, rate = 1),
      ncol = 1
    )
  )

  z_values <- 0:40

  result <- run_convolution(
    family = "Poisson",
    groups = groups,
    point = FALSE,
    z_values = z_values,
    lambda_list = lambda_list
  )

  expect_type(result, "list")

  expect_true(all(
    c("Node", "Z", "Density") %in%
      names(result[[1]])
  ))

  expect_equal(
    nrow(result[[1]]),
    length(z_values)
  )

})

test_that("Beta convolution runs with point estimates", {

  groups <- list(
    State = c("NSW", "VIC"),
    National = c("AUS")
  )

  alpha_list <- list(
    c(2, 5),
    7
  )

  beta_list <- list(
    c(6, 3),
    4
  )

  weights_list <- list(
    c(0.4, 0.6),
    1
  )

  result <- run_convolution(
    family = "Beta",
    groups = groups,
    point = TRUE,
    z_values = NULL,
    alpha_list = alpha_list,
    beta_list = beta_list,
    weights_list = weights_list
  )

  expect_type(result, "list")

})

test_that("ZOIB convolution runs with point estimates", {

  groups <- list(
    State = c("NSW", "VIC"),
    National = c("AUS")
  )

  alpha_list <- list(
    c(2, 5),
    7
  )

  beta_list <- list(
    c(6, 3),
    4
  )

  zoi_list <- list(
    c(0.1, 0.2),
    0.05
  )

  coi_list <- list(
    c(0.05, 0.10),
    0.02
  )

  weights_list <- list(
    c(0.4, 0.6),
    1
  )

  result <- run_convolution(
    family = "ZOIB",
    groups = groups,
    point = TRUE,
    z_values = NULL,
    alpha_list = alpha_list,
    beta_list = beta_list,
    zoi_list = zoi_list,
    coi_list = coi_list,
    weights_list = weights_list
  )

  expect_type(result, "list")

})

test_that("run_convolution returns densities with expected columns", {

  groups <- list(
    NSW = c("NSW_Male", "NSW_Female"),
    National = c("AUS")
  )

  lambda_list <- list(
    c(12, 10),
    22
  )

  z_values <- 0:20

  result <- run_convolution(
    family = "Poisson",
    groups = groups,
    point = TRUE,
    z_values = z_values,
    lambda_list = lambda_list
  )

  expect_named(
    result[[1]],
    c("Node", "Z", "Density")
  )

})
