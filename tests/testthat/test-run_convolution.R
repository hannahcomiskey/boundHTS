# tests/testthat/test-run_convolution.R
test_that("run_convolution function fails gracefully on invalid family", {
  set.seed(123)

  S <- rbind(Total = c(1, 1, 1, 1),
             A = c(1, 1, 0, 0),
             B = c(0, 0, 1, 1),
             A1 = c(1, 0, 0, 0),
             A2 = c(0, 1, 0, 0),
             B1 = c(0, 0, 1, 0),
             B2 = c(0, 0, 0, 1))
  colnames(S) <- c("A1", "A2", "B1", "B2")

  lambda_mat <- c(Total = 8, A = 5, B = 6, A1 = 2, A2 = 4, B1 = 3, B2 = 5)

  z_values <- seq(0, 100)

  expect_error(run_convolution(family = "Gaussian",
                               S = S,
                               lambda_mat = lambda_mat,
                               z_values = z_values,
                               point = TRUE,
                               n_sims = 50, n_draws = 100), "not")
})


test_that("Poisson convolution runs with point estimates", {

  S <- rbind(Total = c(1, 1, 1, 1),
             A = c(1, 1, 0, 0),
             B = c(0, 0, 1, 1),
             A1 = c(1, 0, 0, 0),
             A2 = c(0, 1, 0, 0),
             B1 = c(0, 0, 1, 0),
             B2 = c(0, 0, 0, 1))
  colnames(S) <- c("A1", "A2", "B1", "B2")

  lambda_mat <- c(Total = 8, A = 5, B = 6, A1 = 2, A2 = 4, B1 = 3, B2 = 5)

  z_values <- seq(0, 100)

  result <- run_convolution(family = "Poisson",
                            S = S,
                            lambda_mat = lambda_mat,
                            z_values = z_values,
                            point = TRUE,
                            n_sims = 50, n_draws = 100)

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
  S <- rbind(Total = c(1, 1, 1, 1),
             A = c(1, 1, 0, 0),
             B = c(0, 0, 1, 1),
             A1 = c(1, 0, 0, 0),
             A2 = c(0, 1, 0, 0),
             B1 = c(0, 0, 1, 0),
             B2 = c(0, 0, 0, 1))
  colnames(S) <- c("A1", "A2", "B1", "B2")

  lambda_post <- cbind(Total = rgamma(J, 8, 1),
                       A = rgamma(J, 5, 1),
                       B = rgamma(J, 6, 1),
                       A1 = rgamma(J, 2, 1),
                       A2 = rgamma(J, 4, 1),
                       B1 = rgamma(J, 3, 1),
                       B2 = rgamma(J, 5, 1))

  dens_post <- run_convolution(family = "Poisson",
                               S = S,
                               lambda_mat = lambda_post,
                               point = FALSE,
                               n_sims = 50,
                               n_draws = NULL,
                               z_values = seq(0, 100))

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
