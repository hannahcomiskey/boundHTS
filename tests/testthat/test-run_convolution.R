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

  expect_type(dens_post, "list")

  expect_true(all(
    c("Node", "Z", "Density") %in%
      names(dens_post[[1]])
  ))

  expect_equal(
    nrow(dens_post[[1]]),
    length(z_values)
  )

})

test_that("Beta convolution runs with point estimates", {

  set.seed(123)

  S_beta <- rbind(Total = c(1, 1, 1, 1),
                  A = c(1, 1, 0, 0),
                  B = c(0, 0, 1, 1),
                  A1 = c(1, 0, 0, 0),
                  A2 = c(0, 1, 0, 0),
                  B1 = c(0, 0, 1, 0),
                  B2 = c(0, 0, 0, 1))
  colnames(S_beta) <- c("A1", "A2", "B1", "B2")

  weights_list <- list(Total = c(A1 = 0.10, A2 = 0.15, B1 = 0.30, B2 = 0.45),
                       A = c(A1 = 0.40, A2 = 0.60),
                       B = c(B1 = 0.30, B2 = 0.70))

  alpha_mat <- c(Total = 8, A = 5, B = 6, A1 = 2, A2 = 4, B1 = 3, B2 = 5)

  beta_mat <- c(Total = 4, A = 3, B = 2, A1 = 6, A2 = 5, B1 = 7, B2 = 4)

  dens <- run_convolution(
    family = "Beta",
    S = S_beta,
    point = TRUE,
    z_values = NULL,
    alpha_mat = alpha_mat,
    beta_mat = beta_mat,
    weights_list = weights_list,
    n_draws = 500,
    n_sims = 50
  )
  expect_type(dens, "list")

})

test_that("ZOIB convolution runs with point estimates", {


  set.seed(123)

  S_beta <- rbind(Total = c(1, 1, 1, 1),
                  A = c(1, 1, 0, 0),
                  B = c(0, 0, 1, 1),
                  A1 = c(1, 0, 0, 0),
                  A2 = c(0, 1, 0, 0),
                  B1 = c(0, 0, 1, 0),
                  B2 = c(0, 0, 0, 1))
  colnames(S_beta) <- c("A1", "A2", "B1", "B2")

  weights_list <- list(Total = c(A1 = 0.10, A2 = 0.15, B1 = 0.30, B2 = 0.45),
                       A = c(A1 = 0.40, A2 = 0.60),
                       B = c(B1 = 0.30, B2 = 0.70))

  alpha_mat <- c(Total = 8, A = 5, B = 6, A1 = 2, A2 = 4, B1 = 3, B2 = 5)

  beta_mat <- c(Total = 4, A = 3, B = 2, A1 = 6, A2 = 5, B1 = 7, B2 = 4)

  zoi_mat <- c(Total = 0.01, A = 0.03, B = 0.02, A1 = 0.06, A2 = 0.05, B1 = 0.07, B2 = 0.04)

  coi_mat <- c(Total = 0.01, A = 0.01, B = 0.01, A1 = 0.01, A2 = 0.01, B1 = 0.01, B2 = 0.01)

  dens <- run_convolution(
    family = "ZOIB",
    S = S_beta,
    point = TRUE,
    z_values = NULL,
    alpha_mat = alpha_mat,
    beta_mat = beta_mat,
    zoi_mat = zoi_mat,
    coi_mat = coi_mat,
    weights_list = weights_list,
    n_draws = 500,
    n_sims = 50
  )
  expect_type(dens, "list")
})

test_that("run_convolution returns densities with expected columns", {

  set.seed(123)

  S_beta <- rbind(Total = c(1, 1, 1, 1),
                  A = c(1, 1, 0, 0),
                  B = c(0, 0, 1, 1),
                  A1 = c(1, 0, 0, 0),
                  A2 = c(0, 1, 0, 0),
                  B1 = c(0, 0, 1, 0),
                  B2 = c(0, 0, 0, 1))
  colnames(S_beta) <- c("A1", "A2", "B1", "B2")

  weights_list <- list(Total = c(A1 = 0.10, A2 = 0.15, B1 = 0.30, B2 = 0.45),
                       A = c(A1 = 0.40, A2 = 0.60),
                       B = c(B1 = 0.30, B2 = 0.70))

  alpha_mat <- c(Total = 8, A = 5, B = 6, A1 = 2, A2 = 4, B1 = 3, B2 = 5)

  beta_mat <- c(Total = 4, A = 3, B = 2, A1 = 6, A2 = 5, B1 = 7, B2 = 4)

  zoi_mat <- c(Total = 0.01, A = 0.03, B = 0.02, A1 = 0.06, A2 = 0.05, B1 = 0.07, B2 = 0.04)

  coi_mat <- c(Total = 0.01, A = 0.01, B = 0.01, A1 = 0.01, A2 = 0.01, B1 = 0.01, B2 = 0.01)

  dens <- run_convolution(
    family = "ZOIB",
    S = S_beta,
    point = TRUE,
    z_values = NULL,
    alpha_mat = alpha_mat,
    beta_mat = beta_mat,
    zoi_mat = zoi_mat,
    coi_mat = coi_mat,
    weights_list = weights_list,
    n_draws = 500,
    n_sims = 50
  )
  expect_type(dens, "list")

  expect_named(
    dens[[1]],
    c("Node", "Z", "Density")
  )

})
