test_that("run_convolution_with_tilting returns correctly structured output", {
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

   res <- run_convolution_with_tilting(family = "Poisson",
                                       S = S,
                                       mu_theory = c(Total = 9, A = 4, B = 10, A1 = 1, A2 = 4, B1 = 2, B2 = 7),
                                       lambda_mat = lambda_mat,
                                       z_values = z_values,
                                       point = TRUE,
                                       n_sims = 50,
                                       n_draws = 100)


  ## ---- structure checks ----
  expect_type(res, "list")

  expect_true(all(c(
    "base_density",
    "tilted_density",
    "tilted_samples",
    "tilting_parameter"
  ) %in% names(res)))

  expect_length(res$base_density, 3)
  expect_length(res$tilted_density, 3)

  ## ---- no NULL outputs ----
  expect_false(any(vapply(res$base_density, is.null, logical(1))))
  expect_false(any(vapply(res$tilted_density, is.null, logical(1))))
  expect_false(any(vapply(res$tilted_samples, is.null, logical(1))))

  ## ---- tilting parameters exist ----
  expect_true(all(vapply(res$tilting_parameter, is.numeric, logical(1))))
})


test_that("convolution and tilting preserve grid consistency", {

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

  res <- run_convolution_with_tilting(family = "Poisson",
                                      S = S,
                                      mu_theory = c(Total = 9, A = 4, B = 10, A1 = 1, A2 = 4, B1 = 2, B2 = 7),
                                      lambda_mat = lambda_mat,
                                      z_values = z_values,
                                      point = TRUE,
                                      n_sims = 50,
                                      n_draws = 100)

  expect_true(all(vapply(res$base_density, function(df) {
    identical(df$Z, z_values)
  }, logical(1))))

  expect_true(all(vapply(res$tilted_density, function(df) {
    identical(df$Z, z_values)
  }, logical(1))))
})


test_that("density outputs are valid probability distributions", {

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

  res <- run_convolution_with_tilting(family = "Poisson",
                                      S = S,
                                      mu_theory = c(Total = 9, A = 4, B = 10, A1 = 1, A2 = 4, B1 = 2, B2 = 7),
                                      lambda_mat = lambda_mat,
                                      z_values = z_values,
                                      point = TRUE,
                                      n_sims = 50,
                                      n_draws = 100)

  ## ---- non-negativity ----
  expect_true(all(vapply(res$base_density, function(df) {
    all(df$Density >= 0)
  }, logical(1))))

  expect_true(all(vapply(res$tilted_density, function(df) {
    all(df$Density >= 0)
  }, logical(1))))

  ## ---- finiteness ----
  expect_true(all(vapply(res$base_density, function(df) {
    all(is.finite(df$Density))
  }, logical(1))))
})


test_that("tilting meaningfully changes distributions", {
  set.seed(123)

  S <- rbind(Total = c(1, 1, 1, 1),
             A = c(1, 1, 0, 0),
             B = c(0, 0, 1, 1),
             A1 = c(1, 0, 0, 0),
             A2 = c(0, 1, 0, 0),
             B1 = c(0, 0, 1, 0),
             B2 = c(0, 0, 0, 1))
  colnames(S) <- c("A1", "A2", "B1", "B2")

  lambda_mat <- c(Total = 15, A = 7, B = 8, A1 = 3, A2 = 4, B1 = 3, B2 = 5)

  z_values <- seq(0, 100)

  res <- run_convolution_with_tilting(family = "Poisson",
                                      S = S,
                                      mu_theory = c(Total = 20, A = 10, B = 10, A1 = 4, A2 = 6, B1 = 2, B2 = 8),
                                      lambda_mat = lambda_mat,
                                      z_values = z_values,
                                      point = TRUE,
                                      n_sims = 50,
                                      n_draws = 100)

  diffs <- mapply(function(base, tilt) {
    sum(abs(base$Density - tilt$Density))
  }, res$base_density, res$tilted_density)

  expect_true(any(diffs > 1e-6))
})


test_that("tilted samples are valid numeric vectors", {

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

  res <- run_convolution_with_tilting(family = "Poisson",
                                      S = S,
                                      mu_theory = c(Total = 9, A = 4, B = 10, A1 = 1, A2 = 4, B1 = 2, B2 = 7),
                                      lambda_mat = lambda_mat,
                                      z_values = z_values,
                                      point = TRUE,
                                      n_sims = 50,
                                      n_draws = 100)

  expect_true(all(vapply(res$tilted_samples, is.numeric, logical(1))))
  expect_true(all(vapply(res$tilted_samples, length, integer(1)) > 0))
})


test_that("run_convolution_with_tilting function fails gracefully on invalid family", {
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

  expect_error(run_convolution_with_tilting(family = "Gaussian",
                                            S = S,
                                            mu_theory = c(Total = 9, A = 4, B = 10, A1 = 1, A2 = 4, B1 = 2, B2 = 7),
                                            lambda_mat = lambda_mat,
                                            z_values = z_values,
                                            point = TRUE,
                                            n_sims = 50,
                                            n_draws = 100), "not")
})
