test_that("Poisson_convolution returns valid hierarchical densities", {

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

  ## ---------------------------------------------------------
  ## 1. Output is list
  ## ---------------------------------------------------------
  expect_type(res, "list")

  ## ---------------------------------------------------------
  ## 2. No NULL outputs
  ## ---------------------------------------------------------
  expect_true(all(!vapply(res, is.null, logical(1))))

  ## ---------------------------------------------------------
  ## 3. Required columns exist
  ## ---------------------------------------------------------
  expect_true(all(vapply(res$tilted_density, function(df) {
    all(c("Node", "Z", "Density") %in% names(df))
  }, logical(1))))

  ## ---------------------------------------------------------
  ## 5. No NA in Node column
  ## ---------------------------------------------------------
  expect_true(all(vapply(res$tilted_density, function(df) {
    !any(is.na(df$Node))
  }, logical(1))))

  ## ---------------------------------------------------------
  ## 6. Z grid consistency
  ## ---------------------------------------------------------
  expect_true(all(vapply(res$tilted_density, function(df) {
    identical(df$Z, z_values)
  }, logical(1))))

  ## ---------------------------------------------------------
  ## 7. Density sanity checks
  ## ---------------------------------------------------------
  expect_true(all(vapply(res$tilted_density, function(df) {
    all(is.finite(df$Density))
  }, logical(1))))

  expect_true(all(vapply(res$tilted_density, function(df) {
    all(df$Density >= 0)
  }, logical(1))))

  ## ---------------------------------------------------------
  ## 8. No completely empty outputs
  ## ---------------------------------------------------------
  expect_true(all(vapply(res$tilted_density, nrow, integer(1)) > 0))
})

