test_that("Poisson_convolution returns valid hierarchical densities", {

  set.seed(123)

  ## ---------------------------------------------------------
  ## Simple 2-level hierarchy
  ## ---------------------------------------------------------
  groups <- list(
    State = c("NSW_Male", "NSW_Female"),
    National = c("AUS")
  )

  z_values <- 0:20

  lambda_list <- list(
    c(10, 12),
    22
  )

  ## ---------------------------------------------------------
  ## Run function
  ## ---------------------------------------------------------
  res <- Poisson_convolution(
    groups = groups,
    z_values = z_values,
    lambda_list = lambda_list,
    point = TRUE
  )

  ## ---------------------------------------------------------
  ## 1. Output is list
  ## ---------------------------------------------------------
  expect_type(res, "list")

  ## ---------------------------------------------------------
  ## 2. No NULL outputs
  ## ---------------------------------------------------------
  expect_true(all(!vapply(res, is.null, logical(1))))

  ## ---------------------------------------------------------
  ## 3. Each element is a tibble/data.frame
  ## ---------------------------------------------------------
  expect_true(all(vapply(res, inherits, logical(1), "data.frame")))

  ## ---------------------------------------------------------
  ## 4. Required columns exist
  ## ---------------------------------------------------------
  expect_true(all(vapply(res, function(df) {
    all(c("Node", "Z", "Density") %in% names(df))
  }, logical(1))))

  ## ---------------------------------------------------------
  ## 5. No NA in Node column
  ## ---------------------------------------------------------
  expect_true(all(vapply(res, function(df) {
    !any(is.na(df$Node))
  }, logical(1))))

  ## ---------------------------------------------------------
  ## 6. Z grid consistency
  ## ---------------------------------------------------------
  expect_true(all(vapply(res, function(df) {
    identical(df$Z, z_values)
  }, logical(1))))

  ## ---------------------------------------------------------
  ## 7. Density sanity checks
  ## ---------------------------------------------------------
  expect_true(all(vapply(res, function(df) {
    all(is.finite(df$Density))
  }, logical(1))))

  expect_true(all(vapply(res, function(df) {
    all(df$Density >= 0)
  }, logical(1))))

  ## ---------------------------------------------------------
  ## 8. No completely empty outputs
  ## ---------------------------------------------------------
  expect_true(all(vapply(res, nrow, integer(1)) > 0))
})
