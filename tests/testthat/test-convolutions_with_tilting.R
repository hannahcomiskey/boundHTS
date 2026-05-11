test_that("convolution_with_tilting returns correctly structured output", {

  set.seed(123)

  groups <- list(
    State = c("NSW", "VIC"),
    National = "AUS"
  )

  lambda_list <- list(
    c(10, 12),
    22
  )

  z_values <- 0:30
  mu_theory <- list(22, 22)

  res <- convolution_with_tilting(
    family = "Poisson",
    groups = groups,
    point = TRUE,
    mu_theory = mu_theory,
    z_values = z_values,
    lambda_list = lambda_list
  )

  ## ---- structure checks ----
  expect_type(res, "list")

  expect_true(all(c(
    "base_density",
    "tilted_density",
    "tilted_samples",
    "tilting_parameter"
  ) %in% names(res)))

  expect_length(res$base_density, 1)
  expect_length(res$tilted_density, 1)

  ## ---- no NULL outputs ----
  expect_false(any(vapply(res$base_density, is.null, logical(1))))
  expect_false(any(vapply(res$tilted_density, is.null, logical(1))))

  expect_false(any(vapply(res$tilted_samples, is.null, logical(1))))

  ## ---- tilting parameters exist ----
  expect_true(all(vapply(res$tilting_parameter, is.numeric, logical(1))))
})


test_that("convolution and tilting preserve grid consistency", {

  set.seed(123)

  groups <- list(
    State = c("A", "B"),
    National = "C"
  )

  lambda_list <- list(
    c(8, 10),
    18
  )

  z_values <- 0:25
  mu_theory <- list(15, 15)

  res <- convolution_with_tilting(
    family = "Poisson",
    groups = groups,
    point = TRUE,
    mu_theory = mu_theory,
    z_values = z_values,
    lambda_list = lambda_list
  )

  expect_true(all(vapply(res$base_density, function(df) {
    identical(df$Z, z_values)
  }, logical(1))))

  expect_true(all(vapply(res$tilted_density, function(df) {
    identical(df$Z, z_values)
  }, logical(1))))
})


test_that("density outputs are valid probability distributions", {

  set.seed(123)

  groups <- list(
    State = c("NSW", "VIC"),
    National = "AUS"
  )

  lambda_list <- list(
    c(6, 9),
    15
  )

  z_values <- 0:20
  mu_theory <- list(10, 10)

  res <- convolution_with_tilting(
    family = "Poisson",
    groups = groups,
    point = TRUE,
    mu_theory = mu_theory,
    z_values = z_values,
    lambda_list = lambda_list
  )

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

  groups <- list(
    State = c("A", "B"),
    National = "C"
  )

  lambda_list <- list(
    c(3, 4),
    7
  )

  res <- convolution_with_tilting(
    family = "Poisson",
    groups = groups,
    point = TRUE,
    mu_theory = list(15, 15),
    z_values = 0:20,
    lambda_list = lambda_list
  )

  diffs <- mapply(function(base, tilt) {
    sum(abs(base$Density - tilt$Density))
  }, res$base_density, res$tilted_density)

  expect_true(any(diffs > 1e-6))
})


test_that("tilted samples are valid numeric vectors", {

  set.seed(123)

  groups <- list(
    State = c("A", "B"),
    National = "C"
  )

  lambda_list <- list(
    c(5, 6),
    10
  )

  res <- convolution_with_tilting(
    family = "Poisson",
    groups = groups,
    point = TRUE,
    mu_theory = list(10, 10),
    z_values = 0:15,
    lambda_list = lambda_list
  )

  expect_true(all(vapply(res$tilted_samples, is.numeric, logical(1))))
  expect_true(all(vapply(res$tilted_samples, length, integer(1)) > 0))
})


test_that("function fails gracefully on invalid family", {

  expect_error(
    convolution_with_tilting(
      family = "Gaussian",
      groups = list(A = 1),
      point = TRUE,
      mu_theory = list(1),
      z_values = 0:5
    ),
    "not"
  )
})
