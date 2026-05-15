test_that("ZOIB_convolution_density returns a finite scalar", {

  set.seed(123)

  ## Define a simple two-level hierarchy
  groups <- list(
    bottom = c("A", "B"),   # bottom level
    top = c("Total")     # aggregated level
  )

  ## Point estimates for Beta parameters
  alpha_list <- list(
    c(2, 5),
    c(7)
  )

  beta_list <- list(
    c(6, 3),
    c(4)
  )

  zoi_list <- list(
    c(0.1, 0.1),
    c(0.05)
  )

  coi_list <- list(
    c(0.05, 0.05),
    c(0.05)
  )


  ## Aggregation weights
  weights_list <- list(
    c(0.4, 0.6),
    1
  )

  ## Estimate densities
  dens <- ZOIB_convolution(
    groups = groups,
    alpha_list = alpha_list,
    beta_list = beta_list,
    zoi_list =   zoi_list,
    coi_list = coi_list,
    weights_list = weights_list,
    point = TRUE,
    n_draws = 500,
    n_sims = 50
  )

  dens

  # ---- expectations ----
  expect_type(dens, "list")
  expect_length(dens[[1]]$Density, 1000)
  expect_true(all(is.finite(dens[[1]]$Density)))
  expect_true(all(dens[[1]]$Density >= 0))
})

test_that("ZOIB_convolution_density is stable when z is out of support", {

  set.seed(123)

  ## Define a simple two-level hierarchy
  groups <- list(
    bottom = c("A", "B"),   # bottom level
    top = c("Total")     # aggregated level
  )

  ## Point estimates for Beta parameters
  alpha_list <- list(
    c(2, 5),
    c(7)
  )

  beta_list <- list(
    c(6, 3),
    c(4)
  )

  zoi_list <- list(
    c(0.1, 0.1),
    c(0.05)
  )

  coi_list <- list(
    c(0.05, 0.05),
    c(0.05)
  )


  ## Aggregation weights
  weights_list <- list(
    c(0.4, 0.6),
    1
  )

  ## Estimate densities
  expect_error(ZOIB_convolution(
    groups = groups,
    alpha_list = alpha_list,
    beta_list = beta_list,
    zoi_list =   zoi_list,
    coi_list = coi_list,
    weights_list = weights_list,
    point = TRUE,
    n_draws = 500,
    n_sims = 50,
    z_values = seq(-10, 10, length.out=200)
  ))

})

test_that("ZOIB_convolution_density errors on incompatible dimensions", {

  set.seed(123)

  ## Define a simple two-level hierarchy
  groups <- list(
    bottom = c("A", "B"),   # bottom level
    top = c("Total")     # aggregated level
  )

  ## Point estimates for Beta parameters
  alpha_list <- list(
    c(2, 5),
    c(7)
  )

  beta_list <- list(
    c(6, 3, 4),
    c(2,3),
    c(4)
  )

  zoi_list <- list(
    c(0.1, 0.1),
    c(0.05)
  )

  coi_list <- list(
    c(0.05, 0.05),
    c(0.05)
  )

  ## Aggregation weights
  weights_list <- list(
    c(0.4, 0.6),
    1
  )

  ## Estimate densities
  expect_error(ZOIB_convolution(
    groups = groups,
    alpha_list = alpha_list,
    beta_list = beta_list,
    zoi_list =   zoi_list,
    coi_list =   coi_list,
    weights_list = weights_list,
    point = TRUE,
    n_draws = 500,
    n_sims = 50,
    z_values = seq(-10, 10, length.out=200)))

})
