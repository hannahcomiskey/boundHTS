test_that("dZOIB_4p returns a numeric vector of correct length", {

  set.seed(123)

  n_nodes <- 3

  z <- seq(0, 1, length.out=100)

  alpha_point <- runif(n_nodes, min = 2, max = 10)

  beta_point <- runif(n_nodes, min = 2, max = 10)

  zoi_point <- runif(n_nodes, 0.1, 0.3)

  coi_point <- runif(n_nodes, 0.2, 0.8)

  upper <- rep(1, n_nodes)

  dens <- dZOIB_4p(
    x = z,
    alpha_point =  alpha_point[3],
    beta_point = beta_point[3],
    zi_point =  zoi_mc[3],
    coi_point = coi_mc[3],
    weight = upper[3]
  )

  expect_type(dens, "double")
  expect_length(dens, 100)
  expect_true(all(is.finite(dens)))
  expect_true(all(dens >= 0))
})

test_that("dZOIB_4p handles single-node case correctly", {

  set.seed(1)

  z <- 0.4

  alpha_point <- runif(1, min = 2, max = 10)
  beta_point <- runif(1, min = 2, max = 10)
  zoi_point <- runif(1, 0.1, 0.3)
  coi_point <- runif(1, 0.2, 0.8)
  upper <- 1

  dens <- dZOIB_4p(
    x = z,
    alpha_point =  alpha_point,
    beta_point = beta_point,
    zi_point =  zoi_point,
    coi_point = coi_point,
    weight = upper
  )

  expect_length(dens, 1)
  expect_true(all(is.finite(dens)))
})

test_that("dZOIB_4p assigns correct mass at zero and one boundaries", {

  n_nodes <- 2
  z <- seq(0, 1, length.out=10)
  alpha_point <- runif(1, min = 2, max = 10)
  beta_point <- runif(1, min = 2, max = 10)
  zoi_point <- runif(1, 0.1, 0.3)
  coi_point <- runif(1, 0.2, 0.8)
  upper <- 1

  dens0 <- dZOIB_4p(
    x = z[1],
    alpha_point =  alpha_point,
    beta_point = beta_point,
    zi_point =  zoi_point,
    coi_point = coi_point,
    weight = upper
  )

  expect_equal(
    dens0,
    zoi_point * (1 - coi_point)
  )

  dens1 <- dZOIB_4p(
    x = z[10],
    alpha_point =  alpha_point,
    beta_point = beta_point,
    zi_point =  zoi_point,
    coi_point = coi_point,
    weight = upper
  )

  expect_equal(
    dens1,
    zoi_point * coi_point
  )
})

test_that("dZOIB_4p reduces to Beta density when zoi = 0", {

  set.seed(42)

  n_nodes <- 2
  z <- seq(0, 1, length.out=10)
  alpha_point <- runif(1, min = 2, max = 10)
  beta_point <- runif(1, min = 2, max = 10)
  zoi_point <- 0
  coi_point <- runif(1, 0.2, 0.8)
  upper <- 1

  dens <- dZOIB_4p(
    x = z,
    alpha_point =  alpha_point,
    beta_point = beta_point,
    zi_point =  zoi_point,
    coi_point = coi_point,
    weight = upper
  )
  expect_equal(
    dens,
    ExtDist::dBeta_ab(z, alpha_point, beta_point, 0, 1),
    tolerance = 1e-12
  )
})

