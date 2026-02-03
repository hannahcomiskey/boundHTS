test_that("dZOIB_4p returns a numeric vector of correct length", {

  set.seed(123)

  n_mc    <- 50
  n_nodes <- 3

  z <- 0.5

  Y_mc <- matrix(
    runif(n_mc * n_nodes, 0.05, 0.3),
    nrow = n_mc,
    ncol = n_nodes
  )

  phi_mc <- matrix(
    runif(n_mc * n_nodes, 5, 20),
    nrow = n_mc,
    ncol = n_nodes
  )

  zoi_mc <- matrix(
    runif(n_mc * n_nodes, 0.1, 0.3),
    nrow = n_mc,
    ncol = n_nodes
  )

  coi_mc <- matrix(
    runif(n_mc * n_nodes, 0.2, 0.8),
    nrow = n_mc,
    ncol = n_nodes
  )

  upper <- rep(1, n_nodes)

  dens <- dZOIB_4p(
    z = z,
    Y_mc = Y_mc,
    phi_mc = phi_mc,
    zoi_mc = zoi_mc,
    coi_mc = coi_mc,
    upper = upper
  )

  expect_type(dens, "double")
  expect_length(dens, n_mc)
  expect_true(all(is.finite(dens)))
  expect_true(all(dens >= 0))
})

test_that("dZOIB_4p handles single-node case correctly", {

  set.seed(1)

  n_mc <- 20
  z    <- 0.4

  Y_mc <- matrix(runif(n_mc, 0.1, 0.5), ncol = 1)
  phi_mc <- matrix(runif(n_mc, 5, 10), ncol = 1)
  zoi_mc <- matrix(runif(n_mc, 0.2, 0.4), ncol = 1)
  coi_mc <- matrix(runif(n_mc, 0.3, 0.6), ncol = 1)
  upper  <- 1

  dens <- dZOIB_4p(
    z = z,
    Y_mc = Y_mc,
    phi_mc = phi_mc,
    zoi_mc = zoi_mc,
    coi_mc = coi_mc,
    upper = upper
  )

  expect_length(dens, n_mc)
  expect_true(all(is.finite(dens)))
})

test_that("dZOIB_4p assigns correct mass at zero and one boundaries", {

  n_mc <- 4
  n_nodes <- 2

  # Construct deterministic example
  Y_mc <- matrix(
    c(0.1, 0.2,
      0.1, 0.2,
      0.1, 0.2,
      0.1, 0.2),
    ncol = n_nodes,
    byrow = TRUE
  )

  phi_mc <- matrix(10, n_mc, n_nodes)

  zoi_mc <- matrix(
    c(0.6, 0.7,
      0.6, 0.7,
      0.6, 0.7,
      0.6, 0.7),
    n_mc, n_nodes, byrow = TRUE
  )

  coi_mc <- matrix(
    c(0.2, 0.3,
      0.2, 0.3,
      0.2, 0.3,
      0.2, 0.3),
    n_mc, n_nodes, byrow = TRUE
  )

  upper <- c(1, 1)

  z0 <- rowSums(Y_mc[, -n_nodes, drop=FALSE])

  dens0 <- dZOIB_4p(
    z = z0,
    Y_mc = Y_mc,
    phi_mc = phi_mc,
    zoi_mc = zoi_mc,
    coi_mc = coi_mc,
    upper = upper
  )

  expect_equal(
    dens0,
    zoi_mc[, n_nodes] * (1 - coi_mc[, n_nodes])
  )

  # z so that x_vals >= upper
  z1 <- rowSums(Y_mc[, -n_nodes, drop=FALSE]) + upper[n_nodes]

  dens1 <- dZOIB_4p(
    z = z1,
    Y_mc = Y_mc,
    phi_mc = phi_mc,
    zoi_mc = zoi_mc,
    coi_mc = coi_mc,
    upper = upper
  )

  expect_equal(
    dens1,
    zoi_mc[, n_nodes] * coi_mc[, n_nodes]
  )
})

test_that("dZOIB_4p reduces to Beta density when zoi = 0", {

  set.seed(42)

  n_mc    <- 10
  n_nodes <- 2

  z <- 0.4

  Y_mc <- matrix(runif(n_mc * n_nodes, 0.05, 0.2), n_mc, n_nodes)
  phi_mc <- matrix(runif(n_mc * n_nodes, 10, 20), n_mc, n_nodes)

  zoi_mc <- matrix(0, n_mc, n_nodes)  # no zero-one inflation
  coi_mc <- matrix(runif(n_mc * n_nodes, 0.3, 0.7), n_mc, n_nodes)

  upper <- rep(1, n_nodes)

  dens <- dZOIB_4p(
    z = z,
    Y_mc = Y_mc,
    phi_mc = phi_mc,
    zoi_mc = zoi_mc,
    coi_mc = coi_mc,
    upper = upper
  )

  child_sum <- rowSums(Y_mc[, -n_nodes, drop=FALSE])
  x_vals <- z - child_sum

  alpha <- pmax(Y_mc[, n_nodes] * phi_mc[, n_nodes], 1e-4)
  beta  <- pmax((1 - Y_mc[, n_nodes]) * phi_mc[, n_nodes], 1e-4)

  expect_equal(
    dens,
    ExtDist::dBeta_ab(x_vals, alpha, beta, 0, 1),
    tolerance = 1e-12
  )
})

