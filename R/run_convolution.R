#' Wrapper for generating convolution
#'
#' @param family The name of the probability distribution used in the convolution. Must be a string of one of these options: c("Poisson", "NB", "Beta", "ZIB", "ZOIB", "Gamma").
#' @param groups List of nodes at each level in series. The length of the list is the number of tiers in the hierarchy.
#' @param point A true/false indicator to denote whether you are using point estimates (point=TRUE) or posterior samples (point=FALSE) of the Poisson parameters.
#' @param z_values A grid of evaluation points. If the family is in the Beta-family, this can be left NULL and it defaults to a vector of length 1000 over the unit interval.
#' @param ... Parameters specific to the parameter family. See "Poisson_convolution", "Beta_convolution", "ZIB_convolution", "ZOIB_convolution" and "Gamma_convolution" for family-specific details.
#' @details
#' A wrapper function for calculating Monte Carlo estimates for the aggregate
#' density Z using point estimates or posterior samples of lambda.
#' @examples
#' ## ---------------------------------------------------------
#' ## Example: Point-estimate Poisson convolution
#' ## ---------------------------------------------------------
#' set.seed(123)
#'
#' S <- rbind(Total = c(1, 1, 1, 1),
#' A = c(1, 1, 0, 0),
#' B = c(0, 0, 1, 1),
#' A1 = c(1, 0, 0, 0),
#' A2 = c(0, 1, 0, 0),
#' B1 = c(0, 0, 1, 0),
#' B2 = c(0, 0, 0, 1))
#' colnames(S) <- c("A1", "A2", "B1", "B2")
#'
#' lambda_mat <- c(Total = 8,
#' A = 5,
#' B = 6,
#' A1 = 2,
#' A2 = 4,
#' B1 = 3,
#' B2 = 5)
#'
#' z_values <- seq(0, 100)
#'
#' dens <- run_convolution(family = "Poisson", S = S, lambda_mat = lambda_mat,
#' z_values = z_values, point = TRUE, n_sims = 50, n_draws = 100)
#'
#' head(dens[[1]])
#'
#' ## ---------------------------------------------------------
#' ## Example: Point-estimate Beta convolution
#' ## ---------------------------------------------------------
#'
#' set.seed(123)
#'
#' S_beta <- rbind(
#'   Total = c(1, 1, 1, 1),
#'   A     = c(1, 1, 0, 0),
#'   B     = c(0, 0, 1, 1),
#'   A1    = c(1, 0, 0, 0),
#'   A2    = c(0, 1, 0, 0),
#'   B1    = c(0, 0, 1, 0),
#'   B2    = c(0, 0, 0, 1)
#' )
#' colnames(S_beta) <- c("A1", "A2", "B1", "B2")
#'
#' local_weights <- list(
#'   Total = c(A1 = 0.10, A2 = 0.15, B1 = 0.30, B2 = 0.45),
#'   A     = c(A1 = 0.40, A2 = 0.60),
#'   B     = c(B1 = 0.30, B2 = 0.70)
#' )
#'
#' alpha_mat <- c(
#'   Total = 8,
#'   A = 5,
#'   B = 6,
#'   A1 = 2,
#'   A2 = 4,
#'   B1 = 3,
#'   B2 = 5
#' )
#'
#' beta_mat <- c(
#'   Total = 4,
#'   A = 3,
#'   B = 2,
#'   A1 = 6,
#'   A2 = 5,
#'   B1 = 7,
#'   B2 = 4
#' )
#'
#' dens <- run_convolution(
#'   family = "Beta",
#'   S = S_beta,
#'   point = TRUE,
#'   z_values = NULL,
#'   alpha_mat = alpha_mat,
#'   beta_mat = beta_mat,
#'   weights_list = weights_list,
#'   n_draws = 500,
#'   n_sims = 50
#' )
#'
#' ## Inspect results
#' head(dens[[1]])
#'
#'
#' ## ---------------------------------------------------------------
#' ## Example: Zero-one-inflated Beta convolution
#' ## ---------------------------------------------------------------
#'
#' set.seed(123)
#'
#' S_beta <- rbind(
#'   Total = c(1, 1, 1, 1),
#'   A = c(1, 1, 0, 0),
#'   B = c(0, 0, 1, 1),
#'   A1 = c(1, 0, 0, 0),
#'   A2 = c(0, 1, 0, 0),
#'   B1 = c(0, 0, 1, 0),
#'   B2 = c(0, 0, 0, 1)
#' )
#' colnames(S_beta) <- c("A1", "A2", "B1", "B2")
#'
#' local_weights <- list(
#'   Total = c(A1 = 0.10, A2 = 0.15, B1 = 0.30, B2 = 0.45),
#'   A = c(A1 = 0.40, A2 = 0.60),
#'   B = c(B1 = 0.30, B2 = 0.70)
#' )
#'
#' alpha_mat <- c(
#'   Total = 8,
#'   A = 5,
#'   B = 6,
#'   A1 = 2,
#'   A2 = 4,
#'   B1 = 3,
#'   B2 = 5
#' )
#'
#' beta_mat <- c(
#'   Total = 4,
#'   A = 3,
#'   B = 2,
#'   A1 = 6,
#'   A2 = 5,
#'   B1 = 7,
#'   B2 = 4
#' )
#'
#' zoi_mat <- c(
#'   Total = 0.05,
#'   A = 0.01,
#'   B = 0.03,
#'   A1 = 0.04,
#'   A2 = 0.03,
#'   B1 = 0.02,
#'   B2 = 0.01
#' )
#'
#' coi_mat <- c(
#'   Total = 0.01,
#'   A = 0.02,
#'   B = 0.02,
#'   A1 = 0.01,
#'   A2 = 0.01,
#'   B1 = 0.01,
#'   B2 = 0.01
#' )
#'
#' dens_zoib <- run_convolution(
#'   family = "ZOIB",
#'   S = S_beta,
#'   point = TRUE,
#'   z_values = NULL,
#'   alpha_mat = alpha_mat,
#'   beta_mat = beta_mat,
#'   zoi_mat = zoi_mat,
#'   coi_mat = coi_mat,
#'   weights_list = weights_list,
#'   n_draws = 500,
#'   n_sims = 50
#' )
#'
#' ## Inspect results
#' head(dens_zoib[[1]])
#'
#' @return The aggregate density Z over a grid of values.
#' @export
#'

run_convolution <- function(family, S, point, z_values, ...) {

  args <- list(...)

  if (!family %in% c("Poisson", "Beta", "ZIB", "ZOIB", "Gamma")) {
    stop("Family name not recognised. Please refer to help file for naming.", call. = FALSE)
  }

  if (family == "Poisson") {

    dens_list <- Poisson_convolution(S = S,
                                     z_values = z_values,
                                     lambda_list = args$lambda_list,
                                     point = point
                                     )
  } else if (family == "Beta") {

    dens_list <- Beta_convolution(S = S,
                                  alpha_mat = args$alpha_mat,
                                  beta_mat = args$beta_mat,
                                  weights_list = args$weights_list,
                                  point = point,
                                  n_draws = default_settings(args$n_draws, 2000),
                                  n_sims = default_settings(args$n_sims, 100),
                                  z_values = z_values
                                  )
  } else if (family == "ZIB") {

    dens_list <- ZIB_convolution(S = S,
                                 alpha_mat = args$alpha_mat,
                                 beta_mat = args$beta_mat,
                                 zi_mat = args$zi_mat,
                                 weights_list = args$weights_list,
                                 point = point,
                                 n_draws = default_settings(args$n_draws, 2000),
                                 n_sims = default_settings(args$n_sims, 100),
                                 z_values = z_values
                                 )
  } else if (family == "ZOIB") {

    dens_list <- ZOIB_convolution(S = S,
                                  alpha_mat = args$alpha_mat,
                                  beta_mat = args$beta_mat,
                                  zoi_mat = args$zoi_mat,
                                  coi_mat = args$coi_mat,
                                  weights_list = args$weights_list,
                                  point = point,
                                  n_draws = default_settings(args$n_draws, 2000),
                                  n_sims = default_settings(args$n_sims, 100),
                                  z_values = z_values
                                  )
  } else if (family == "Gamma") {

    dens_list <- Gamma_convolution(S = S,
                                   shape_mat = args$shape_mat,
                                   rate_mat = args$rate_mat,
                                   point = point,
                                   n_draws = default_settings(args$n_draws, 2000),
                                   n_sims = default_settings(args$n_sims, 100),
                                   z_values = z_values
                                   )
  }

  return(dens_list)
}
