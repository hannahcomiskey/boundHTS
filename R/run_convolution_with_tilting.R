#' Convolution with Exponential Tilting Wrapper
#'
#' @description
#' This function combines hierarchical convolution of predictive densities
#' with exponential tilting to align the resulting distributions with
#' externally specified target means.
#'
#' It is designed for probabilistic hierarchical forecasting workflows
#' where base predictive distributions (e.g. Poisson, Beta, ZOIB)
#' are first constructed via convolution and then post-processed using
#' exponential tilting for moment correction.
#'
#' @param family Character. Distribution family used in convolution.
#' Must be one of: `"Poisson"`, `"Beta"`, `"ZIB"`, `"ZOIB"`.
#' @param S Named list defining hierarchical structure from bottom
#' to top level. Each element contains node names at that level.
#' @param point Logical. If TRUE, uses point estimates of parameters.
#' If FALSE, uses posterior samples.
#' @param mu_theory Numeric list. Target means for each hierarchy level
#' used in exponential tilting.
#' @param z_values Numeric vector. Grid of evaluation points.
#' If NULL and Beta-family is used internally, defaults to a unit grid.
#' @param ... Additional arguments passed to `run_convolution()`.
#' These include:
#' \itemize{
#'   \item lambda_list (Poisson)
#'   \item alpha_list, beta_list, weights_list (Beta)
#'   \item alpha_list, beta_list, zoi_list, coi_list, weights_list (ZOIB)
#' }
#'
#' @details
#' The function proceeds in two stages:
#'
#' 1. Hierarchical convolution is performed using `run_convolution()`
#'    to generate base predictive densities.
#'
#' 2. Each level is independently adjusted using exponential tilting
#'    via `tilt_density()` to match a target mean.
#'
#' The output preserves hierarchical structure and returns both
#' raw and tilted densities.
#'
#' @return A list with components:
#' \describe{
#'   \item{base_density}{List of original convolution densities per level.}
#'   \item{tilted_density}{List of tilted densities per level.}
#'   \item{tilted_samples}{Samples drawn from tilted densities.}
#'   \item{tilting_parameter}{Estimated tilting parameters (nu) per level.}
#' }
#'
#' @examples
#'
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
#' lambda_mat <- c(Total = 8, A = 5, B = 6, A1 = 2, A2 = 4, B1 = 3, B2 = 5)
#'
#' z_values <- seq(0, 100)
#'
#' dens <- run_convolution_with_tilting(family = "Poisson",
#'  S = S,
#'  mu_theory = c(Total = 9, A = 4, B = 10, A1 = 1, A2 = 4, B1 = 2, B2 = 7),
#'  lambda_mat = lambda_mat,
#'  z_values = z_values,
#'  point = TRUE,
#'  n_sims = 50,
#'  n_draws = 100)
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
#' weights_list <- list( Total = c(A1 = 0.10, A2 = 0.15, B1 = 0.30, B2 = 0.45),
#' A = c(A1 = 0.40, A2 = 0.60), B = c(B1 = 0.30, B2 = 0.70))
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
#' mu_theory <- c( Total = 0.55, A = 0.18, B = 0.37, A1 = 0.05, A2 = 0.11, B1 = 0.14, B2 = 0.31)
#'
#' dens <- run_convolution_with_tilting(family = "Beta",
#'   S = S_beta,
#'   point = TRUE,
#'   mu_theory = mu_theory,
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
#' weights_list <- list(
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
#' mu_theory <- c( Total = 0.60, A = 0.22, B = 0.40, A1 = 0.07, A2 = 0.14, B1 = 0.16, B2 = 0.34)
#'
#' dens_zoib <- run_convolution_with_tilting(family = "ZOIB",
#'   S = S_beta,
#'   point = TRUE,
#'   mu_theory = mu_theory,
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
#'
#' @export
#'

run_convolution_with_tilting <- function(family, S, point, mu_theory, z_values, ...) {

  ## ---------------------------------------------------------------
  ## 1. Run convolution step
  ## ---------------------------------------------------------------

  base_list <- run_convolution(family = family,
                               S = S,
                               point = point,
                               z_values = z_values,
                               ...
                               )

  n_levels <- length(base_list)

  ## Storage objects
  tilted_density <- vector("list", n_levels)
  tilting_param <- vector("list", n_levels)
  tilted_samples <- vector("list", n_levels)

  ## ---------------------------------------------------------------
  ## 2. Apply exponential tilting per level
  ## ---------------------------------------------------------------
  for (i in seq_len(n_levels)) {

    base_df <- base_list[[i]]

    node_name <- unique(base_df$Node)
    z_vals <- base_df$Z
    f_y <- base_df$Density

    ## Discrete for Poisson, continuous otherwise
    discrete_flag <- identical(family, "Poisson") # T/F

    tilt_res <- tilt_density(
      mu_theory = mu_theory[node_name],
      z_vals = z_vals,
      f_y = f_y,
      discrete = discrete_flag
    )

    tilted_density[[i]] <- data.frame(
      Node = node_name,
      Z = z_vals,
      Density = tilt_res$f_tilted
    )

    tilting_param[[i]] <- tilt_res$nu_star
    tilted_samples[[i]] <- tilt_res$tilted_samps
  }

  ## ---------------------------------------------------------------
  ## 3. Return structured output
  ## ---------------------------------------------------------------
  list(base_density = base_list,
       tilted_density = tilted_density,
       tilted_samples = tilted_samples,
       tilting_parameter = tilting_param
       )
}
