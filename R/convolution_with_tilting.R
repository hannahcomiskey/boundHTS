#' Convolution + Exponential Tilting Wrapper
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
#' @param groups Named list defining hierarchical structure from bottom
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
#' set.seed(123)
#'
#' groups <- list(
#'   State = c("NSW", "VIC"),
#'   National = c("AUS")
#' )
#'
#' lambda_list <- list(
#'   c(12, 10),
#'   22
#' )
#'
#' z_values <- 0:60
#'
#' mu_theory <- list(22, 22)
#'
#' res <- convolution_with_tilting(
#'   family = "Poisson",
#'   groups = groups,
#'   point = TRUE,
#'   mu_theory = mu_theory,
#'   z_values = z_values,
#'   lambda_list = lambda_list
#' )
#'
#' head(res$base_density[[1]])
#'
#' @export
#'

convolution_with_tilting <- function(family,
                                     groups,
                                     point,
                                     mu_theory,
                                     z_values,
                                     ...) {

  ## ---------------------------------------------------------------
  ## 1. Run convolution step
  ## ---------------------------------------------------------------
  base_list <- run_convolution(
    family = family,
    groups = groups,
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

    z_vals <- base_df$Z
    f_y    <- base_df$Density

    ## Discrete for Poisson, continuous otherwise
    discrete_flag <- identical(family, "Poisson")

    tilt_res <- tilt_density(
      mu_theory = mu_theory[[i]],
      z_vals = z_vals,
      f_y = f_y,
      discrete = discrete_flag
    )

    tilted_density[[i]] <- data.frame(
      Z = z_vals,
      Density = tilt_res$f_tilted
    )

    tilting_param[[i]] <- tilt_res$nu_star
    tilted_samples[[i]] <- tilt_res$tilted_samps
  }

  ## ---------------------------------------------------------------
  ## 3. Return structured output
  ## ---------------------------------------------------------------
  list(
    base_density = base_list,
    tilted_density = tilted_density,
    tilted_samples = tilted_samples,
    tilting_parameter = tilting_param
  )
}
