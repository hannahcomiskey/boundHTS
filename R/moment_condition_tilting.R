#' The moment condition via the cumulant generating function (CGF) for continuous densities.
#'
#' @param f_y calculate density of aggregated
#' @param y_vals evaluation points over which the density is calculated
#' @param nu The tilting parameter
#' @param mu_theory The theorised (target) mean
#' @details
#' A loss function to optimise over to evaluate the tilting parameter (nu) used in exponential tilting.
#'
#' @return The residuals in means (estimated - target) after tilting with a given nu
#' @export

moment_condition_tilting <- function(nu, f_y, y_vals, mu_theory) {

  # Grid
  dz <- mean(diff(y_vals))

  # Stabilise
  log_w <- nu * y_vals + log(f_y)
  w <- exp(log_w)

  # MGF terms
  mgf_stable <- sum(w * dz)
  deriv_cgf <- sum(y_vals * w * dz) / mgf_stable

  # Residual
  residual <- deriv_cgf - mu_theory
  return(residual)
}
