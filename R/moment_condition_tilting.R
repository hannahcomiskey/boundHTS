#' Moment Condition via the Cumulant Generating Function (CGF) for Continuous Densities
#'
#' Computes the moment condition used in exponential tilting by comparing the
#' expected value under a tilted density to a target theoretical mean. This is
#' useful when solving for the tilting parameter \code{nu} that shifts the
#' density to match a desired mean.
#'
#' @param nu Numeric tilting parameter.
#' @param f_y Numeric vector giving the density values of the aggregated variable
#'   evaluated at \code{y_vals}.
#' @param y_vals Numeric vector of evaluation points corresponding to \code{f_y}.
#' @param mu_theory Numeric scalar specifying the target mean.
#'
#' @details
#' The function computes a loss or residual based on the difference between the
#' expected value under the tilted density and the theoretical target mean:
#'
#' \deqn{R(\nu) = \frac{\sum y_i \exp(\nu y_i) f_y \Delta y}{\sum \exp(\nu y_i) f_y \Delta y} - \mu_{\rm theory}}
#'
#' where \eqn{\Delta y} is the spacing of \code{y_vals}. This residual can be
#' used in a root-finding procedure to solve for \code{nu}.
#'
#' @return Numeric scalar: the residual between the mean of the tilted density
#'   and the target mean.
#'
#' @examples
#' # Simple Gaussian example
#' y_vals <- seq(-3, 3, length.out = 100)
#' f_y <- dnorm(y_vals)
#' mu_theory <- 0.5
#' nu_start <- 0
#'
#' # Evaluate residual at nu = 0
#' moment_condition_tilting(nu_start, f_y, y_vals, mu_theory)
#'
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
