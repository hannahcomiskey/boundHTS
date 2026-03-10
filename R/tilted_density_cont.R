#' Tilt the Density of f(y) Using the Optimized Nu Parameter
#'
#' @param nu Numeric; the tilting parameter.
#' @param f_y Numeric vector; the density values of the aggregated series
#'   evaluated at `y_vals`.
#' @param y_vals Numeric vector; evaluation points corresponding to `f_y`.
#'
#' @details
#' Exponential tilting adjusts the density to shift the mean (or other moments)
#' according to a tilting parameter `nu`. The tilted density is computed as:
#'
#' \deqn{f_{\rm tilt}(y) = \frac{\exp(\nu y) f_y}{\int \exp(\nu y) f_y dy}}
#'
#' The integral in the denominator is approximated numerically using the
#' trapezoidal rule via \code{pracma::trapz}.
#'
#' @return Numeric vector of the same length as `y_vals`, representing the
#'   normalized tilted density.
#'
#' @examples
#' y_vals <- seq(-3, 3, length.out = 100)
#' f_y <- dnorm(y_vals)
#' nu <- 0.5
#' f_tilt <- tilted_density_cont(nu, f_y, y_vals)
#' sum(f_tilt) * mean(diff(y_vals))  # approx. 1
#'
#' @export

tilted_density_cont <- function(nu, f_y, y_vals) {
  log_weights <- nu * y_vals + log(f_y) # log MGF
  w <- exp(log_weights)
  f_tilt <- w / pracma::trapz(y_vals, w) # normalise
  return(f_tilt)
}
