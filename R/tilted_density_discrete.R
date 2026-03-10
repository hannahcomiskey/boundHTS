#' Tilt the Density of f(y) for Discrete Outcomes Using the Optimized Nu Parameter
#'
#'
#' @param nu Numeric; the tilting parameter.
#' @param f_y Numeric vector; the discrete density values of the aggregated
#'   series evaluated at `y_vals`.
#' @param y_vals Numeric vector; evaluation points corresponding to `f_y`.
#'
#' @details
#' Exponential tilting adjusts the density to shift the mean (or other moments)
#' according to a tilting parameter `nu`. For discrete densities, the tilted
#' density is computed as:
#'
#' \deqn{f_{\rm tilt}(y_i) = \frac{\exp(\nu y_i) f_y[i]}{\sum_j \exp(\nu y_j) f_y[j]}}
#'
#' @return Numeric vector of the same length as `y_vals`, representing the
#'   normalized tilted discrete density.
#'
#' @examples
#' y_vals <- 0:5
#' f_y <- dbinom(y_vals, size = 5, prob = 0.5)
#' nu <- 0.5
#' f_tilt <- tilted_density_discrete(nu, f_y, y_vals)
#' sum(f_tilt)  # should be 1
#'
#' @export

tilted_density_discrete <- function(nu, f_y, y_vals) {
  log_weights <- nu * y_vals + log(f_y) # log MGF
  w <- exp(log_weights)
  f_tilt <- w / sum(w) # normalise
  return(f_tilt)
}
