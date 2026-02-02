#' Tilting the density of f(y) using the optimised nu parameter.
#'
#' @param f_y calculate density of aggregated series
#' @param y_vals evaluation points over which the density is calculated
#' @param nu The tilting parameter
#' @return The tilted density of f(y).
#' @export

tilted_density_cont <- function(nu, f_y, y_vals) {
  log_weights <- nu * y_vals + log(f_y) # log MGF
  w <- exp(log_weights)
  f_tilt <- w / pracma::trapz(y_vals, w) # normalise
  return(f_tilt)
}
