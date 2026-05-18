#' Calculate alpha and beta params for beta distribution from mean and variance
#'
#' @param mean estimated mean
#' @param sd estimated standard deviation
#' @details
#' The corresponding alpha and beta shape parameters
#'
#' @return A vector of shape parameters
#' @export

beta_params <- function(mean, sd) {
  var <- sd^2
  alpha <- ((1 - mean) / var - 1 / mean) * mean^2
  beta  <- alpha * (1 / mean - 1)
  return(c(alpha, beta))
}

#' Define safe exponential (to avoid numbers blowing up)
#' @param x vector of values to exponentiate
#' @export

safe_exp <- function(x) {
  exp(pmin(x, 700))
}

#' To handle ... in functions
#'
#' @description
#' Returns `y` if `x` is `NULL`, otherwise returns `x`.
#'
#' @param x vector of values
#' @param y vector of values
#' @export
default_settings <- function(x, y) {
  if (is.null(x)) y else x
}
