#' Random Generator for a Zero-Inflated Four-Parameter Beta Distribution
#'
#' @param n_mc Integer; number of Monte Carlo samples to generate.
#' @param alpha Numeric; first shape parameter of the Beta distribution.
#' @param beta Numeric; second shape parameter of the Beta distribution.
#' @param zoi Numeric; probability of zero/one inflation.
#' @param lower Numeric; lower bound of the Beta distribution (default = 0).
#' @param weight Numeric; upper bound of the Beta distribution (default = 1).
#'
#' @details
#' This function generates random samples from a zero-inflated Beta
#' distribution defined on the interval \eqn{[0, weight]}. For each draw:
#'
#' \itemize{
#' \item With probability \code{zoi}, the observation is drawn from a
#'   zero inflation process.
#' \item With probability \eqn{1 - zoi}, the observation is drawn from a
#'   Beta distribution with parameters \code{alpha} and \code{beta},
#'   scaled to the interval \eqn{[0, weight]} using
#'   \code{ExtDist::rBeta_ab()}.
#' }
#'
#' The implementation is vectorised for efficient Monte Carlo simulation.
#'
#' @return A numeric vector of length \code{n_mc} containing random draws from
#' the zero-inflated Beta distribution.
#'
#' @examples
#' set.seed(1)
#'
#' draws <- rZIB_4p(
#'   n_mc = 1000,
#'   alpha = 2,
#'   beta = 10,
#'   zoi = 0.1,
#'   weight = 1
#' )
#'
#' hist(draws, breaks = 40)
#'
#' @export

rZIB_4p <- function(n_mc, alpha_point, beta_point, zoi_point, lower=0, weight = 1) {

  Y <- numeric(n_mc)

  for(i in seq_len(n_mc)) {

    inflate <- stats::runif(1) < zoi_point

    if (inflate) {
      # zero inflation
      Y[i] <- 0
    } else {
      # weighted beta component
      Y[i] <- ExtDist::rBeta_ab(1, alpha_point, beta_point, lower, weight)
    }
  }

  return(Y)
}
