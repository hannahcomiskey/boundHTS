#' Random Generator for a Zero-One-Inflated Four-Parameter Beta Distribution
#'
#' @param n_mc Integer; number of Monte Carlo samples to generate.
#' @param alpha_point Numeric; first shape parameter of the Beta distribution.
#' @param beta_point Numeric; second shape parameter of the Beta distribution.
#' @param zoi_point Numeric; probability of zero/one inflation.
#' @param coi_point Numeric; conditional probability of one given inflation.
#' @param lower Numeric; lower bound of the Beta distribution (default = 0).
#' @param weight Numeric; upper bound of the Beta distribution (default = 1).
#'
#' @details
#' This function generates random samples from a zero-one-inflated Beta
#' distribution defined on the interval \eqn{[0, weight]}. For each draw:
#'
#' \itemize{
#' \item With probability \code{zoi}, the observation is drawn from a
#'   zero-one inflation process.
#' \item Conditional on inflation, the value is \eqn{weight} with probability
#'   \code{coi} and \eqn{0} otherwise.
#' \item With probability \eqn{1 - zoi}, the observation is drawn from a
#'   Beta distribution with parameters \code{alpha} and \code{beta},
#'   scaled to the interval \eqn{[0, weight]} using
#'   \code{ExtDist::rBeta_ab()}.
#' }
#'
#' The implementation is vectorised for efficient Monte Carlo simulation.
#'
#' @return A numeric vector of length \code{n_mc} containing random draws from
#' the zero-one-inflated Beta distribution.
#'
#' @examples
#' set.seed(1)
#'
#' draws <- rZOIB_4p(
#'   n_mc = 1000,
#'   alpha_point = 2,
#'   beta_point = 10,
#'   zoi_point = 0.1,
#'   coi_point = 0.2,
#'   weight = 1
#' )
#'
#' hist(draws, breaks = 40)
#' @importFrom ExtDist rBeta_ab
#' @export

rZOIB_4p <- function(n_mc, alpha_point, beta_point, zoi_point, coi_point, lower=0, weight = 1) {

  Y <- numeric(n_mc)

  for(i in seq_len(n_mc)) {

    inflate <- stats::runif(1) < zoi_point

    if (inflate) {
      # zero or one inflation constrained by weight
      Y[i] <- stats::rbinom(1, 1, prob = coi_point) * weight
    } else {
      # weighted beta component
      Y[i] <- ExtDist::rBeta_ab(1, alpha_point, beta_point, lower, weight)
    }
  }

  return(Y)
}
