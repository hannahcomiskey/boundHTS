#' Logit Transformation of Variance-Covariance Matrix via Delta Method
#'
#' Computes the variance-covariance matrix of the logit-transformed means
#' using the Delta method. This is useful when you have estimated means
#' constrained between 0 and 1 and wish to approximate the covariance of
#' their logit-transformed values.
#'
#' @param mu Numeric vector of estimated means (values between 0 and 1).
#' @param Sigma Numeric variance-covariance matrix corresponding to \code{mu}.
#'
#' @details
#' The function applies the Delta method to obtain an approximate variance-covariance
#' matrix of the logit-transformed means. The Jacobian of the logit transformation
#' is computed elementwise and then applied to the original variance-covariance matrix:
#'
#' \deqn{Var(logit(\mu)) \approx J \, \Sigma \, J, \quad
#' J = diag\left( \frac{1}{\mu_i (1 - \mu_i)} \right)}
#'
#' @return A numeric matrix representing the approximate variance-covariance
#' matrix of the logit-transformed means.
#'
#' @examples
#' mu <- c(0.2, 0.5, 0.8)
#' Sigma <- matrix(c(0.01, 0.002, 0.001,
#'                   0.002, 0.02, 0.003,
#'                   0.001, 0.003, 0.015), nrow = 3, byrow = TRUE)
#'
#' logit_vcov_delta(mu, Sigma)
#'
#' @export

logit_vcov_delta <- function(mu, Sigma) {
  stopifnot(
    is.numeric(mu),
    is.matrix(Sigma),
    length(mu) == nrow(Sigma),
    nrow(Sigma) == ncol(Sigma)
  )

  # Jacobian of elementwise logit
  J <- diag(1 / (mu * (1 - mu))) # derivative of logit(mu/1-mu)

  # Delta-method variance–covariance
  J %*% Sigma %*% J
}
