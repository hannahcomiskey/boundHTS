#' Get logit transformation of Sigma variance-covariance matrix via Delta method
#'
#' @param mu estimated mean
#' @param Sigma estimated variance-covariance matrix
#' @details
#' The logit transformation of Sigma variance-covariance matrix via the Delta method
#'
#' @return logit transformed variance-covariance matrix
#' @export
#'
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
