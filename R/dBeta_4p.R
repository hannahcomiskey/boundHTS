#' Calculate 4 parameter Beta PDF
#'
#' @param x evaluation points
#' @param alpha estimated shape parameter of Beta distribution
#' @param beta estimated shape parameter of Beta distribution
#' @param lower lower bound
#' @param upper upper bound
#' @details
#' A wrapper function for `ExtDist::dBeta_ab`.
#'
#' @return The density of the 4 parameter Beta distribution
#' @export
#' @examples
#' \dontrun{
#' dBeta_4p(x = 0.07, lower=0, upper=0.5, alpha=5, beta=2)
#' }

dBeta_4p <- function(x, lower=0, upper, alpha, beta) {
  x_scaled <- x / upper  # Scale to [0,1]
  ifelse(x_scaled < 0 | x_scaled > 1, 0,
         ifelse( alpha < 0 | beta < 0, 0,
                 ExtDist::dBeta_ab(x, alpha, beta, lower, upper)
         )
  )
}
