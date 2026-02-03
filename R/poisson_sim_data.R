#' Simulated bounded hierarchical time series data
#'
#' A simulated dataset used to illustrate hierarchical coherent
#' forecasting for bounded time series.
#'
#' @format A data.frame with columns:
#' \describe{
#'   \item{X}{Covariate}
#'   \item{Tot}{Top-level count}
#'   \item{AA}{Bottom-level count (node AA)}
#'   \item{AB}{Bottom-level count (node AB)}
#'   \item{BA}{Bottom-level count (node BA)}
#'   \item{BB}{Bottom-level count (node BB)}
#' }
#' @source Simulated by the package authors
"poisson_sim_data"
