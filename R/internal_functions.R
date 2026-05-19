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


#' Recursive weight propagation
#'
#' @description Returns weights matrix
#' @param node name of nodes
#' @param W empty weights matrix
#' @param weights_list List of aggregation weights
#' @param bottom_names Names of bottom nodes
#' @export
propagate_weights <- function(weights_list, node, W, bottom_names) {

  ## Already computed
  if (sum(W[node, ]) > 0) {
    return(W[node, ])
  }

  ## Local children weights
  w_local <- weights_list[[node]]

  if (is.null(w_local)) {
    stop(paste0(
      "No local weights supplied for parent node '",
      node,
      "'."
    ))
  }

  ## Weight validation
  if (!is.numeric(w_local)) {
    stop(paste0(
      "Weights for node '",
      node,
      "' must be numeric."
    ))
  }

  if (abs(sum(w_local) - 1) > 1e-8) {
    stop(paste0(
      "Weights for node '",
      node,
      "' must sum to 1."
    ))
  }

  child_nodes <- names(w_local)

  if (is.null(child_nodes)) {
    stop(paste0(
      "Weights for node '",
      node,
      "' must be named."
    ))
  }

  ## Aggregate recursively
  node_weights <- numeric(length(bottom_names))
  names(node_weights) <- bottom_names

  for (k in seq_along(child_nodes)) {

    child <- child_nodes[k]
    child_weight <- w_local[k]

    ## Recursively compute child weights
    child_bottom_weights <- propagate_weights(weights_list = weights_list,
                                              node = child,
                                              W = W,
                                              bottom_names = bottom_names)

    node_weights <- node_weights +
      child_weight * child_bottom_weights
  }

  return(node_weights)
}
