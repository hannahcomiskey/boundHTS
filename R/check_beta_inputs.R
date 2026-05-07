#' Validate inputs for Beta convolution functions
#'
#' Internal helper used to validate inputs supplied to the
#' Beta convolution density estimation functions.
#'
#' @param groups A list of node names at each hierarchical level.
#' @param alpha_list A list of alpha parameter vectors or matrices.
#' @param beta_list A list of beta parameter vectors or matrices.
#' @param weights_list A list of aggregation weight vectors.
#' @param point Logical indicating whether point estimates are supplied.
#' @param n_draws Positive integer giving the number of draws.
#' @param n_sims Positive integer giving the number of simulations.
#' @param z_values Optional numeric vector of evaluation points.
#'
#' @return Invisibly returns TRUE if all checks pass.
#'
#' @keywords internal
#' @noRd

check_beta_convolution_inputs <- function(groups, alpha_list, beta_list, weights_list,
                                          point, n_draws, n_sims, z_values) {

  ## -----------------------------
  ## Global checks
  ## -----------------------------

  if (!is.list(groups)) {
    stop("'groups' must be a list.")
  }

  if (!is.list(alpha_list)) {
    stop("'alpha_list' must be a list.")
  }

  if (!is.list(beta_list)) {
    stop("'beta_list' must be a list.")
  }

  if (!is.list(weights_list)) {
    stop("'weights_list' must be a list.")
  }

  if (length(groups) != length(alpha_list)) {
    stop("'groups' and 'alpha_list' must have the same length.")
  }

  if (length(groups) != length(beta_list)) {
    stop("'groups' and 'beta_list' must have the same length.")
  }

  if (length(groups) != length(weights_list)) {
    stop("'groups' and 'weights_list' must have the same length.")
  }

  if (!is.logical(point) || length(point) != 1) {
    stop("'point' must be TRUE or FALSE.")
  }

  if (!is.numeric(n_draws) ||
      length(n_draws) != 1 ||
      n_draws <= 0 ||
      n_draws %% 1 != 0) {

    stop("'n_draws' must be a positive integer.")
  }

  if (!is.numeric(n_sims) ||
      length(n_sims) != 1 ||
      n_sims <= 0 ||
      n_sims %% 1 != 0) {

    stop("'n_sims' must be a positive integer.")
  }

  ## z grid
  if (!is.null(z_values)) {

    if (!is.numeric(z_values)) {
      stop("'z_values' must be numeric.")
    }

    if (any(z_values < 0 | z_values > 1)) {
      stop("'z_values' must lie in [0, 1].")
    }
  }

  ## -----------------------------
  ## Per-level checks
  ## -----------------------------

  for (x in seq_along(groups)) {

    node_names  <- groups[[x]]
    alpha_input <- alpha_list[[x]]
    beta_input  <- beta_list[[x]]
    weights     <- weights_list[[x]]

    ## weights length
    if (length(weights) != length(node_names)) {

      stop(
        paste0(
          "Level ", x,
          ": length(weights) must equal number of nodes."
        )
      )
    }

    ## Positive weights
    if (any(weights <= 0)) {

      stop(
        paste0(
          "Level ", x,
          ": weights must all be positive."
        )
      )
    }

    ## Point estimate case
    if (point) {

      if (!is.numeric(alpha_input) || !is.vector(alpha_input)) {

        stop(
          paste0(
            "Level ", x,
            ": alpha_input must be a numeric vector when point = TRUE."
          )
        )
      }

      if (!is.numeric(beta_input) || !is.vector(beta_input)) {

        stop(
          paste0(
            "Level ", x,
            ": beta_input must be a numeric vector when point = TRUE."
          )
        )
      }

      if (length(alpha_input) != length(node_names)) {

        stop(
          paste0(
            "Level ", x,
            ": alpha vector length must equal number of nodes."
          )
        )
      }

      if (length(beta_input) != length(node_names)) {

        stop(
          paste0(
            "Level ", x,
            ": beta vector length must equal number of nodes."
          )
        )
      }

    } else {

      ## Posterior matrix case

      if (!is.matrix(alpha_input)) {

        stop(
          paste0(
            "Level ", x,
            ": alpha_input must be a matrix when point = FALSE."
          )
        )
      }

      if (!is.matrix(beta_input)) {
        stop(
          paste0(
            "Level ", x,
            ": beta_input must be a matrix when point = FALSE."
          )
        )
      }

      if (!all(dim(alpha_input) == dim(beta_input))) {

        stop(
          paste0(
            "Level ", x,
            ": alpha and beta matrices must have identical dimensions."
          )
        )
      }

      if (ncol(alpha_input) != length(node_names)) {

        stop(
          paste0(
            "Level ", x,
            ": number of matrix columns must equal number of nodes."
          )
        )
      }
    }

    ## Positive beta parameters
    if (any(alpha_input <= 0, na.rm = TRUE)) {

      stop(
        paste0(
          "Level ", x,
          ": alpha parameters must all be positive."
        )
      )
    }

    if (any(beta_input <= 0, na.rm = TRUE)) {

      stop(
        paste0(
          "Level ", x,
          ": beta parameters must all be positive."
        )
      )
    }
  }

  invisible(TRUE)
}
