#' Validate inputs for Gamma convolution functions
#'
#' Internal helper used to validate inputs supplied to the Gamma convolution density estimation functions.
#'
#' @param groups A list of node names at each hierarchical level.
#' @param shape_list A list of shape parameter vectors or matrices.
#' @param rate_list A list of rate parameter vectors or matrices.
#' @param point Logical indicating whether point estimates are supplied.
#' @param n_draws Positive integer giving the number of draws.
#' @param n_sims Positive integer giving the number of simulations.
#' @param z_values Optional numeric vector of evaluation points.
#'
#' @return Invisibly returns TRUE if all checks pass.
#'
#' @keywords internal
#' @noRd

check_Gamma_convolution_inputs <- function(groups, shape_list, rate_list, point, n_draws, n_sims, z_values) {

  ## -----------------------------
  ## Global checks
  ## -----------------------------
  if (!isTRUE(point) && !is.null(n_draws)) {
    stop("n_draws must be NULL when point = FALSE. The function automatically uses the number of posterior draws.")
  }

  if (!is.list(groups)) {
    stop("'groups' must be a list.")
  }

  if (!is.list(shape_list)) {
    stop("'shape_list' must be a list.")
  }

  if (!is.list(rate_list)) {
    stop("'rate_list' must be a list.")
  }

  if (length(groups) != length(shape_list)) {
    stop("'groups' and 'shape_list' must have the same length.")
  }

  if (length(groups) != length(rate_list)) {
    stop("'groups' and 'rate_list' must have the same length.")
  }

  if (!is.logical(point) || length(point) != 1) {
    stop("'point' must be TRUE or FALSE.")
  }

  if (isTRUE(point)) {
    if (!is.numeric(n_draws) || length(n_draws) != 1 ||
        is.na(n_draws) || !is.finite(n_draws) ||
        n_draws <= 0 || n_draws %% 1 != 0) {
      stop("'n_draws' must be a positive integer.")
    }
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

    if (any(z_values < 0)) {
      stop("'z_values' must be greater than 0.")
    }
  }

  ## -----------------------------
  ## Per-level checks
  ## -----------------------------

  for (x in seq_along(groups)) {

    node_names  <- groups[[x]]
    shape_input <- shape_list[[x]]
    rate_input  <- rate_list[[x]]

    ## Point estimate case
    if (point) {

      if (!is.numeric(shape_input) || !is.vector(shape_input)) {

        stop(
          paste0(
            "Level ", x,
            ": shape_input must be a numeric vector when point = TRUE."
          )
        )
      }

      if (!is.numeric(rate_input) || !is.vector(rate_input)) {

        stop(
          paste0(
            "Level ", x,
            ": rate_input must be a numeric vector when point = TRUE."
          )
        )
      }

      if (length(shape_input) != length(node_names)) {

        stop(
          paste0(
            "Level ", x,
            ": alpha vector length must equal number of nodes."
          )
        )
      }

      if (length(rate_input) != length(node_names)) {

        stop(
          paste0(
            "Level ", x,
            ": beta vector length must equal number of nodes."
          )
        )
      }

    } else {

      ## Posterior matrix case

      if (!is.matrix(shape_input)) {

        stop(
          paste0(
            "Level ", x,
            ": shape_input must be a matrix when point = FALSE."
          )
        )
      }

      if (!is.matrix(rate_input)) {
        stop(
          paste0(
            "Level ", x,
            ": rate_input must be a matrix when point = FALSE."
          )
        )
      }

      if (!all(dim(shape_input) == dim(rate_input))) {

        stop(
          paste0(
            "Level ", x,
            ": shape and rate matrices must have identical dimensions."
          )
        )
      }

      if (ncol(shape_input) != length(node_names)) {

        stop(
          paste0(
            "Level ", x,
            ": number of matrix columns must equal number of nodes."
          )
        )
      }
    }

    ## Positive beta parameters
    if (any(shape_input <= 0, na.rm = TRUE)) {

      stop(
        paste0(
          "Level ", x,
          ": shape parameters must all be positive."
        )
      )
    }

    if (any(rate_input <= 0, na.rm = TRUE)) {

      stop(
        paste0(
          "Level ", x,
          ": rate parameters must all be positive."
        )
      )
    }
  }

  invisible(TRUE)
}
