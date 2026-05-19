#' Validate inputs for Gamma convolution functions
#'
#' Internal helper used to validate inputs supplied to Gamma convolution
#' density estimation functions using S + shape_mat + rate_mat structure.
#'
#' @param S Binary summation matrix defining hierarchy structure.
#'
#' @param shape_mat Numeric vector or matrix of Gamma shape parameters.
#'
#' @param rate_mat Numeric vector or matrix of Gamma rate parameters.
#'
#' @param point Logical indicating whether point estimates are supplied.
#'
#' @param n_draws Number of posterior draws (or NULL in posterior mode).
#'
#' @param n_sims Number of Monte Carlo simulations per draw.
#'
#' @param z_values Optional numeric vector of evaluation points.
#'
#' @return Invisibly returns TRUE if all checks pass.
#'
#' @keywords internal
#' @noRd

check_Gamma_convolution_inputs <- function(
    S,
    shape_mat,
    rate_mat,
    point,
    n_draws,
    n_sims,
    z_values
) {

  ## -----------------------------
  ## Global checks
  ## -----------------------------

  if (!is.matrix(S) && !is.data.frame(S)) {
    stop("'S' must be a binary matrix or data.frame.")
  }

  if (!is.logical(point) || length(point) != 1) {
    stop("'point' must be TRUE or FALSE.")
  }

  if (!isTRUE(point) && !is.null(n_draws)) {
    stop("n_draws must be NULL when point = FALSE (auto-inferred from matrix).")
  }

  if (isTRUE(point)) {

    if (!is.numeric(n_draws) || length(n_draws) != 1 ||
        is.na(n_draws) || !is.finite(n_draws) ||
        n_draws <= 0 || n_draws %% 1 != 0) {
      stop("'n_draws' must be a positive integer.")
    }
  }

  if (!is.numeric(n_sims) || length(n_sims) != 1 ||
      n_sims <= 0 || n_sims %% 1 != 0) {
    stop("'n_sims' must be a positive integer.")
  }

  if (!is.numeric(z_values)) {
    stop("'z_values' must be numeric.")
  }

  if (any(z_values < 0, na.rm = TRUE)) {
    stop("'z_values' must be greater than 0.")
  }

  ## -----------------------------
  ## Node structure checks
  ## -----------------------------

  bottom_nodes <- colnames(S)

  if (is.null(bottom_nodes)) {
    stop("'S' must have column names corresponding to bottom nodes.")
  }

  n_nodes <- length(bottom_nodes)

  ## -----------------------------
  ## Parameter checks
  ## -----------------------------

  if (isTRUE(point)) {

    if (!is.numeric(shape_mat) || !is.vector(shape_mat)) {
      stop("shape_mat must be a numeric vector when point = TRUE.")
    }

    if (!is.numeric(rate_mat) || !is.vector(rate_mat)) {
      stop("rate_mat must be a numeric vector when point = TRUE.")
    }

    if (length(shape_mat) != nrow(S)) {
      stop("Length of shape_mat must equal number of nodes in S.")
    }

    if (length(rate_mat) != nrow(S)) {
      stop("Length of rate_mat must equal number of nodes in S.")
    }

  } else {

    if (!is.matrix(shape_mat)) {
      stop("shape_mat must be a matrix when point = FALSE.")
    }

    if (!is.matrix(rate_mat)) {
      stop("rate_mat must be a matrix when point = FALSE.")
    }

    if (!all(dim(shape_mat) == dim(rate_mat))) {
      stop("shape_mat and rate_mat must have identical dimensions.")
    }

    if (ncol(shape_mat) != nrow(S)) {
      stop("Number of columns in shape_mat must equal number of nodes in S.")
    }

    if (ncol(rate_mat) != nrow(S)) {
      stop("Number of columns in rate_mat must equal number of nodes in S.")
    }
  }

  ## -----------------------------
  ## Positivity checks
  ## -----------------------------

  if (any(shape_mat <= 0, na.rm = TRUE)) {
    stop("All shape parameters must be positive.")
  }

  if (any(rate_mat <= 0, na.rm = TRUE)) {
    stop("All rate parameters must be positive.")
  }

  invisible(TRUE)
}
