#' Validate inputs for Beta convolution functions
#'
#' Internal helper used to validate inputs supplied to the
#' Beta convolution density estimation functions.
#'
#' @param S Binary summation matrix with rows corresponding to all
#' hierarchy nodes and columns corresponding to bottom-level nodes.
#'
#' @param W Aggregation weight matrix.
#'
#' @param alpha_mat Numeric vector (point = TRUE) or matrix
#' (point = FALSE) of Beta shape1/alpha parameters for all hierarchy nodes.
#'
#' @param beta_mat Numeric vector (point = TRUE) or matrix
#' (point = FALSE) of Beta shape2/beta parameters for all hierarchy nodes.
#'
#' @param point Logical indicating whether point estimates are supplied.
#'
#' @param n_draws Positive integer giving the number of draws.
#'
#' @param n_sims Positive integer giving the number of simulations.
#'
#' @param z_values Optional numeric vector of evaluation points.
#'
#' @return Invisibly returns TRUE if all checks pass.
#'
#' @keywords internal
#' @noRd

check_beta_convolution_inputs <- function(
    S,
    W,
    alpha_mat,
    beta_mat,
    point,
    n_draws,
    n_sims,
    z_values
) {

  ## --------------------------------------------------------
  ## S matrix checks
  ## --------------------------------------------------------

  if (!is.matrix(S)) {
    stop("'S' must be a matrix.")
  }

  if (!all(S %in% c(0, 1))) {
    stop("'S' must contain only 0/1 entries.")
  }

  if (is.null(rownames(S))) {
    stop("'S' must contain row names.")
  }

  if (is.null(colnames(S))) {
    stop("'S' must contain column names.")
  }

  ## --------------------------------------------------------
  ## W matrix checks
  ## --------------------------------------------------------

  if (!is.matrix(W)) {
    stop("'W' must be a matrix.")
  }

  if (!all(dim(S) == dim(W))) {
    stop("'S' and 'W' must have identical dimensions.")
  }

  if (is.null(rownames(W))) {
    stop("'W' must contain row names.")
  }

  if (is.null(colnames(W))) {
    stop("'W' must contain column names.")
  }

  ## --------------------------------------------------------
  ## Global checks
  ## --------------------------------------------------------

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

  ## --------------------------------------------------------
  ## z grid checks
  ## --------------------------------------------------------

  if (!is.null(z_values)) {

    if (!is.numeric(z_values)) {
      stop("'z_values' must be numeric.")
    }

    if (any(z_values < 0 | z_values > 1)) {
      stop("'z_values' must lie in [0, 1].")
    }
  }

  ## --------------------------------------------------------
  ## Point estimate case
  ## --------------------------------------------------------

  if (isTRUE(point)) {

    if (!is.numeric(alpha_mat) ||
        !is.vector(alpha_mat)) {

      stop(
        paste0(
          "'alpha_mat' must be a numeric vector ",
          "when point = TRUE."
        )
      )
    }

    if (!is.numeric(beta_mat) ||
        !is.vector(beta_mat)) {

      stop(
        paste0(
          "'beta_mat' must be a numeric vector ",
          "when point = TRUE."
        )
      )
    }

    if (length(alpha_mat) != nrow(S)) {

      stop(
        "'alpha_mat' length must equal ",
        "the number of hierarchy nodes."
      )
    }

    if (length(beta_mat) != nrow(S)) {

      stop(
        "'beta_mat' length must equal ",
        "the number of hierarchy nodes."
      )
    }

    if (is.null(names(alpha_mat))) {

      stop(
        "'alpha_mat' must contain names ",
        "matching hierarchy nodes."
      )
    }

    if (is.null(names(beta_mat))) {

      stop(
        "'beta_mat' must contain names ",
        "matching hierarchy nodes."
      )
    }

  } else {

    ## ------------------------------------------------------
    ## Posterior matrix case
    ## ------------------------------------------------------

    if (!is.matrix(alpha_mat)) {

      stop(
        "'alpha_mat' must be a matrix ",
        "when point = FALSE."
      )
    }

    if (!is.matrix(beta_mat)) {

      stop(
        "'beta_mat' must be a matrix ",
        "when point = FALSE."
      )
    }

    if (!all(dim(alpha_mat) == dim(beta_mat))) {

      stop(
        "'alpha_mat' and 'beta_mat' ",
        "must have identical dimensions."
      )
    }

    if (nrow(alpha_mat) != n_draws) {

      stop(
        "Number of rows in 'alpha_mat' ",
        "must equal 'n_draws'."
      )
    }

    if (nrow(beta_mat) != n_draws) {

      stop(
        "Number of rows in 'beta_mat' ",
        "must equal 'n_draws'."
      )
    }

    if (ncol(alpha_mat) != nrow(S)) {

      stop(
        "Number of columns in 'alpha_mat' ",
        "must equal the number of hierarchy nodes."
      )
    }

    if (ncol(beta_mat) != nrow(S)) {

      stop(
        "Number of columns in 'beta_mat' ",
        "must equal the number of hierarchy nodes."
      )
    }

    if (is.null(colnames(alpha_mat))) {

      stop(
        "'alpha_mat' must contain column names ",
        "matching hierarchy nodes."
      )
    }

    if (is.null(colnames(beta_mat))) {

      stop(
        "'beta_mat' must contain column names ",
        "matching hierarchy nodes."
      )
    }
  }

  ## --------------------------------------------------------
  ## Positive Beta parameters
  ## --------------------------------------------------------

  if (any(alpha_mat <= 0, na.rm = TRUE)) {

    stop(
      "'alpha_mat' parameters must all be positive."
    )
  }

  if (any(beta_mat <= 0, na.rm = TRUE)) {

    stop(
      "'beta_mat' parameters must all be positive."
    )
  }

  invisible(TRUE)
}
