#' Validate inputs for Poisson convolution functions
#'
#' Internal helper function used to validate inputs supplied to
#' Poisson convolution density estimation functions.
#'
#' @param groups A named list of nodes at each hierarchical level.
#' @param z_values Numeric vector of evaluation points.
#' @param lambda_list A list of lambda vectors or matrices.
#' @param point Logical indicating whether point estimates are supplied.
#'
#' @return Invisibly returns TRUE if all checks pass.
#'
#' @keywords internal
#' @noRd
check_poisson_convolution_inputs <- function(groups,
                                             z_values,
                                             lambda_list,
                                             point) {

  ## ---------------------------------------------------------------
  ## groups checks
  ## ---------------------------------------------------------------

  if (!is.list(groups)) {
    stop(
      "'groups' must be a list.",
      call. = FALSE
    )
  }

  if (length(groups) == 0L) {
    stop(
      "'groups' must contain at least one hierarchy level.",
      call. = FALSE
    )
  }

  if (is.null(names(groups)) || any(names(groups) == "")) {
    stop(
      "'groups' must be a named list.",
      call. = FALSE
    )
  }

  ## check node names
  for (x in seq_along(groups)) {

    node_names <- groups[[x]]

    if (!is.character(node_names)) {
      stop(
        paste0(
          "groups[[", x, "]] must be a character vector."
        ),
        call. = FALSE
      )
    }

    if (length(node_names) == 0L) {
      stop(
        paste0(
          "groups[[", x, "]] must contain at least one node."
        ),
        call. = FALSE
      )
    }

    if (anyNA(node_names)) {
      stop(
        paste0(
          "groups[[", x, "]] contains missing node names."
        ),
        call. = FALSE
      )
    }
  }

  ## ---------------------------------------------------------------
  ## z_values checks
  ## ---------------------------------------------------------------

  if (!is.numeric(z_values)) {
    stop(
      "'z_values' must be numeric.",
      call. = FALSE
    )
  }

  if (length(z_values) == 0L) {
    stop(
      "'z_values' must contain at least one value.",
      call. = FALSE
    )
  }

  if (anyNA(z_values)) {
    stop(
      "'z_values' must not contain missing values.",
      call. = FALSE
    )
  }

  if (any(z_values < 0)) {
    stop(
      "'z_values' must be non-negative.",
      call. = FALSE
    )
  }

  ## ---------------------------------------------------------------
  ## lambda_list checks
  ## ---------------------------------------------------------------

  if (!is.list(lambda_list)) {
    stop(
      "'lambda_list' must be a list.",
      call. = FALSE
    )
  }

  if (length(lambda_list) != length(groups)) {
    stop(
      "'lambda_list' must have the same length as 'groups'.",
      call. = FALSE
    )
  }

  ## ---------------------------------------------------------------
  ## point checks
  ## ---------------------------------------------------------------

  if (!is.logical(point) || length(point) != 1L) {
    stop(
      "'point' must be a single logical value.",
      call. = FALSE
    )
  }

  ## ---------------------------------------------------------------
  ## Per-level checks
  ## ---------------------------------------------------------------

  for (x in seq_along(groups)) {

    node_names  <- groups[[x]]
    lambda_input <- lambda_list[[x]]

    ## ------------------------------------------------------------
    ## Point-estimate case
    ## ------------------------------------------------------------

    if (isTRUE(point)) {

      if (!is.numeric(lambda_input) || !is.vector(lambda_input)) {
        stop(
          paste0(
            "lambda_list[[", x,
            "]] must be a numeric vector when point = TRUE."
          ),
          call. = FALSE
        )
      }

      if (length(lambda_input) != length(node_names)) {
        stop(
          paste0(
            "lambda_list[[", x,
            "]] must have length equal to the number of nodes."
          ),
          call. = FALSE
        )
      }
    }

    ## ------------------------------------------------------------
    ## Posterior-sample case
    ## ------------------------------------------------------------

    if (!isTRUE(point)) {

      if (!is.matrix(lambda_input)) {
        stop(
          paste0(
            "lambda_list[[", x,
            "]] must be a matrix when point = FALSE."
          ),
          call. = FALSE
        )
      }

      if (!is.numeric(lambda_input)) {
        stop(
          paste0(
            "lambda_list[[", x,
            "]] must be numeric."
          ),
          call. = FALSE
        )
      }

      if (ncol(lambda_input) != length(node_names)) {
        stop(
          paste0(
            "The number of columns in lambda_list[[", x,
            "]] must equal the number of nodes."
          ),
          call. = FALSE
        )
      }
    }

    ## ------------------------------------------------------------
    ## Shared checks
    ## ------------------------------------------------------------

    if (anyNA(lambda_input)) {
      stop(
        paste0(
          "lambda_list[[", x,
          "]] contains missing values."
        ),
        call. = FALSE
      )
    }

    if (any(lambda_input <= 0)) {
      stop(
        paste0(
          "lambda_list[[", x,
          "]] must contain positive values only."
        ),
        call. = FALSE
      )
    }
  }

  invisible(TRUE)
}
