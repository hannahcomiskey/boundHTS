#' Tilt a Base Density to Match a Target Mean
#'
#' @param mu_theory Numeric. The target mean to match with the tilted density.
#' @param y_vals Numeric vector. Support points of the base density.
#' @param f_y Numeric vector. Base density values corresponding to `y_vals`.
#' @param z_values Numeric vector. Points at which the tilted density is interpolated. Defaults to `y_vals`.
#' @param discrete Logical. If TRUE, calculate density of discrete random variable. Defaults to FALSE.
#'
#' @return A list with components:
#' \describe{
#'   \item{f_tilted}{Numeric vector. Tilted and normalized density evaluated at `z_values`.}
#'   \item{tilted_samps}{Numeric vector. Random samples drawn from the tilted density.}
#'   \item{tilting_parameter}{Numeric. Estimated exponential tilting parameter used.}
#' }
#'
#' @details
#' The function searches for the tilting parameter `nu` that satisfies the
#' moment condition using a grid search followed by `uniroot` to refine the
#' solution. If no sign change is detected in the moment condition, no tilting
#' is applied (i.e., `nu = 0`). The method is suitable for discrete or bounded
#' densities and can be used for post-hoc adjustment in hierarchical
#' probabilistic forecasting.
#'
#' @examples
#' Example: tilt a simple discrete density
#' y <- 0:5
#' f <- dbinom(y, size = 5, prob = 0.3) # E(Y) = 5*0.3 = 1.5
#' mu_target <- 2.5
#' res <- tilt_density(mu_theory = mu_target,
#'                     y_vals = y,
#'                     f_y = f, discrete=TRUE)
#' plot(y, f, col = "blue", type='h',
#' main = "Original vs Tilted Density",
#' xlab = "y", ylab = "Density")
#' lines(y, res$f_tilted, col = "red", type='h')
#' legend("topright", legend = c("Original", "Tilted"),
#' col = c("blue", "red"), lwd = 2)
#'
#' Example: tilt a simple continuous density
#' y <- seq(0, 1, length.out = 200)
#' f <- dbeta(y, shape1 = 2, shape2 = 5)   # E(Y) = 2/(2+5) = 0.2857
#' mu_target <- 0.5
#' res <- tilt_density(mu_theory = mu_target,
#'                     y_vals = y,
#'                     f_y = f,
#'                     discrete = FALSE)
#' plot(y, f, type = "l", col = "blue", lwd = 2,
#' main = "Original vs Tilted Density",
#' xlab = "y", ylab = "Density")
#' lines(y, res$f_tilted, col = "red", lwd = 2)
#' legend("topright", legend = c("Original", "Tilted"),
#'         col = c("blue", "red"), lwd = 2)
#'
#' @export

tilt_density <- function(mu_theory, y_vals, f_y, discrete=FALSE) {

  # Extract and normalise convolution base density
  if(discrete==FALSE) {
    f_y <- f_y / pracma::trapz(y_vals, f_y)
  } else{
    f_y <- f_y / sum(f_y)
  }

  # Root finding bracket search
  nu_grid <- seq(-2000, 2000, length.out = 4000)
  vals <- sapply(nu_grid, moment_condition_tilting,
                 f_y = f_y,
                 y_vals = y_vals,
                 mu_theory = mu_theory)

  idx <- which(diff(sign(vals)) != 0)

  if (length(idx) == 0) {
    warning("No sign change for node — using no tilting (nu=0).")
    nu_star <- 0
  } else {
    lower <- nu_grid[idx[1]]
    upper <- nu_grid[idx[1] + 1]
    nu_star <- stats::uniroot(moment_condition_tilting,
                              lower = lower,
                              upper = upper,
                              f_y = f_y,
                              y_vals = y_vals,
                              mu_theory = mu_theory)$root
  }

  # Tilted density
  if(discrete==FALSE){
    f_tilt_i <- tilted_density_cont(nu_star, f_y, y_vals) # tilt for continuous RV
    f_tilt_i <- f_tilt_i / pracma::trapz(y_vals, f_tilt_i) # normalise
    f_tilted <- stats::approx(y_vals, f_tilt_i, xout = y_vals, rule = 2)$y
  } else{
    f_tilt_i <- tilted_density_discrete(nu_star, f_y, y_vals) # tilt for discrete RV
    f_tilted <- f_tilt_i / sum(f_tilt_i) # normalise
  }

  # Sample from tilted density
  tilted_samps <- sample(y_vals, 5000, replace = TRUE, prob = f_tilted)

  return(list(nu_star = nu_star,
              f_tilted = f_tilted,
              tilted_samps = tilted_samps
              )
         )
}
