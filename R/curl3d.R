#' Numerical curl of a 3D vector field
#'
#' Computes the curl \eqn{\nabla \times \mathbf F} of a vector field
#' \eqn{\mathbf F(x,y,z)} (or \eqn{\mathbf F(x,y,z,t)}) at a given point using
#' second-order \strong{central finite differences}.
#'
#' For \eqn{\mathbf F=(F_x,F_y,F_z)}, the curl is
#' \deqn{
#' \nabla \times \mathbf F =
#' \left(
#'   \frac{\partial F_z}{\partial y} - \frac{\partial F_y}{\partial z},\;
#'   \frac{\partial F_x}{\partial z} - \frac{\partial F_z}{\partial x},\;
#'   \frac{\partial F_y}{\partial x} - \frac{\partial F_x}{\partial y}
#' \right).
#' }
#'
#' If \code{F} accepts a time argument \code{t}, the curl is evaluated at a fixed time
#' \code{tval}, i.e., \eqn{\nabla\times \mathbf F(\cdot,\cdot,\cdot,t=\texttt{tval})}.
#'
#' @param F Field function: \code{function(x,y,z)} or \code{function(x,y,z,t)} returning
#'   a numeric vector \code{c(Fx, Fy, Fz)}.
#' @param x,y,z Coordinates of the evaluation point.
#' @param h Finite-difference step size(s). One of:
#'   \itemize{
#'     \item Scalar (same step for \eqn{x,y,z});
#'     \item Numeric vector of length 3 \code{c(hx, hy, hz)};
#'     \item \code{NULL} (default): chosen as \code{1e-4 * (1 + abs(coord))} per axis.
#'   }
#' @param tval Time value if \code{F} depends on \code{t}. Default \code{0}.
#' @param method Differencing scheme. Currently \code{"central"} (2nd order).
#'
#' @return Named numeric vector of length 3 with the curl components:
#'   \code{c(omega_x, omega_y, omega_z)}.
#'
#' @examples
#' # Example 1: F(x,y,z) = (-y, x, 0.6)  => curl = (0, 0, 2)
#' F1 <- function(x,y,z) c(-y, x, 0.6)
#' curl3d(F1, x = 0.1, y = -0.3, z = 2)  # ~ c(0, 0, 2)
#'
#' # Example 2 (with t): F(x,y,z,t) = (-y, x + t, z)
#' # curl = (0, 0, 2) (independent of t)
#' F2 <- function(x,y,z,t) c(-y, x + t, z)
#' curl3d(F2, x = 1, y = 2, z = 3, tval = 5)
#'
#' # Tip: for highly varying fields, try smaller steps:
#' curl3d(F1, x=1, y=1, z=1, h = 1e-5)
#'
#' @export
curl3d <- function(F, x, y, z, h = NULL, tval = 0, method = c("central")) {
  method <- match.arg(method)
  if (!is.function(F)) stop("'F' must be function(x,y,z) or function(x,y,z,t).", call. = FALSE)
  if (!is.numeric(x) || !is.numeric(y) || !is.numeric(z) ||
      length(x)!=1L || length(y)!=1L || length(z)!=1L ||
      !is.finite(x) || !is.finite(y) || !is.finite(z)) {
    stop("'x','y','z' must be finite numeric scalars.", call. = FALSE)
  }

  # step sizes
  if (is.null(h)) {
    hx <- 1e-4 * (1 + abs(x))
    hy <- 1e-4 * (1 + abs(y))
    hz <- 1e-4 * (1 + abs(z))
  } else if (length(h) == 1L) {
    hx <- hy <- hz <- as.numeric(h)
  } else if (length(h) == 3L) {
    hx <- as.numeric(h[1]); hy <- as.numeric(h[2]); hz <- as.numeric(h[3])
  } else {
    stop("'h' must be a scalar, a length-3 vector, or NULL.", call. = FALSE)
  }
  if (any(!is.finite(c(hx,hy,hz))) || any(c(hx,hy,hz) <= 0)) {
    stop("Step sizes 'h' must be positive and finite.", call. = FALSE)
  }

  # wrapper: evaluate F with or without time
  F_eval <- function(xx, yy, zz, tt) {
    out <- if (length(formals(F)) >= 4L) F(xx, yy, zz, tt) else F(xx, yy, zz)
    if (!is.numeric(out) || length(out) != 3L || any(!is.finite(out))) {
      stop("F(x,y,z[,t]) must return a finite numeric vector of length 3.", call. = FALSE)
    }
    out
  }

  # central differences (2nd order)
  # dFz/dy
  Fp <- F_eval(x, y + hy, z, tval)[3]; Fm <- F_eval(x, y - hy, z, tval)[3]
  dFz_dy <- (Fp - Fm) / (2 * hy)

  # dFy/dz
  Fp <- F_eval(x, y, z + hz, tval)[2]; Fm <- F_eval(x, y, z - hz, tval)[2]
  dFy_dz <- (Fp - Fm) / (2 * hz)

  # dFx/dz
  Fp <- F_eval(x, y, z + hz, tval)[1]; Fm <- F_eval(x, y, z - hz, tval)[1]
  dFx_dz <- (Fp - Fm) / (2 * hz)

  # dFz/dx
  Fp <- F_eval(x + hx, y, z, tval)[3]; Fm <- F_eval(x - hx, y, z, tval)[3]
  dFz_dx <- (Fp - Fm) / (2 * hx)

  # dFy/dx
  Fp <- F_eval(x + hx, y, z, tval)[2]; Fm <- F_eval(x - hx, y, z, tval)[2]
  dFy_dx <- (Fp - Fm) / (2 * hx)

  # dFx/dy
  Fp <- F_eval(x, y + hy, z, tval)[1]; Fm <- F_eval(x, y - hy, z, tval)[1]
  dFx_dy <- (Fp - Fm) / (2 * hy)

  # curl
  omega_x <- dFz_dy - dFy_dz
  omega_y <- dFx_dz - dFz_dx
  omega_z <- dFy_dx - dFx_dy
  c(omega_x = omega_x, omega_y = omega_y, omega_z = omega_z)
}
