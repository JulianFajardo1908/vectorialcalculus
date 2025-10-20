#' Plot a smooth parametric surface r(u,v) and estimate its area
#'
#' @description
#' \deqn{\text{Area} \approx \iint \left\| \mathbf{r}_u(u,v) \times \mathbf{r}_v(u,v) \right\| \, du \, dv}
#' via the composite trapezoid rule.
#'
#' @details
#' \deqn{\mathbf{r}(u,v) = \big(x(u,v),\, y(u,v),\, z(u,v)\big).}
#'
#' @param xfun,yfun,zfun Functions of \code{(u, v)} returning numeric scalars.
#' @param urange,vrange Numeric length-2 ranges for \code{u} and \code{v}.
#' @param nu,nv Integers \eqn{\ge 20}, grid sizes along \code{u} and \code{v}.
#' @param h_u,h_v Steps for finite differences in \code{u} and \code{v} (central differences).
#' @param title_prefix Character used in plot title (default \code{"r(u,v)"}).
#'
#' @return
#' \itemize{
#'   \item \code{plot}: plotly object.
#'   \item \code{area}: numeric estimate of the surface area.
#'   \item \code{grid}: list with \code{U}, \code{V}, \code{X}, \code{Y}, \code{Z}.
#' }
#'
#' @export
surface_parametric_area <- function(
    xfun, yfun, zfun,
    urange = c(0, 2*pi),
    vrange = c(0, 2*pi),
    nu = 160, nv = 160,
    h_u = NULL, h_v = NULL,
    title_prefix = "r(u,v)"
) {
  stopifnot(is.function(xfun), is.function(yfun), is.function(zfun))
  stopifnot(length(urange) == 2, length(vrange) == 2, nu >= 20, nv >= 20)

  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("Package 'plotly' is required. Install with: install.packages('plotly')")
  }

  # Parameter grid
  u <- seq(urange[1], urange[2], length.out = nu)
  v <- seq(vrange[1], vrange[2], length.out = nv)
  U <- matrix(rep(u, each = nv), nrow = nv, ncol = nu)
  V <- matrix(rep(v, times = nu), nrow = nv, ncol = nu)

  # Vectorized coord functions
  xv <- Vectorize(function(uu, vv) xfun(uu, vv))
  yv <- Vectorize(function(uu, vv) yfun(uu, vv))
  zv <- Vectorize(function(uu, vv) zfun(uu, vv))

  # Evaluate and FORCE matrix shape (nv x nu)
  X <- matrix(xv(U, V), nrow = nv, ncol = nu)
  Y <- matrix(yv(U, V), nrow = nv, ncol = nu)
  Z <- matrix(zv(U, V), nrow = nv, ncol = nu)

  # Finite-difference steps
  if (is.null(h_u)) h_u <- max(1e-6, (urange[2] - urange[1]) / (nu - 1))
  if (is.null(h_v)) h_v <- max(1e-6, (vrange[2] - vrange[1]) / (nv - 1))

  # Safe derivatives on matrices
  d_du <- function(M, delta) {
    if (is.null(dim(M))) stop("Expected a matrix in d_du; got a vector.")
    out <- M * 0
    nc <- ncol(M)
    if (nc >= 3) out[, 2:(nc-1)] <- (M[, 3:nc] - M[, 1:(nc-2)]) / (2 * delta)
    out[, 1] <- (M[, 2] - M[, 1]) / delta
    out[, nc] <- (M[, nc] - M[, nc-1]) / delta
    out
  }
  d_dv <- function(M, delta) {
    if (is.null(dim(M))) stop("Expected a matrix in d_dv; got a vector.")
    out <- M * 0
    nr <- nrow(M)
    if (nr >= 3) out[2:(nr-1), ] <- (M[3:nr, ] - M[1:(nr-2), ]) / (2 * delta)
    out[1, ] <- (M[2, ] - M[1, ]) / delta
    out[nr, ] <- (M[nr, ] - M[nr-1, ]) / delta
    out
  }

  # r_u and r_v
  Xu <- d_du(X, h_u); Yu <- d_du(Y, h_u); Zu <- d_du(Z, h_u)
  Xv <- d_dv(X, h_v); Yv <- d_dv(Y, h_v); Zv <- d_dv(Z, h_v)

  # Cross product and Jacobian magnitude
  Cx <- Yu * Zv - Zu * Yv
  Cy <- Zu * Xv - Xu * Zv
  Cz <- Xu * Yv - Yu * Xv
  mag <- sqrt(Cx^2 + Cy^2 + Cz^2)

  # Composite trapezoid weights
  w_u <- rep(1, nu); w_u[c(1, nu)] <- 0.5
  w_v <- rep(1, nv); w_v[c(1, nv)] <- 0.5
  W <- outer(w_v, w_u)
  du <- (urange[2] - urange[1]) / (nu - 1)
  dv <- (vrange[2] - vrange[1]) / (nv - 1)

  area_est <- sum(W * mag) * du * dv

  # Plot
  p <- plotly::plot_ly() |>
    plotly::add_surface(x = ~X, y = ~Y, z = ~Z, showscale = FALSE, name = "Surface") |>
    plotly::layout(
      title = paste0("Parametric surface area (estimate): ", signif(area_est, 8)),
      scene = list(xaxis = list(title = "x"),
                   yaxis = list(title = "y"),
                   zaxis = list(title = "z"))
    )

  message(sprintf("Estimated area \u2248 %.10g", area_est))
  list(plot = p, area = area_est, grid = list(U = U, V = V, X = X, Y = Y, Z = Z))
}

