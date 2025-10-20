# Surface integral over an oriented surface z = g(x,y) with optional plot
# Supports:
#  - Scalar surface integral: \u222C_S phi dS
#  - Flux (vector surface integral): \u222C_S F \u00B7 n dS
# Orientation: "up" (default) or "down"
#
# Numeric method: composite trapezoid rule on a uniform (x,y) grid.
#
# Arguments:
#  gfun     : function(x, y) -> z
#  xlim,ylim: numeric length-2 ranges for x and y
#  nx, ny   : grid sizes (>= 20 recommended)
#  scalar_phi: optional function(x, y, z) -> numeric (for \u222C phi dS)
#  vector_F : optional function(x, y, z) -> c(Fx, Fy, Fz) (for \u222C F\u00B7n dS)
#             If both provided, both integrals are computed.
#  orientation: "up" or "down" (affects the vector/flux integral)
#  plot     : logical; if TRUE, returns a plotly surface colored by the integrand
#  title    : character title prefix for plots
#
# Returns: list with fields (area_density_integral, flux_integral, plot, grids, fields)
surface_integral_z <- function(
    gfun,
    xlim, ylim,
    nx = 160, ny = 160,
    scalar_phi = NULL,
    vector_F  = NULL,
    orientation = c("up", "down"),
    plot = TRUE,
    title = "Surface integral over z=g(x,y)"
) {
  stopifnot(is.function(gfun), length(xlim) == 2, length(ylim) == 2, nx >= 20, ny >= 20)
  orientation <- match.arg(orientation)

  if (plot && !requireNamespace("plotly", quietly = TRUE)) {
    stop("Package 'plotly' is required for plotting. Install with: install.packages('plotly')")
  }

  # Grid in (x,y)
  x <- seq(xlim[1], xlim[2], length.out = nx)
  y <- seq(ylim[1], ylim[2], length.out = ny)
  X <- matrix(rep(x, each = ny), nrow = ny, ncol = nx)
  Y <- matrix(rep(y, times = nx), nrow = ny, ncol = nx)
  # Evaluate surface z = g(x,y)
  gv <- Vectorize(function(xx, yy) gfun(xx, yy))
  Z <- matrix(gv(X, Y), nrow = ny, ncol = nx)

  # Finite differences for gx, gy (central inside, one-sided on edges)
  d_dx <- function(M, dx) {
    out <- M * 0
    nc <- ncol(M)
    if (nc >= 3) out[, 2:(nc-1)] <- (M[, 3:nc] - M[, 1:(nc-2)]) / (2 * dx)
    out[, 1]  <- (M[, 2] - M[, 1]) / dx
    out[, nc] <- (M[, nc] - M[, nc-1]) / dx
    out
  }
  d_dy <- function(M, dy) {
    out <- M * 0
    nr <- nrow(M)
    if (nr >= 3) out[2:(nr-1), ] <- (M[3:nr, ] - M[1:(nr-2), ]) / (2 * dy)
    out[1, ]  <- (M[2, ] - M[1, ]) / dy
    out[nr, ] <- (M[nr, ] - M[nr-1, ]) / dy
    out
  }

  dx <- (xlim[2] - xlim[1]) / (nx - 1)
  dy <- (ylim[2] - ylim[1]) / (ny - 1)
  gx <- d_dx(Z, dx)
  gy <- d_dy(Z, dy)

  # Composite trapezoid weights on the (x,y) grid
  wx <- rep(1, nx); wx[c(1, nx)] <- 0.5
  wy <- rep(1, ny); wy[c(1, ny)] <- 0.5
  W  <- outer(wy, wx) # ny x nx

  # 1) Scalar surface integral: \u222C_S phi dS = \u222C_D phi(x,y,g) * sqrt(1+gx^2+gy^2) dx dy
  area_density_integral <- NA_real_
  scalar_integrand <- NULL
  if (is.function(scalar_phi)) {
    phiv <- Vectorize(function(xx, yy, zz) scalar_phi(xx, yy, zz))
    PHI  <- matrix(phiv(X, Y, Z), nrow = ny, ncol = nx)
    J    <- sqrt(1 + gx^2 + gy^2)
    scalar_integrand <- PHI * J
    area_density_integral <- sum(W * scalar_integrand) * dx * dy
  }

  # 2) Flux (vector surface integral): \u222C_S F\u00B7n dS = \u222C_D F \u00B7 N dx dy, with
  #    N = (-gx, -gy, 1) for upward; N = (gx, gy, -1) for downward (sign flip)
  flux_integral <- NA_real_
  flux_integrand <- NULL
  if (is.function(vector_F)) {
    Fmat <- array(0, dim = c(ny, nx, 3))
    Fv <- function(xx, yy, zz) {
      val <- vector_F(xx, yy, zz)
      if (length(val) != 3) stop("vector_F must return numeric length-3 c(Fx,Fy,Fz).")
      val
    }
    # Evaluate F component-wise
    Fx <- matrix(0, ny, nx); Fy <- Fx; Fz <- Fx
    for (j in 1:nx) for (i in 1:ny) {
      vv <- Fv(X[i, j], Y[i, j], Z[i, j])
      Fx[i, j] <- vv[1]; Fy[i, j] <- vv[2]; Fz[i, j] <- vv[3]
    }
    if (orientation == "up") {
      Nx <- -gx; Ny <- -gy; Nz <-  1
    } else { # "down"
      Nx <-  gx; Ny <-  gy; Nz <- -1
    }
    flux_integrand <- Fx * Nx + Fy * Ny + Fz * Nz
    flux_integral  <- sum(W * flux_integrand) * dx * dy
  }

  # Optional plot: color by whichever integrand is present (phi*J first, else flux)
  plt <- NULL
  if (plot) {
    col_mat <- if (!is.null(scalar_integrand)) scalar_integrand else flux_integrand
    if (is.null(col_mat)) col_mat <- Z * 0
    plt <- plotly::plot_ly() |>
      plotly::add_surface(x = ~X, y = ~Y, z = ~Z, surfacecolor = ~col_mat,
                          colorbar = list(title = if (!is.null(scalar_integrand)) "phi\u00B7|ru\u00D7rv|" else "F\u00B7N"),
                          name = "Surface") |>
      plotly::layout(
        title = if (!is.na(area_density_integral)) {
          paste0(title,
                 sprintf(" | scalar \u222C \u2248 %.6g", area_density_integral),
                 if (!is.na(flux_integral)) sprintf(" | flux \u222C \u2248 %.6g", flux_integral) else "")
        } else if (!is.na(flux_integral)) {
          paste0(title, sprintf(" | flux \u222C \u2248 %.6g", flux_integral))
        } else title,
        scene = list(
          xaxis = list(title = "x"),
          yaxis = list(title = "y"),
          zaxis = list(title = "z")
        )
      )
  }

  list(
    area_density_integral = area_density_integral,
    flux_integral = flux_integral,
    plot = plt,
    grids = list(X = X, Y = Y, Z = Z, gx = gx, gy = gy, W = W, dx = dx, dy = dy),
    fields = list(scalar_integrand = scalar_integrand, flux_integrand = flux_integrand)
  )
}
