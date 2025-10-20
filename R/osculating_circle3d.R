#' Osculating discs (or circles) of a spatial curve r(t)
#'
#' For each \code{t0} in \code{t_points}, it computes the Frenet frame
#' \eqn{\mathbf{T}(t_0), \mathbf{N}(t_0), \mathbf{B}(t_0)} and the curvature
#' \eqn{\kappa(t_0) = \| \mathbf{r}'(t_0) \times \mathbf{r}''(t_0) \| \,/\, \| \mathbf{r}'(t_0) \|^3 }.
#' The osculating circle has center \eqn{\mathbf{C} = \mathbf{r}(t_0) + \mathbf{N}(t_0)/\kappa(t_0)}
#' and radius \eqn{R = 1/\kappa(t_0)}.
#'
#' You can draw the \strong{osculating disc} (\code{fill = "disk"}) or just the
#' \strong{osculating circumference} (\code{fill = "ring"}). A parametrization is:
#' \deqn{\mathbf{S}(u,v) = \mathbf{C} + R\,u\,\big(-\mathbf{N}(t_0)\cos v + \mathbf{T}(t_0)\sin v\big).}
#'
#' @param X,Y,Z Functions of \code{t} returning \code{x(t)}, \code{y(t)}, \code{z(t)}.
#' @param a,b Endpoints of the interval \code{[a,b]}.
#' @param t_points Numeric vector of \code{t} where osculating circles/discs are constructed.
#' @param h Step size for centered finite differences.
#' @param plot \code{TRUE}/\code{FALSE}. If \code{TRUE}, draw with \pkg{plotly}.
#' @param n_samples Number of samples of the base curve (only if \code{show_curve = TRUE}).
#' @param fill \code{"disk"} (filled surface) or \code{"ring"} (only circumference).
#' @param ru Radial subdivisions (only if \code{fill = "disk"}).
#' @param rv Angular subdivisions (also used as ring points).
#' @param colorscale Plotly colorscale for the discs (e.g., \code{"Reds"}).
#' @param opacity Opacity of the discs (0–1, only if \code{fill = "disk"}).
#' @param ring_line Line style for the ring (only if \code{fill = "ring"}).
#' @param show_curve,show_points Whether to show the base curve and the points \code{r(t)}.
#' @param curve_line,point_marker Styles for curve and points (Plotly).
#' @param show_radius \code{TRUE}/\code{FALSE} to draw a radius from \code{C} to the boundary.
#' @param radius_phase Angle (in radians) of the radius (default \code{0}).
#' @param radius_line Line style for the radius, e.g., \code{list(color = "orange", width = 5)}.
#' @param scene 3D scene settings (Plotly).
#' @param bg Background colors (\code{paper}, \code{plot}).
#' @param lighting Lighting options for \code{add_surface} (when \code{fill = "disk"}).
#' @param tol Numerical tolerance.
#'
#' @return
#' \describe{
#'   \item{\code{data}}{A \code{tibble} with columns
#'     \code{t, x, y, z, kappa, cx, cy, cz, radius, Tx, Ty, Tz, Nx, Ny, Nz, Bx, By, Bz}.}
#'   \item{\code{plot}}{(if \code{plot = TRUE}) a \pkg{plotly} object.}
#' }
#'
#' @examples
#' X <- function(t) cos(t); Y <- function(t) sin(t); Z <- function(t) 0.2*t
#' osculating_circle3d(X, Y, Z, a = 0, b = 6*pi, t_points = c(pi, 2*pi), plot = FALSE)
#'
#' @export
osculating_circle3d <- function(
    X, Y, Z,
    a, b,
    t_points,
    h = 1e-4,
    plot = FALSE,
    n_samples = 400,
    fill = c("disk","ring"),
    ru = 24,
    rv = 72,
    colorscale = "Reds",
    opacity = 0.6,
    ring_line = list(color = "red", width = 4, dash = "solid"),
    show_curve = TRUE,
    show_points = TRUE,
    curve_line = list(color = "blue", width = 2, dash = "solid"),
    point_marker = list(color = "black", size = 3, symbol = "circle"),
    show_radius = FALSE,
    radius_phase = 0,
    radius_line = list(color = "orange", width = 5, dash = "solid"),
    scene = list(
      aspectmode = "data",
      xaxis = list(title = "x(t)"),
      yaxis = list(title = "y(t)"),
      zaxis = list(title = "z(t)")
    ),
    bg = list(paper = "white", plot = "white"),
    lighting = list(ambient = 1, diffuse = 0.15, specular = 0, roughness = 1, fresnel = 0),
    tol = 1e-10
) {
  fill <- match.arg(fill)

  if (!is.numeric(t_points) || any(!is.finite(t_points)))
    stop("'t_points' must be finite numeric values.", call. = FALSE)
  if (any(t_points < a | t_points > b))
    stop("All 't_points' must lie within [a,b].", call. = FALSE)
  if (b < a) stop("'b' must be >= 'a'.", call. = FALSE)

  d1 <- function(f, t, h) (f(t + h) - f(t - h)) / (2*h)
  d2 <- function(f, t, h) (f(t + h) - 2*f(t) + f(t - h)) / (h*h)
  r  <- function(t) c(X(t), Y(t), Z(t))
  dot   <- function(a,b) sum(a*b)
  norm  <- function(a) sqrt(dot(a,a))
  cross <- function(a,b) c(a[2]*b[3]-a[3]*b[2],
                           a[3]*b[1]-a[1]*b[3],
                           a[1]*b[2]-a[2]*b[1])

  compute_one <- function(t0) {
    r0 <- r(t0)
    r1 <- c(d1(X, t0, h), d1(Y, t0, h), d1(Z, t0, h))
    r2 <- c(d2(X, t0, h), d2(Y, t0, h), d2(Z, t0, h))
    n1 <- norm(r1)
    if (n1 < tol) {
      return(c(t=t0, r0, kappa=NA_real_, rep(NA_real_, 3), radius=NA_real_,
               rep(NA_real_, 9)))
    }
    T <- r1 / n1
    c12 <- cross(r1, r2); n12 <- norm(c12)
    if (n12 < tol) {
      return(c(t=t0, r0, kappa=0, rep(NA_real_, 3), radius=Inf,
               T, rep(NA_real_, 6)))
    }
    B <- c12 / n12
    N <- cross(B, T); N <- N / norm(N)
    kappa <- n12 / (n1^3)
    C <- r0 + N / kappa
    R <- 1 / kappa
    c(t=t0, r0, kappa=kappa, C, radius=R, T, N, B)
  }

  M <- t(vapply(t_points, compute_one, numeric(1 + 3 + 1 + 3 + 1 + 3 + 3 + 3)))
  colnames(M) <- c("t","x","y","z","kappa","cx","cy","cz","radius",
                   "Tx","Ty","Tz","Nx","Ny","Nz","Bx","By","Bz")
  out <- tibble::as_tibble(M)

  if (isTRUE(plot)) {
    if (!requireNamespace("plotly", quietly = TRUE)) {
      warning("You need 'plotly' installed for plotting.", call. = FALSE)
    } else {
      curve_df <- NULL
      if (isTRUE(show_curve)) {
        ts <- seq(a, b, length.out = n_samples)
        curve_df <- tibble::tibble(x = X(ts), y = Y(ts), z = Z(ts))
      }
      plt <- if (!is.null(curve_df)) {
        plotly::plot_ly(curve_df, x=~x, y=~y, z=~z,
                        type="scatter3d", mode="lines",
                        line=curve_line, hoverinfo="none", showlegend=FALSE)
      } else plotly::plot_ly()

      if (isTRUE(show_points)) {
        plt <- plt |>
          plotly::add_trace(
            x = out$x, y = out$y, z = out$z,
            type = "scatter3d", mode = "markers",
            marker = point_marker, showlegend = FALSE, hoverinfo = "none"
          )
      }

      V <- seq(0, 2*pi, length.out = rv)
      cosV <- cos(V); sinV <- sin(V)

      for (i in seq_len(nrow(out))) {
        if (!is.finite(out$kappa[i]) || out$kappa[i] <= tol) next
        C  <- c(out$cx[i], out$cy[i], out$cz[i])
        R  <- out$radius[i]
        T  <- c(out$Tx[i], out$Ty[i], out$Tz[i])
        N  <- c(out$Nx[i], out$Ny[i], out$Nz[i])

        if (fill == "disk") {
          U <- seq(0, 1, length.out = ru)
          Xmat <- matrix(NA_real_, ru, rv)
          Ymat <- matrix(NA_real_, ru, rv)
          Zmat <- matrix(NA_real_, ru, rv)
          for (j in seq_len(rv)) {
            dir <- -N * cosV[j] + T * sinV[j]
            Xmat[, j] <- C[1] + R * U * dir[1]
            Ymat[, j] <- C[2] + R * U * dir[2]
            Zmat[, j] <- C[3] + R * U * dir[3]
          }
          plt <- plt |>
            plotly::add_surface(
              x = Xmat, y = Ymat, z = Zmat,
              colorscale = colorscale, showscale = FALSE,
              opacity = opacity, lighting = lighting
            )
        } else { # ring
          x_ring <- C[1] + R * (-N[1] * cosV + T[1] * sinV)
          y_ring <- C[2] + R * (-N[2] * cosV + T[2] * sinV)
          z_ring <- C[3] + R * (-N[3] * cosV + T[3] * sinV)
          x_ring <- c(x_ring, x_ring[1])
          y_ring <- c(y_ring, y_ring[1])
          z_ring <- c(z_ring, z_ring[1])

          plt <- plt |>
            plotly::add_trace(
              x = x_ring, y = y_ring, z = z_ring,
              type = "scatter3d", mode = "lines",
              line = ring_line, showlegend = FALSE, hoverinfo = "none"
            )
        }

        if (isTRUE(show_radius)) {
          v0 <- radius_phase
          dir0 <- -N * cos(v0) + T * sin(v0)
          p_end <- C + R * dir0
          plt <- plt |>
            plotly::add_trace(
              x = c(C[1], p_end[1]),
              y = c(C[2], p_end[2]),
              z = c(C[3], p_end[3]),
              type = "scatter3d", mode = "lines",
              line = radius_line, showlegend = FALSE, hoverinfo = "none"
            )
        }
      }

      plt <- plt |>
        plotly::layout(
          title = if (fill == "disk") "Osculating Discs" else "Osculating Circles",
          scene = scene,
          paper_bgcolor = bg$paper,
          plot_bgcolor  = bg$plot
        )
      print(plt)
    }
  }

  out
}
