#' Directional derivative and local visualization on a surface z = F(x,y)
#'
#' Numerical translation of \code{der_dir(...)} from WxMaxima.
#' Given \eqn{F:\mathbb{R}^2\to\mathbb{R}}, a point \eqn{(a,b)}, and a direction
#' \eqn{\mathbf d=(d_x,d_y)}, this computes \eqn{D=\nabla F(a,b)\cdot \hat{\mathbf d}} and,
#' optionally, draws a strip of the surface near the line along \eqn{\hat{\mathbf d}},
#' the directional curve on \eqn{F}, and the vectors \eqn{\mathbf u}, \eqn{\mathbf w}.
#'
#' @param F Function \code{F(x,y)} returning \code{z}.
#' @param point Numeric vector \code{c(a,b)}.
#' @param dir In-plane direction vector \code{c(dx,dy)} (non-zero).
#' @param h Step size for centered finite differences (default \code{1e-4}).
#' @param plot Logical; if \code{TRUE}, draw with \pkg{plotly}.
#' @param x_window,y_window Half-widths of the local window in \code{x} and \code{y}.
#' @param n_s,n_r Strip resolutions (parameters analogous to \code{x} and \code{r} in Maxima).
#' @param show_strip \code{TRUE}/\code{FALSE}, draw the strip over the surface (instead of the whole surface).
#' @param strip_colorscale Plotly colorscale for the strip (e.g. \code{"Blues"}).
#' @param strip_opacity Strip opacity (0–1).
#' @param show_surface_grid \code{TRUE}/\code{FALSE} to show a grid over the strip.
#' @param surface_grid_color Grid color (e.g. \code{"rgba(60,80,200,0.25)"}).
#' @param surface_grid_width Grid line width.
#' @param curve_line Line style for the directional curve (e.g. \code{list(color="red", width=0.8)}).
#' @param point_marker Style of the base point \eqn{(a,b,f(a,b))}.
#' @param u_line,w_line Line styles for the vectors \eqn{\mathbf u} and \eqn{\mathbf w}.
#' @param t_range Range of \code{t} for the directional curve (default \code{c(-2,2)}).
#' @param scene 3D scene settings (default \code{aspectmode="data"}).
#' @param bg Background colors: \code{list(paper="white", plot="white")}.
#' @param tol Numerical tolerance.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{fx}, \code{fy}: numerical gradient at \code{(a,b)}.
#'   \item \code{D}: directional derivative \eqn{D=\nabla F\cdot \hat d}.
#'   \item \code{dir_hat}: unit direction \eqn{\hat d}.
#'   \item \code{u}, \code{w}: 3D vectors used in the figure.
#'   \item \code{fig}: plotly object if \code{plot=TRUE}, otherwise \code{NULL}.
#' }
#'
#' @examples
#' F <- function(x,y) (x^2 + y^2)/8 + sin(x/2)
#' res <- directional_derivative3d(
#'   F, point = c(1, -1), dir = c(1, 0.5),
#'   plot = TRUE, strip_colorscale = "Blues",
#'   show_surface_grid = TRUE
#' )
#' res$D
#'
#' @export
directional_derivative3d <- function(
    F,
    point,
    dir,
    h = 1e-4,
    plot = FALSE,
    x_window = 2, y_window = 2,
    n_s = 180, n_r = 50,
    show_strip = TRUE,
    strip_colorscale = "Blues",
    strip_opacity = 0.30,
    show_surface_grid = TRUE,
    surface_grid_color = "rgba(60,80,200,0.25)",
    surface_grid_width = 1,
    curve_line = list(color = "red",   width = 0.8),
    point_marker = list(color = "black", size = 3, symbol = "circle"),
    u_line = list(color = "black", width = 1.0),
    w_line = list(color = "black", width = 0.5),
    t_range = c(-2, 2),
    scene = list(aspectmode = "data",
                 xaxis = list(title = "x"),
#' @noRd
                 yaxis = list(title = "y"),
                 zaxis = list(title = "z")),
    bg = list(paper = "white", plot = "white"),
    tol = 1e-12
#' @noRd
) {
  if (!is.numeric(point) || length(point) != 2L || any(!is.finite(point)))
    stop("'point' must be a finite numeric vector c(a,b).", call. = FALSE)
  if (!is.numeric(dir) || length(dir) != 2L || any(!is.finite(dir)))
    stop("'dir' must be a finite numeric vector c(dx,dy).", call. = FALSE)

  a <- point[1]; b <- point[2]
  dx <- dir[1]; dy <- dir[2]
  nd <- sqrt(dx*dx + dy*dy); if (nd < tol) stop("'dir' cannot be the zero vector.", call. = FALSE)
  d_hat <- c(dx, dy) / nd

  # Centered partial derivatives
  fx <- (F(a + h, b) - F(a - h, b)) / (2*h)
  fy <- (F(a, b + h) - F(a, b - h)) / (2*h)

  # Directional derivative
  D  <- fx * d_hat[1] + fy * d_hat[2]

  # Vectors u and w (as in your original script)
  u3 <- c(d_hat[1], d_hat[2], D)
  nu3 <- sqrt(sum(u3*u3)); if (nu3 < tol) nu3 <- 1
  u <- u3 / nu3
  w <- c(u[1], u[2], 0)

  fig <- NULL
  if (isTRUE(plot)) {
    if (!requireNamespace("plotly", quietly = TRUE)) {
      warning("To plot you need 'plotly' installed.", call. = FALSE)
    } else {
      xs_line <- seq(a - x_window, a + x_window, length.out = n_s)

      # Build a narrow surface strip along the chosen direction
      build_strip <- function() {
        # line through (a,b) with slope dy/dx if dx != 0
        if (abs(d_hat[1]) >= 1e-14) {
          slope <- d_hat[2] / d_hat[1]
          s_seq <- xs_line - a
          r_seq <- seq(0, 1, length.out = n_r)
          S <- matrix(rep(s_seq, each = n_r), nrow = n_r)
          R <- matrix(rep(r_seq, times = length(s_seq)), nrow = n_r)
          Xmat <- matrix(rep(xs_line, each = n_r), nrow = n_r)
          Ymat <- R * (b + y_window) + (1 - R) * (slope * S + b)
        } else {
          # "vertical" direction: vary y and mix x between a and a + x_window
          y_line <- seq(b - y_window, b + y_window, length.out = n_s)
          r_seq  <- seq(0, 1, length.out = n_r)
          Ymat <- matrix(rep(y_line, each = n_r), nrow = n_r)
          Xmat <- matrix(
            rep((r_seq * (a + x_window) + (1 - r_seq) * a), times = length(y_line)),
            nrow = n_r
          )
        }
        Zmat <- matrix(mapply(F, as.numeric(Xmat), as.numeric(Ymat)),
                       nrow = nrow(Xmat), ncol = ncol(Xmat))
        list(X=Xmat, Y=Ymat, Z=Zmat)
      }

      plt <- plotly::plot_ly()

      if (isTRUE(show_strip)) {
        strip <- build_strip()
        contours_arg <- if (isTRUE(show_surface_grid)) list(
          x = list(show = TRUE, color = surface_grid_color, width = surface_grid_width),
          y = list(show = TRUE, color = surface_grid_color, width = surface_grid_width),
          z = list(show = FALSE)
        ) else NULL

        plt <- plt |>
          plotly::add_surface(
            x = strip$X, y = strip$Y, z = strip$Z,
            colorscale = strip_colorscale,
            showscale  = FALSE,
            opacity    = strip_opacity,
            contours   = contours_arg
          )
      }

      # Base point and vectors u, w (as segments)
      f0 <- F(a, b)
      p0 <- c(a, b, f0)

      plt <- plt |>
        plotly::add_trace(
          x = p0[1], y = p0[2], z = p0[3],
          type = "scatter3d", mode = "markers",
          marker = point_marker, showlegend = FALSE, hoverinfo = "none"
        )

      add_seg <- function(p, v, line) {
        plotly::add_trace(
          plt,
          x = c(p[1], p[1] + v[1]), y = c(p[2], p[2] + v[2]), z = c(p[3], p[3] + v[3]),
          type = "scatter3d", mode = "lines",
          line = line, showlegend = FALSE, hoverinfo = "none"
        )
      }
      # Visual segment lengths
      L_u <- 1.0; L_w <- 1.0
      plt <- add_seg(p0, L_u * u, u_line)
      plt <- add_seg(p0, L_w * w, w_line)

      # Directional curve (parametric line through (a,b) along u[1:2])
      ts <- seq(t_range[1], t_range[2], length.out = 400)
      x_dir <- a + u[1]*ts
      y_dir <- b + u[2]*ts
      z_dir <- vapply(seq_along(ts), function(i) F(x_dir[i], y_dir[i]), numeric(1))

      plt <- plt |>
        plotly::add_trace(
          x = x_dir, y = y_dir, z = z_dir,
          type = "scatter3d", mode = "lines",
          line = curve_line, showlegend = FALSE, hoverinfo = "none"
        )

      # Layout
      ttl <- sprintf("Directional derivative D = ∇F(a,b)·d̂ = %.6g", D)
      plt <- plt |>
        plotly::layout(
          title = ttl,
          scene = scene,
          paper_bgcolor = bg$paper,
          plot_bgcolor  = bg$plot
        )

      print(plt)
      fig <- plt
    }
  }

  list(
    fx = fx, fy = fy,
    D  = D,
    dir_hat = d_hat,
    u = u, w = w,
    fig = fig
  )
}
