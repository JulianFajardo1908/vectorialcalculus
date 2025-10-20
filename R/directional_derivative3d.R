#' Directional derivative (nD) and 2D visualization for z = f(x, y)
#'
#' Computes the directional derivative \eqn{D_{\mathbf{\hat v}} f(\mathbf{x}_0)}
#' via centered differences in any dimension. If \code{length(x0) == 2} and
#' \code{plot = TRUE}, it draws a local strip of the surface \eqn{z = f(x,y)},
#' the base point, the direction, and the directional curve.
#'
#' @param f Function returning a numeric scalar. It may be defined either as
#'   \code{f(x, y, ...)} with as many numeric arguments as dimensions, or as
#'   \code{f(x_vec)} taking a numeric vector.
#' @param x0 Numeric vector; evaluation point (length = dimension, \eqn{\ge 2}).
#' @param v  Numeric vector; direction (same length as \code{x0}, nonzero).
#' @param h  Numeric step for finite differences (default \code{1e-6}).
#' @param plot Logical; if \code{TRUE} and \code{length(x0) == 2}, draw with \pkg{plotly}.
#' @param x_window,y_window Numeric half-widths for the 2D strip window (plot only).
#' @param n_s Integer; samples along the direction line (plot only).
#' @param n_r Integer; radial samples across the strip (plot only).
#' @param show_strip Logical; draw the surface strip (plot only).
#' @param strip_colorscale Character; Plotly colorscale for the strip.
#' @param strip_opacity Numeric in \eqn{[0,1]}; strip opacity.
#' @param show_surface_grid Logical; show strip grid (plot only).
#' @param surface_grid_color Character; grid color (plot only).
#' @param surface_grid_width Numeric; grid width (plot only).
#' @param curve_line,point_marker,u_line,w_line Plotly styles (plot only).
#' @param t_range Numeric length-2; parameter range for the directional curve (plot only).
#' @param scene,bg Plotly scene/background (plot only).
#' @param tol Numeric tolerance for zero checks.
#'
#' @return A list with
#' \describe{
#'   \item{\code{D}}{Directional derivative at \code{x0} along \code{v_hat}.}
#'   \item{\code{v_hat}}{Unit direction.}
#'   \item{\code{fx}, \code{fy}}{(only in 2D) centered partials at \code{x0}.}
#'   \item{\code{fig}}{Plotly object if plotted, else \code{NULL}.}
#' }
#'
#' @examples
#' # Works for nD (no plot):
#' f3 <- function(x, y, z) x^2 + y^2 + z
#' directional_derivative3d(f3, x0 = c(1, 0, 0), v = c(1, 1, 0))
#'
#' # 2D with plot (z = f(x,y)):
#' f2 <- function(x, y) x^2 + y^2
#' # directional_derivative3d(f2, x0 = c(0, 0), v = c(1, 2), plot = TRUE)
#'
#' @export
directional_derivative3d <- function(
    f,
    x0,
    v,
    h = 1e-6,
    plot = FALSE,
    x_window = 2, y_window = 2,
    n_s = 180, n_r = 50,
    show_strip = TRUE,
    strip_colorscale = "Blues",
    strip_opacity = 0.30,
    show_surface_grid = TRUE,
    surface_grid_color = "rgba(60,80,200,0.25)",
    surface_grid_width = 1,
    curve_line   = list(color = "red",   width = 1),
    point_marker = list(color = "black", size = 3, symbol = "circle"),
    u_line = list(color = "black", width = 1.0),
    w_line = list(color = "black", width = 0.5),
    t_range = c(-2, 2),
    scene = list(aspectmode = "data",
                 xaxis = list(title = "x"),
                 yaxis = list(title = "y"),
                 zaxis = list(title = "z")),
    bg = list(paper = "white", plot = "white"),
    tol = 1e-12
) {
  # --- Validaciones base ---
  stopifnot(is.numeric(h), length(h) == 1L, is.finite(h), h > 0)
  if (!is.numeric(x0) || any(!is.finite(x0))) stop("'x0' must be finite numeric.", call. = FALSE)
  if (!is.numeric(v)  || any(!is.finite(v))  || length(v) != length(x0))
    stop("'v' must be finite numeric and same length as 'x0'.", call. = FALSE)
  nv <- sqrt(sum(v * v)); if (nv < tol) stop("'v' cannot be the zero vector.", call. = FALSE)
  v_hat <- as.numeric(v / nv)

  # Evaluador robusto: acepta f(x,y,...) o f(c(...))
  eval_f <- function(xx) {
    xx <- as.numeric(xx)
    out <- tryCatch(do.call(f, as.list(xx)),
                    error = function(e) f(xx))
    if (!is.numeric(out) || length(out) != 1L || !is.finite(out))
      stop("f(...) did not return a finite numeric scalar.", call. = FALSE)
    out
  }

  # Derivada direccional
  g   <- function(t) eval_f(x0 + t * v_hat)
  D   <- (g(h) - g(-h)) / (2 * h)

  fx <- fy <- NA_real_
  fig <- NULL

  # Si es 2D, computa fx,fy y (opcional) plotea
  if (length(x0) == 2L) {
    a <- x0[1]; b <- x0[2]
    fx <- (eval_f(c(a + h, b)) - eval_f(c(a - h, b))) / (2 * h)
    fy <- (eval_f(c(a, b + h)) - eval_f(c(a, b - h))) / (2 * h)

    if (isTRUE(plot)) {
      if (!requireNamespace("plotly", quietly = TRUE)) {
        warning("To plot you need 'plotly' installed.", call. = FALSE)
      } else {
        # Linea base en x (si dirección no es vertical), o en y si lo es
        xs_line <- seq(a - x_window, a + x_window, length.out = n_s)

        build_strip <- function() {
          if (abs(v_hat[1]) >= 1e-14) {
            slope <- v_hat[2] / v_hat[1]
            s_seq <- xs_line - a
            r_seq <- seq(0, 1, length.out = n_r)
            S <- matrix(rep(s_seq, each = n_r), nrow = n_r)
            Xmat <- matrix(rep(xs_line, each = n_r), nrow = n_r)
            Ymat <- r_seq %o% rep(b + y_window, length(xs_line)) +
              (1 - r_seq) %o% (slope * as.numeric(S) + b)
            Ymat <- matrix(Ymat, nrow = n_r)
          } else {
            y_line <- seq(b - y_window, b + y_window, length.out = n_s)
            r_seq  <- seq(0, 1, length.out = n_r)
            Ymat <- matrix(rep(y_line, each = n_r), nrow = n_r)
            Xmat <- matrix(rep(r_seq * (a + x_window) + (1 - r_seq) * a,
                               times = length(y_line)), nrow = n_r)
          }
          Zmat <- matrix(mapply(function(xx, yy) eval_f(c(xx, yy)),
                                as.numeric(Xmat), as.numeric(Ymat)),
                         nrow = nrow(Xmat), ncol = ncol(Xmat))
          list(X = Xmat, Y = Ymat, Z = Zmat)
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
              colorscale = strip_colorscale, showscale = FALSE,
              opacity = strip_opacity, contours = contours_arg
            )
        }

        f0 <- eval_f(c(a, b))
        p0 <- c(a, b, f0)

        # vectores
        u3 <- c(v_hat[1], v_hat[2], D)
        nu3 <- sqrt(sum(u3 * u3)); if (nu3 < tol) nu3 <- 1
        u <- u3 / nu3
        w <- c(u[1], u[2], 0)

        add_seg <- function(plt, p, v, line)
          plotly::add_trace(plt,
                            x = c(p[1], p[1] + v[1]), y = c(p[2], p[2] + v[2]), z = c(p[3], p[3] + v[3]),
                            type = "scatter3d", mode = "lines", line = line,
                            showlegend = FALSE, hoverinfo = "none"
          )

        plt <- plt |>
          plotly::add_trace(
            x = p0[1], y = p0[2], z = p0[3],
            type = "scatter3d", mode = "markers",
            marker = point_marker, showlegend = FALSE, hoverinfo = "none"
          )
        plt <- add_seg(plt, p0, u, u_line)
        plt <- add_seg(plt, p0, w, w_line)

        # Curva direccional (x(t), y(t), z(t) = f(x(t), y(t)))
        ts <- seq(t_range[1], t_range[2], length.out = 400)
        x_dir <- a + v_hat[1] * ts
        y_dir <- b + v_hat[2] * ts
        z_dir <- vapply(seq_along(ts), function(i) eval_f(c(x_dir[i], y_dir[i])), numeric(1))

        plt <- plt |>
          plotly::add_trace(
            x = x_dir, y = y_dir, z = z_dir,
            type = "scatter3d", mode = "lines",
            line = curve_line, showlegend = FALSE, hoverinfo = "none"
          ) |>
          plotly::layout(
            title = sprintf("Directional derivative D = %.6g", D),
            scene = scene, paper_bgcolor = bg$paper, plot_bgcolor = bg$plot
          )

        fig <- plt
        print(plt)
      }
    }
  }

  out <- list(D = D, v_hat = v_hat, fx = fx, fy = fy, fig = fig)
  class(out) <- c("directional_derivative", class(out))
  out
}
