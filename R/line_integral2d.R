#' Line integral of a scalar field in \eqn{\mathbb{R}^2} with 3D visualization
#'
#' Numerically computes \eqn{\int_a^b f(\mathbf r(t))\,\|\mathbf r'(t)\|\,dt} and draws:
#' \itemize{
#'   \item The surface \eqn{z = f(x,y)} over a rectangle defined by the
#'         \emph{endpoints} of the curve (\eqn{x\in[x_a,x_b]}, \eqn{y\in[y_a,y_b]}).
#'   \item The curve \eqn{\mathbf r(t)} in the plane \eqn{z=0} and its “lifted”
#'         version \eqn{(x(t),y(t),f(x(t),y(t)))}.
#'   \item A vertical \strong{curtain} over the curve:
#'         \eqn{(x(t),y(t),v\,f(x(t),y(t)))}, \eqn{v\in[0,1]}.
#' }
#'
#' Derivatives are computed with centered finite differences. The surface domain
#' mirrors the WxMaxima behavior: it uses only the curve \emph{endpoints} \code{r(a)} and \code{r(b)}.
#'
#' @param f Scalar field \code{function(x,y)} returning a numeric \eqn{f(x,y)}.
#' @param r Parametric curve \code{function(t)} returning \code{c(x,y)}.
#' @param a,b Parameter limits \eqn{t} with \code{b > a}.
#' @param plot Logical. If \code{TRUE}, plot with \pkg{plotly}.
#' @param n_curve Number of samples along the curve (for curtain and traces).
#' @param n_curtain_v Vertical resolution of the curtain (\eqn{v}).
#' @param n_surf_x,n_surf_y Grid resolution for the surface \eqn{z=f(x,y)}.
#' @param colorscale Colorscale for the surface (Plotly name, a single color,
#'   or a vector like \code{c("white","#2a9d8f")}).
#' @param surface_opacity Surface opacity in \code{[0,1]}.
#' @param show_surface_grid,surface_grid_color,surface_grid_width Grid overlay on the surface.
#' @param curtain Curtain opacity in \code{[0,1]} (0 disables the curtain).
#' @param curtain_colorscale Colorscale for the curtain (same formats as \code{colorscale}).
#' @param curve_color,curve_width Style for the curve in \eqn{z=0}.
#' @param lifted_color,lifted_width Style for the lifted curve.
#' @param h Step for finite differences in \eqn{t}. If \code{NULL}, it is chosen automatically.
#' @param method Integration method: \code{"adaptive"} uses \code{stats::integrate};
#'   \code{"simpson"} uses composite Simpson with \code{n_simpson} subintervals (even).
#' @param n_simpson Number of subintervals for Simpson (adjusted to be even if necessary).
#' @param scene,bg 3D scene and background settings.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{value}: numeric value of the line integral,
#'   \item \code{fig}: a \pkg{plotly} object if \code{plot=TRUE}, otherwise \code{NULL}.
#' }
#'
#' @examples
#' \dontshow{if (interactive()) \{}
#' # Scalar field and a planar “spiral-ish” curve
#' f <- function(x,y) x^2 + y^2
#' r <- function(t) c(t*cos(t), t*sin(3*t))
#' line_integral2d(
#'   f, r, a = 0, b = 2*pi,
#'   plot = TRUE,
#'   colorscale = c("white","#a8dadc"), surface_opacity = 0.35,
#'   curtain = 0.4, curtain_colorscale = c("white","#1d3557"),
#'   curve_color = "black", lifted_color = "red"
#' )
#' \dontshow{\}}
#'
#' @export
line_integral2d <- function(
    f, r, a, b,
    plot = TRUE,
    n_curve = 400,
    n_curtain_v = 40,
    n_surf_x = 80, n_surf_y = 80,
    colorscale = "Viridis",
    surface_opacity = 0.65,
    show_surface_grid  = TRUE,
    surface_grid_color = "rgba(80,80,80,0.25)",
    surface_grid_width = 1,
    curtain = 0.4,
    curtain_colorscale = c("white","#2a9d8f"),
    curve_color = "black",
    curve_width = 3,
    lifted_color = "red",
    lifted_width = 2,
    h = NULL,
    method = c("adaptive","simpson"),
    n_simpson = 1000,
    scene = list(
      aspectmode = "data",
      xaxis = list(title = "x"),
      yaxis = list(title = "y"),
      zaxis = list(title = "z")
    ),
    bg = list(paper = "white", plot = "white")
) {
  method <- match.arg(method)
  # --- checks
  if (!is.function(f)) stop("'f' must be function(x,y) -> scalar.", call. = FALSE)
  if (!is.function(r)) stop("'r' must be function(t) -> c(x,y).", call. = FALSE)
  if (!is.numeric(a) || !is.numeric(b) || length(a)!=1L || length(b)!=1L ||
      !is.finite(a) || !is.finite(b) || b <= a)
    stop("'a' and 'b' must be scalars with b > a.", call. = FALSE)

  # --- color helpers (same style used elsewhere in the package)
  is_color <- function(x) {
    if (!is.character(x) || length(x) != 1) return(FALSE)
    if (grepl("^rgba?\\(", x, ignore.case = TRUE)) return(TRUE)
    if (grepl("^#([0-9A-Fa-f]{6}|[0-9A-Fa-f]{8})$", x)) return(TRUE)
    ok <- TRUE; tryCatch(grDevices::col2rgb(x), error = function(...) ok <<- FALSE); ok
  }
  as_rgba <- function(col, alpha = NULL) {
    if (grepl("^rgba?\\(", col, ignore.case = TRUE)) return(col)
    rgb <- grDevices::col2rgb(col) / 255
    a_  <- if (is.null(alpha)) 1 else max(0, min(1, alpha))
    sprintf("rgba(%g,%g,%g,%g)", 255*rgb[1], 255*rgb[2], 255*rgb[3], a_)
  }
  as_colorscale <- function(x, alpha = NULL) {
    if (is.null(x)) return(NULL)
    if (is.list(x) && length(x) >= 2 && is.numeric(x[[1]][[1]])) return(x)
    if (is.character(x) && length(x) == 1) {
      if (is_color(x)) {
        ccol <- as_rgba(x, alpha); return(list(list(0, ccol), list(1, ccol)))
      } else return(x)
    }
    if (is.character(x) && length(x) > 1) {
      cols <- vapply(x, as_rgba, character(1), alpha = alpha)
      pos  <- seq(0, 1, length.out = length(cols))
      return(lapply(seq_along(cols), function(i) list(pos[i], cols[i])))
    }
    stop("Unrecognized 'colorscale' format.", call. = FALSE)
  }

  # --- numerical derivative of r(t)
  if (is.null(h)) h <- (b - a) * 1e-4 + 1e-6
  rprime <- function(t) {
    # protect boundaries
    tt1 <- pmax(a, pmin(b, t - h))
    tt2 <- pmax(a, pmin(b, t + h))
    (r(tt2) - r(tt1)) / (tt2 - tt1)
  }
  speed <- function(t) {
    v <- rprime(t); sqrt(sum(v*v))
  }

  # --- integrand: f(r(t)) * ||r'(t)||
  integrand <- function(t) {
    p <- r(t); fx <- f(p[1], p[2])
    fx * speed(t)
  }
  integrand_vec <- Vectorize(integrand, "t")

  # --- integral
  value <- if (method == "adaptive") {
    stats::integrate(integrand_vec, lower = a, upper = b, rel.tol = 1e-6)$value
  } else {
    n <- as.integer(n_simpson); if (n %% 2L == 1L) n <- n + 1L
    tt <- seq(a, b, length.out = n + 1L)
    yy <- integrand_vec(tt)
    hS <- (b - a) / n
    hS * (yy[1] + yy[length(yy)] + 4 * sum(yy[seq(2, n, by = 2)]) +
            2 * sum(yy[seq(3, n-1, by = 2)])) / 3
  }

  # --- plotting data
  fig <- NULL
  if (isTRUE(plot)) {
    if (!requireNamespace("plotly", quietly = TRUE)) {
      warning("To plot you need 'plotly' installed.", call. = FALSE)
    } else {
      # sampled curve
      t_curve <- seq(a, b, length.out = n_curve)
      Rxy <- t(vapply(t_curve, r, numeric(2)))
      cx <- Rxy[,1]; cy <- Rxy[,2]; cz0 <- rep(0, length(t_curve))
      cf <- vapply(seq_along(t_curve), function(i) f(cx[i], cy[i]), numeric(1))
      cz <- cf

      # surface z=f(x,y) over rectangle [r(a), r(b)] (as in Maxima)
      ra <- r(a); rb <- r(b)
      xs <- seq(ra[1], rb[1], length.out = n_surf_x)
      ys <- seq(ra[2], rb[2], length.out = n_surf_y)
      X <- matrix(rep(xs, each = n_surf_y), nrow = n_surf_y)
      Y <- matrix(rep(ys, times = n_surf_x), nrow = n_surf_y)
      Z <- matrix(NA_real_, nrow = n_surf_y, ncol = n_surf_x)
      for (j in seq_len(n_surf_x)) {
        Z[, j] <- vapply(ys, function(y) f(xs[j], y), numeric(1))
      }

      # curtain: (r(u), v * f(r(u)))  u in [a,b], v in [0,1]
      plt <- plotly::plot_ly()
      if (is.numeric(curtain) && curtain > 0) {
        vseq <- seq(0, 1, length.out = n_curtain_v)
        Xc <- matrix(rep(cx, each = n_curtain_v), nrow = n_curtain_v)
        Yc <- matrix(rep(cy, each = n_curtain_v), nrow = n_curtain_v)
        Zc <- matrix(rep(cz, each = n_curtain_v), nrow = n_curtain_v) *
          matrix(rep(vseq, times = n_curve), nrow = n_curtain_v)
        plt <- plt |>
          plotly::add_surface(
            x = Xc, y = Yc, z = Zc,
            colorscale = as_colorscale(curtain_colorscale),
            showscale  = FALSE,
            opacity    = curtain
          )
      }

      # surface z=f(x,y)
      contours_arg <- if (isTRUE(show_surface_grid)) list(
        x = list(show = TRUE, color = surface_grid_color, width = surface_grid_width),
        y = list(show = TRUE, color = surface_grid_color, width = surface_grid_width),
        z = list(show = FALSE)
      ) else NULL

      plt <- plt |>
        plotly::add_surface(
          x = X, y = Y, z = Z,
          colorscale = as_colorscale(colorscale),
          showscale  = FALSE,
          opacity    = surface_opacity,
          contours   = contours_arg
        )

      # curve at z=0 and lifted curve at z=f
      plt <- plt |>
        plotly::add_trace(
          x = cx, y = cy, z = cz0,
          type = "scatter3d", mode = "lines",
          line = list(color = curve_color, width = curve_width),
          showlegend = FALSE
        ) |>
        plotly::add_trace(
          x = cx, y = cy, z = cz,
          type = "scatter3d", mode = "lines",
          line = list(color = lifted_color, width = lifted_width),
          showlegend = FALSE
        ) |>
        plotly::layout(
          title = sprintf("Line integral (value \u2248 %.6g)", value),
          scene = scene,
          paper_bgcolor = bg$paper,
          plot_bgcolor  = bg$plot
        )

      fig <- plt
      print(fig)
    }
  }

  list(value = value, fig = fig)
}
