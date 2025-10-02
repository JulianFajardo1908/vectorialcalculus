#' Principal normal vectors N along a 3D curve (numeric)
#'
#' Computes the \strong{principal normal} \eqn{N(t)} at parameter values
#' \code{t_points}. With \eqn{r'(t)} and \eqn{r''(t)}, we form
#' \eqn{T = r'/\|r'\|}, then \eqn{B = \mathrm{normalize}(r' \times r'')}, and finally
#' \eqn{N = \mathrm{normalize}(B \times T)}. If \(\kappa \approx 0\),
#' \eqn{N} is undefined (returns \code{NA}).
#'
#' @param X Function \code{x(t)}.
#' @param Y Function \code{y(t)}.
#' @param Z Function \code{z(t)}.
#' @param a Lower endpoint of \code{[a,b]}.
#' @param b Upper endpoint of \code{[a,b]}.
#' @param t_points Numeric vector of \code{t} where \eqn{N} is evaluated/plotted.
#' @param h Step size for centered finite differences.
#' @param plot Logical; if \code{TRUE}, plots with \pkg{plotly}.
#' @param n_samples Number of points to sample the curve (plotting only).
#' @param vec_scale Base scale for \eqn{N} segments. If \code{NULL}, estimated (5% of span).
#' @param vec_factor Multiplicative factor applied to \code{vec_scale} (default \code{1}).
#' @param curve_line Plotly style for the curve.
#' @param N_line Plotly style for the \eqn{N} segments.
#' @param show_curve Logical; show the base curve.
#' @param show_points Logical; mark \eqn{r(t)} at \code{t_points}.
#' @param point_marker Plotly marker style for base points.
#' @param scene Plotly 3D scene settings.
#' @param bg Background colors for canvas (\code{paper}, \code{plot}).
#' @param tol Tolerance for singularities and \(\kappa \approx 0\).
#'
#' @return A \code{tibble} with columns \code{t, x, y, z, Nx, Ny, Nz}.
#'
#' @examples
#' X <- function(t) t*cos(t); Y <- function(t) t*sin(3*t); Z <- function(t) t
#' Ntab <- normal3d(X, Y, Z, a = 0, b = 2*pi, t_points = c(pi/3, pi, 5*pi/3))
#' # \donttest{ if (requireNamespace("plotly", quietly = TRUE)) {
#' #   normal3d(X, Y, Z, 0, 2*pi, t_points = c(pi/3, pi, 5*pi/3),
#' #            plot = TRUE, n_samples = 300,
#' #            vec_factor = 1.4,
#' #            N_line = list(color="green", width=5))
#' # } }
#' @noRd

.assert_fun <- function(f, nm) {
  if (!is.function(f)) stop(sprintf("'%s' must be a function of a single variable t.", nm), call. = FALSE)
}
.assert_numeric_scalar <- function(x, nm) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x)) {
    stop(sprintf("'%s' must be a finite numeric scalar.", nm), call. = FALSE)
  }
}

normal3d <- function(
    X, Y, Z,
    a, b,
    t_points,
    h = 1e-4,
    plot = FALSE,
    n_samples = 400,
    vec_scale = NULL,
    vec_factor = 1,
    curve_line = list(color = "blue",  width = 2, dash = "solid"),
    N_line    = list(color = "green", width = 5, dash = "solid"),
    show_curve  = TRUE,
    show_points = TRUE,
    point_marker = list(color = "black", size = 3, symbol = "circle"),
    scene = list(
      aspectmode = "data",
      xaxis = list(title = "x(t)"),
      yaxis = list(title = "y(t)"),
#' @noRd
#' @noRd
      zaxis = list(title = "z(t)")
    ),
    bg = list(paper = "white", plot = "white"),
    tol = 1e-10
) {
#' @noRd
  if (b < a) stop("'b' must be >= 'a'.", call. = FALSE)
  if (!is.numeric(t_points) || any(!is.finite(t_points)))
    stop("'t_points' must be finite numeric.", call. = FALSE)
  if (any(t_points < a | t_points > b))
    stop("All 't_points' must lie within [a, b].", call. = FALSE)
  if (vec_factor <= 0) stop("'vec_factor' must be > 0.", call. = FALSE)
  if (!is.numeric(h) || h <= 0) stop("'h' must be positive.", call. = FALSE)

  d1 <- function(f, t, h) (f(t + h) - f(t - h)) / (2*h)
  d2 <- function(f, t, h) (f(t + h) - 2*f(t) + f(t - h)) / (h*h)
  r  <- function(t) c(X(t), Y(t), Z(t))
  dot  <- function(a,b) sum(a*b); norm <- function(a) sqrt(dot(a,a))
  cross <- function(a,b) c(a[2]*b[3]-a[3]*b[2],
                           a[3]*b[1]-a[1]*b[3],
                           a[1]*b[2]-a[2]*b[1])

  compute_one <- function(t0) {
    r1 <- c(d1(X, t0, h), d1(Y, t0, h), d1(Z, t0, h))
    r2 <- c(d2(X, t0, h), d2(Y, t0, h), d2(Z, t0, h))
    n1 <- norm(r1)
    if (n1 < tol) return(c(t = t0, r(t0), Nx = NA_real_, Ny = NA_real_, Nz = NA_real_))
    T  <- r1 / n1
    c12 <- cross(r1, r2); n12 <- norm(c12)
    if (n12 < tol) return(c(t = t0, r(t0), Nx = NA_real_, Ny = NA_real_, Nz = NA_real_))
    B <- c12 / n12
    N <- cross(B, T); N <- N / norm(N)
    c(t = t0, r(t0), N)
  }

  M <- t(vapply(t_points, compute_one, numeric(1 + 3 + 3)))
  colnames(M) <- c("t","x","y","z","Nx","Ny","Nz")
  out <- tibble::as_tibble(M)

  if (isTRUE(plot)) {
    if (!requireNamespace("plotly", quietly = TRUE)) {
      warning("To plot you need 'plotly' installed.", call. = FALSE)
    } else {
      data_curve <- if (isTRUE(show_curve) || is.null(vec_scale))
        curve_sample3d(X, Y, Z, a, b, n_samples = n_samples) else NULL

      if (is.null(vec_scale)) {
        rx <- if (!is.null(data_curve)) diff(range(data_curve$x)) else 1
        ry <- if (!is.null(data_curve)) diff(range(data_curve$y)) else 1
        rz <- if (!is.null(data_curve)) diff(range(data_curve$z)) else 1
        vec_scale <- 0.05 * max(rx, ry, rz)
      }
      vec_scale <- vec_scale * vec_factor

      make_segment_df <- function(p0, v, s) tibble::tibble(
        x = c(p0[1], p0[1] + s*v[1]),
        y = c(p0[2], p0[2] + s*v[2]),
        z = c(p0[3], p0[3] + s*v[3])
      )

      plt <- if (isTRUE(show_curve) && !is.null(data_curve)) {
        plotly::plot_ly(
          data_curve, x = ~x, y = ~y, z = ~z,
          type = "scatter3d", mode = "lines",
          line = curve_line, hoverinfo = "none", showlegend = FALSE
        )
      } else plotly::plot_ly()

      if (isTRUE(show_points)) {
        plt <- plt |>
          plotly::add_trace(
            x = out$x, y = out$y, z = out$z,
            type = "scatter3d", mode = "markers",
            marker = point_marker, showlegend = FALSE, hoverinfo = "none"
          )
      }

      for (i in seq_len(nrow(out))) {
        v  <- as.numeric(out[i, c("Nx","Ny","Nz")]); if (any(is.na(v))) next
        p0 <- as.numeric(out[i, c("x","y","z")])
        seg <- make_segment_df(p0, v, vec_scale)
        plt <- plt |>
          plotly::add_trace(
            data = seg, x = ~x, y = ~y, z = ~z,
            type = "scatter3d", mode = "lines",
            line = N_line, hoverinfo = "none", showlegend = FALSE
          )
      }

      plt <- plt |>
        plotly::layout(
          title = "Principal normal vectors N",
          scene = scene, paper_bgcolor = bg$paper, plot_bgcolor = bg$plot
        )
      print(plt)
    }
  }

  out
}
