# -------------------------------------------------------------------
#' Curvature and torsion of a 3D parametric curve (numeric)
#'
#' Computes at a point \code{t0} the \strong{curvature}
#' \deqn{\kappa(t0) = \frac{\lVert r'(t0) \times r''(t0)\rVert}{\lVert r'(t0)\rVert^3}}
#' and the \strong{torsion}
#' \deqn{\tau(t0) = \frac{\det(r'(t0), r''(t0), r'''(t0))}{\lVert r'(t0) \times r''(t0)\rVert^2}}
#' using centered finite differences for derivatives (1st, 2nd and 3rd).
#' Optionally, it can plot a \emph{local segment} of the curve around \code{t0}
#' and mark the point \code{r(t0)} in 3D with \pkg{plotly}.
#'
#' @param X,Y,Z Functions of \code{t} returning scalars (coordinates of \eqn{r(t)}).
#' @param t0 Parameter value at which curvature and torsion are evaluated.
#' @param h Step size for finite differences (recommended between \code{1e-6} and \code{1e-3}, depending on scale).
#' @param plot Logical. If \code{TRUE}, shows a curve segment around \code{t0}.
#' @param window Window length in parameter \code{t} for plotting (interval \code{[t0 - window/2, t0 + window/2]}).
#' @param n_samples Number of sample points in the window if \code{plot = TRUE}.
#' @param line Line style for the curve, e.g. \code{list(color="red", width=2, dash="solid")}.
#' @param point Marker style for \code{r(t0)}, e.g. \code{list(size=4, color="black", symbol="circle")}.
#' @param scene 3D axes settings (see examples); by default uses \code{aspectmode="data"} for proportional axes.
#' @param bg Background colors: \code{list(paper="white", plot="white")}.
#' @param tol Numeric tolerance to detect degeneracies (\eqn{\lVert r'(t0)\rVert \approx 0} or \eqn{\lVert r' \times r''\rVert \approx 0}).
#'
#' @return A list with:
#' \itemize{
#'   \item \code{kappa}: curvature at \code{t0}.
#'   \item \code{tau}: torsion at \code{t0} (NA if \eqn{r'(t0) \times r''(t0)} is nearly zero).
#'   \item \code{t0}: the evaluation point.
#'   \item \code{point}: vector \code{c(x,y,z)} with \eqn{r(t0)}.
#'   \item \code{r1}, \code{r2}, \code{r3}: numeric derivatives \eqn{r'(t0)}, \eqn{r''(t0)}, \eqn{r'''(t0)}.
#' }
#'
#' @examples
#' # Example curve
#' X <- function(t) t^2*cos(t)
#' Y <- function(t) t^3*sin(3*t)
#' Z <- function(t) t
#' res <- curvature_torsion3d(X, Y, Z, t0 = pi)
#' res$kappa; res$tau
#'
#' # \donttest{ if (requireNamespace("plotly", quietly = TRUE)) {
#' #   curvature_torsion3d(
#' #     X, Y, Z, t0 = pi, plot = TRUE,
#' #     window = 1.0, n_samples = 200,
#' #     line  = list(color="red",  width=2),
#' #     point = list(color="black", size=5, symbol="circle")
#' #   )
#' # } }
#' @export
curvature_torsion3d <- function(
    X, Y, Z, t0,
    h = 1e-4,
    plot = FALSE,
    window = 1.0,
    n_samples = 200,
    line  = list(color = "red",   width = 2, dash = "solid"),
    point = list(color = "black", size  = 5, symbol = "circle"),
    scene = list(
      aspectmode = "data",
      xaxis = list(title = "x(t)"),
      yaxis = list(title = "y(t)"),
      zaxis = list(title = "z(t)")
    ),
    bg = list(paper = "white", plot = "white"),
    tol = 1e-10
) {
  # --------- minimal validation ----------
  if (!is.function(X) || !is.function(Y) || !is.function(Z)) {
    stop("X, Y, Z must be functions of one variable t.", call. = FALSE)
  }
  if (!is.numeric(t0) || length(t0) != 1L || !is.finite(t0)) {
    stop("t0 must be a finite numeric scalar.", call. = FALSE)
  }
  if (!is.numeric(h) || h <= 0) stop("h must be positive.", call. = FALSE)

  # --------- numeric derivatives (centered) ----------
  d1 <- function(f, t, h) (f(t + h) - f(t - h)) / (2*h)
  d2 <- function(f, t, h) (f(t + h) - 2*f(t) + f(t - h)) / (h*h)
  # Third derivative, centered finite difference (step h):
  # f'''(t) ≈ (f(t-2h) - 2f(t-h) + 2f(t+h) - f(t+2h)) / (2 h^3)
  d3 <- function(f, t, h) (f(t - 2*h) - 2*f(t - h) + 2*f(t + h) - f(t + 2*h)) / (2*h^3)

  r  <- function(t) c(X(t), Y(t), Z(t))
  r1 <- c(d1(X, t0, h), d1(Y, t0, h), d1(Z, t0, h))
  r2 <- c(d2(X, t0, h), d2(Y, t0, h), d2(Z, t0, h))
  r3 <- c(d3(X, t0, h), d3(Y, t0, h), d3(Z, t0, h))

  # --------- vector algebra ----------
  dot  <- function(a,b) sum(a*b)
  norm <- function(a) sqrt(dot(a,a))
  cross <- function(a,b) c(a[2]*b[3] - a[3]*b[2],
                           a[3]*b[1] - a[1]*b[3],
                           a[1]*b[2] - a[2]*b[1])

  r1_norm <- norm(r1)
  if (r1_norm < tol) {
    stop("||r'(t0)|| ≈ 0. Curvature/torsion undefined at t0 (possible singular point).", call. = FALSE)
  }

  c12   <- cross(r1, r2)
  n_c12 <- norm(c12)

  kappa <- n_c12 / (r1_norm^3)

  # torsion: triple product / ||r' x r''||^2
  # det([r1 r2 r3]) = r1 · (r2 × r3)
  triple <- dot(r1, cross(r2, r3))
  tau <- if (n_c12 < tol) {
    warning("||r'(t0) × r''(t0)|| ≈ 0. Torsion undefined at t0; returning NA.", call. = FALSE)
    NA_real_
  } else {
    triple / (n_c12^2)
  }

  # --------- optional plot ----------
  if (isTRUE(plot)) {
    if (!requireNamespace("plotly", quietly = TRUE)) {
      warning("To plot you need 'plotly' installed.", call. = FALSE)
    } else {
      a <- t0 - window/2
      b <- t0 + window/2
      # Reuse curve_sample3d if available; otherwise sample directly:
      ts <- seq(a, b, length.out = n_samples)
      data <- tibble::tibble(
        t = ts,
        x = vapply(ts, X, numeric(1)),
        y = vapply(ts, Y, numeric(1)),
        z = vapply(ts, Z, numeric(1))
      )
      p0 <- r(t0)

      plt <- plotly::plot_ly(
        data = data,
        x = ~x, y = ~y, z = ~z,
        type = "scatter3d", mode = "lines",
        line = line
      ) |>
        plotly::add_trace(
          x = p0[1], y = p0[2], z = p0[3],
          type = "scatter3d", mode = "markers",
          marker = point,
          showlegend = FALSE,
          hoverinfo = "none"
        ) |>
        plotly::layout(
          title = paste0("κ ≈ ", signif(kappa, 6),
                         "   ·   τ ≈ ", ifelse(is.na(tau), "NA", signif(tau, 6))),
          scene = scene,
          paper_bgcolor = bg$paper,
          plot_bgcolor  = bg$plot
        )

      print(plt)
    }
  }

  list(
    kappa  = kappa,
    tau    = tau,
    t0     = t0,
    point  = r(t0),
    r1 = r1, r2 = r2, r3 = r3
  )
}
