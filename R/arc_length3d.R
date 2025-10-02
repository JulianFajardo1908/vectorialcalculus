# vector_calculus.R
# Package: vectorcalculus
# Single file with internal helpers, sampling, plotting, and arc_length3d()
# All computations are numeric (no symbolic math). Optional visualization with plotly.

# -------------------------------------------------------------------
#' @title Sampling, 3D Visualization, and Arc Length (numeric)
#' @description
#' Collection of utilities for 3D parametric curves \eqn{(X(t), Y(t), Z(t))}:
#' - \code{arc_length3d()}: arc length by numerical integration.
#' - \code{curve_sample3d()}: sample the curve on \eqn{[a,b]}.
#' - \code{plot_curve3d()}: plot the curve with plotly (if available).
#'
#' Derivatives are computed by centered finite differences and the integral
#' can be solved either with Romberg (\pkg{pracma}) or base R \code{integrate()}.
#'
#' @details
#' This file also defines the internal utility \code{.d_central()} (not exported)
#' for numeric derivatives.
#'
#' @keywords calculus vector-calculus arc-length plotly numerical
#' @importFrom tibble tibble
#' @importFrom stats integrate
#' @importFrom utils head tail
#' @name vector_calculus
NULL

# -------------------------------------------------------------------
# Internal utilities (not exported)
# Centered finite difference for numeric derivative of f at t with step h
.d_central <- function(f, t, h) (f(t + h) - f(t - h)) / (2*h)

# -------------------------------------------------------------------
#' Sample a 3D parametric curve
#'
#' Generates a tibble with columns \code{t, x, y, z} by sampling
#' \eqn{(X(t), Y(t), Z(t))} on \eqn{[a,b]}.
#'
#' @param X,Y,Z Functions of one variable \code{t}, e.g. \code{function(t) 2*cos(t)}.
#' @param a,b Parameter limits for \code{t}.
#' @param n_samples Number of sample points (density).
#' @return A tibble with columns \code{t, x, y, z}.
#' @examples
#' X <- function(t) 2*cos(t); Y <- function(t) 3*sin(t); Z <- function(t) t/5
#' curve_sample3d(X, Y, Z, 0, 2*pi, n_samples = 100)
#' @export
curve_sample3d <- function(X, Y, Z, a, b, n_samples = 400) {
  ts <- seq(a, b, length.out = n_samples)
  tibble::tibble(
    t = ts,
    x = vapply(ts, X, numeric(1)),
    y = vapply(ts, Y, numeric(1)),
    z = vapply(ts, Z, numeric(1))
  )
}

# -------------------------------------------------------------------
#' Plot a 3D curve with plotly
#'
#' Accepts a tibble with columns \code{t, x, y, z} (as produced by \code{curve_sample3d()}).
#' Requires \pkg{plotly} to be installed.
#'
#' @param data Tibble with columns \code{t, x, y, z}.
#' @param mode Plotly trace mode: \code{"lines"}, \code{"lines+markers"}, etc.
#' @param line Line style list (e.g. \code{list(color="blue", width=3, dash="solid")}).
#' @param marker Marker style list (or \code{NULL} to omit).
#' @param title Plot title.
#' @param scene List of 3D axes, e.g. \code{list(xaxis=list(title="x(t)"), ...)}.
#' @param bg Background colors: \code{list(paper="white", plot="white")}.
#' @return A plotly object (printed).
#' @examples
#' data <- curve_sample3d(function(t) 2*cos(t),
#'                        function(t) 3*sin(t),
#'                        function(t) t/5, 0, 2*pi, 100)
#' # \donttest{ if (requireNamespace("plotly", quietly = TRUE)) {
#' #   plot_curve3d(data, line = list(color="red", width=4))
#' # } }
#' @export
plot_curve3d <- function(
    data,
    mode  = "lines",
    line  = list(color = "blue", width = 3, dash = "solid"),
    marker = NULL,
    title  = NULL,
    scene  = list(
      xaxis = list(title = "x(t)"),
      yaxis = list(title = "y(t)"),
      zaxis = list(title = "z(t)")
    ),
    bg = list(paper = "white", plot = "white")
) {
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("You need to install 'plotly' to use plot_curve3d().")
  }
  plt <- plotly::plot_ly(
    data = data, x = ~x, y = ~y, z = ~z,
    type = "scatter3d", mode = mode, line = line
  )
  if (!is.null(marker)) {
    plt <- plt |>
      plotly::style(marker = marker)
  }
  plotly::layout(
    plt,
    title = title,
    scene = scene,
    paper_bgcolor = bg$paper,
    plot_bgcolor  = bg$plot
  )
}

# -------------------------------------------------------------------
#' Arc length of a 3D parametric curve (numeric)
#'
#' Computes \eqn{\int_a^b \|\mathbf{r}'(t)\|\,dt} using centered finite
#' differences for derivatives and numerical integration (Romberg or \code{integrate}).
#' Optionally, it can also plot the curve with \pkg{plotly}.
#'
#' @param X,Y,Z Functions of \code{t}.
#' @param a,b Parameter limits for \code{t}.
#' @param h Step size for centered finite differences (numeric derivatives).
#' @param method_int \code{"romberg"} (requires \pkg{pracma}) or \code{"integrate"} (base R).
#' @param n_samples Number of sample points for plotting (if \code{plot = TRUE}).
#' @param plot \code{TRUE}/\code{FALSE}. If \code{TRUE}, plots with \code{plotly}.
#' @param plot_mode,plot_line,plot_marker,plot_title,plot_scene,plot_bg
#'   Styling parameters passed to \code{plot_curve3d()}.
#'
#' @return Numeric value: arc length.
#' @examples
#' X <- function(t) t^2*cos(t)
#' Y <- function(t) t^3*sin(3*t)
#' Z <- function(t) t
#' arc_length3d(X, Y, Z, 0, 2*pi)
#' # \donttest{ if (requireNamespace("plotly", quietly = TRUE)) {
#' #   arc_length3d(X, Y, Z, 0, 2*pi, plot = TRUE,
#' #              plot_line = list(color="red", width=3),
#' #              n_samples = 300)
#' # } }
#' @export
arc_length3d <- function(
    X, Y, Z, a, b,
    h = 1e-6,
    method_int = c("romberg", "integrate"),
    n_samples = 400,
    plot = FALSE,
    plot_mode  = "lines",
    plot_line  = list(color = "blue", width = 3, dash = "solid"),
    plot_marker = NULL,
    plot_title  = NULL,
    plot_scene  = list(
      xaxis = list(title = "x(t)"),
      yaxis = list(title = "y(t)"),
      zaxis = list(title = "z(t)")
    ),
    plot_bg = list(paper = "white", plot = "white")
) {
  method_int <- match.arg(method_int)

  speed <- function(t) {
    dx <- .d_central(X, t, h)
    dy <- .d_central(Y, t, h)
    dz <- .d_central(Z, t, h)
    sqrt(dx*dx + dy*dy + dz*dz)
  }

  length_val <- switch(
    method_int,
    romberg   = {
      if (!requireNamespace("pracma", quietly = TRUE)) {
        stop("For method_int='romberg' you need to install 'pracma'.")
      }
      pracma::romberg(speed, a, b)$value
    },
    integrate = stats::integrate(speed, lower = a, upper = b, rel.tol = 1e-8)$value
  )

  if (isTRUE(plot)) {
    data <- curve_sample3d(X, Y, Z, a, b, n_samples = n_samples)
    plot_curve3d(
      data,
      mode  = plot_mode,
      line  = plot_line,
      marker = plot_marker,
      title = if (is.null(plot_title)) paste0("Arc length ∈ ", signif(length_val, 6)) else plot_title,
      scene  = plot_scene,
      bg     = plot_bg
    ) |> print()
  }

  length_val
}
