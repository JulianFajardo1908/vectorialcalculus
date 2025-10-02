#' Gradient of a scalar field (with optional plot)
#'
#' Computes the gradient of a scalar field f:R^n→R using central differences.
#' If plot=TRUE, it will also plot the gradient vector at the evaluation point
#' (supports 2D and 3D only).
#'
#' @param f Scalar field: function(x) with x numeric vector of length n.
#' @param x0 Point where the gradient is evaluated (numeric vector).
#' @param h Step size for finite differences. Default NULL → automatic.
#' @param plot Logical. If TRUE and n=2 or 3, plot the gradient with plotly.
#'
#' @return Numeric vector with gradient components.
#' @export
gradient_scalar <- function(f, x0, h = NULL, plot = FALSE) {
  if (!is.function(f)) stop("'f' must be a function(x).")
  if (!is.numeric(x0)) stop("'x0' must be numeric.")
  x0 <- as.numeric(x0)
  n  <- length(x0)

  if (is.null(h)) {
    h <- 1e-4 * (1 + abs(x0))
  } else if (length(h) == 1L) {
    h <- rep(h, n)
  } else if (length(h) != n) {
    stop("h must be scalar, vector of length n, or NULL.")
  }

  grad <- numeric(n)
  for (i in seq_len(n)) {
    ei <- rep(0, n); ei[i] <- h[i]
    grad[i] <- (f(x0 + ei) - f(x0 - ei)) / (2 * h[i])
  }

  # Optional plot
  if (plot && (n %in% c(2,3))) {
    if (!requireNamespace("plotly", quietly = TRUE)) {
      warning("Install 'plotly' for plotting.")
    } else {
      if (n == 2) {
        p <- plotly::plot_ly() |>
          plotly::add_trace(
            x = x0[1], y = x0[2],
            type = "scatter", mode = "markers",
            marker = list(color = "blue", size = 6),
            name = "Point"
          ) |>
          plotly::add_trace(
            x = c(x0[1], x0[1] + grad[1]),
            y = c(x0[2], x0[2] + grad[2]),
            type = "scatter", mode = "lines+markers",
            line = list(color = "red", width = 3),
            marker = list(size = 4, color = "red"),
            name = "Gradient"
          ) |>
          plotly::layout(
            title = "Gradient in 2D",
            xaxis = list(title="x"), yaxis = list(title="y"),
            showlegend = TRUE
          )
        print(p)
      }
      if (n == 3) {
        p <- plotly::plot_ly() |>
          plotly::add_markers(
            x = x0[1], y = x0[2], z = x0[3],
            marker = list(color = "blue", size = 4),
            name = "Point"
          ) |>
          plotly::add_trace(
            x = c(x0[1], x0[1] + grad[1]),
            y = c(x0[2], x0[2] + grad[2]),
            z = c(x0[3], x0[3] + grad[3]),
            type = "scatter3d", mode = "lines+markers",
            line = list(color = "red", width = 6),
            marker = list(size = 3, color = "red"),
            name = "Gradient"
          ) |>
          plotly::layout(
            title = "Gradient in 3D",
            scene = list(
              xaxis = list(title="x"),
              yaxis = list(title="y"),
              zaxis = list(title="z")
            )
          )
        print(p)
      }
    }
  }

  grad
}
