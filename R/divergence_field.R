#' Divergence of a vector field (with optional visualization)
#'
#' Computes the divergence of a vector field F:R^n→R^n using central differences.
#' If plot=TRUE, shows a visualization (2D quiver or 3D starburst).
#'
#' @param F Vector field: function(x) with x numeric vector, returning numeric vector of same length.
#' @param x0 Point of evaluation (numeric vector).
#' @param h Step size(s) for finite differences. Default NULL → automatic.
#' @param plot Logical. If TRUE and n=2 or 3, creates a visualization with plotly.
#'
#' @return Numeric scalar: divergence at x0.
#' @export
divergence_field <- function(F, x0, h = NULL, plot = FALSE) {
  if (!is.function(F)) stop("'F' must be a function(x).")
  if (!is.numeric(x0)) stop("'x0' must be numeric.")
  x0 <- as.numeric(x0)
  n  <- length(x0)

  Fx0 <- F(x0)
  if (!is.numeric(Fx0) || length(Fx0) != n) stop("'F(x)' must return numeric vector of length n.")

  if (is.null(h)) {
    h <- 1e-4 * (1 + abs(x0))
  } else if (length(h) == 1L) {
    h <- rep(h, n)
  } else if (length(h) != n) {
    stop("h must be scalar, vector of length n, or NULL.")
  }

  div <- 0
  for (i in seq_len(n)) {
    ei <- rep(0, n); ei[i] <- h[i]
    Fp <- F(x0 + ei)
    Fm <- F(x0 - ei)
    div <- div + (Fp[i] - Fm[i]) / (2 * h[i])
  }

  # Optional plot
  if (plot && (n %in% c(2,3))) {
    if (!requireNamespace("plotly", quietly = TRUE)) {
      warning("Install 'plotly' for visualization.")
    } else {
      if (n == 2) {
        # Small grid around x0
        xs <- seq(x0[1]-1, x0[1]+1, length.out=5)
        ys <- seq(x0[2]-1, x0[2]+1, length.out=5)
        grid <- expand.grid(x=xs, y=ys)
        vecs <- t(apply(grid, 1, F))

        fig <- plotly::plot_ly()
        fig <- fig |>
          plotly::add_trace(
            x = grid$x, y = grid$y,
            u = vecs[,1], v = vecs[,2],
            type = "cone", sizemode="absolute",
            sizeref=0.5,
            colorscale = if (div >= 0) "Reds" else "Blues",
            showscale=FALSE
          ) |>
          plotly::layout(
            title = sprintf("Divergence in 2D at x0 = (%.2f, %.2f): %.4f", x0[1], x0[2], div),
            xaxis = list(title="x"),
            yaxis = list(title="y")
          )
        print(fig)
      }
      if (n == 3) {
        # Starburst at point
        directions <- rbind(
          c(1,0,0), c(-1,0,0),
          c(0,1,0), c(0,-1,0),
          c(0,0,1), c(0,0,-1)
        )
        arrows <- sweep(directions, 2, x0, "+")

        fig <- plotly::plot_ly()
        for (i in 1:nrow(directions)) {
          fig <- fig |>
            plotly::add_trace(
              x = c(x0[1], arrows[i,1]),
              y = c(x0[2], arrows[i,2]),
              z = c(x0[3], arrows[i,3]),
              type = "scatter3d", mode = "lines",
              line = list(
                color = if (div >= 0) "red" else "blue",
                width = 6
              ),
              showlegend = FALSE
            )
        }
        fig <- fig |>
          plotly::layout(
            title = sprintf("Divergence in 3D at x0=(%.2f,%.2f,%.2f): %.4f", x0[1],x0[2],x0[3],div),
            scene = list(
              xaxis=list(title="x"),
              yaxis=list(title="y"),
              zaxis=list(title="z")
            )
          )
        print(fig)
      }
    }
  }

  div
}
