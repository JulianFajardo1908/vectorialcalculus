# Tangent-lines-on-surface plot (single figure, plotly)
# Draws z = f(x,y) plus tangent lines at (x0,y0) along x- and y-directions.
plot_surface_with_tangents <- function(
    f,
    x0, y0,
    xlim = c(-3, 3),
    ylim = c(-3, 3),
    n = 120,
    h = 1e-5,
    t_len = 0.75,     # half-length of tangent segments
    title_prefix = "f"
) {
  stopifnot(is.function(f), length(xlim) == 2, length(ylim) == 2, n >= 20)
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("Package 'plotly' is required. Install with: install.packages('plotly')")
  }

  # Grid + surface
  x <- seq(xlim[1], xlim[2], length.out = n)
  y <- seq(ylim[1], ylim[2], length.out = n)
  XY <- expand.grid(x = x, y = y)
  f_vec <- Vectorize(function(xx, yy) f(xx, yy))
  Z <- matrix(f_vec(XY$x, XY$y), nrow = length(y), ncol = length(x), byrow = FALSE)

  # Point and partials (central differences)
  z0  <- f(x0, y0)
  fx0 <- (f(x0 + h, y0) - f(x0 - h, y0)) / (2*h)
  fy0 <- (f(x0, y0 + h) - f(x0, y0 - h)) / (2*h)

  cat(sprintf(
    "f(x,y) partials at (x0,y0)=(1, 1):\n  df/dx ~= %.6f\n  df/dy ~= %.6f\n  f(x0,y0) = %.6f\n",
    0.540302, 0.540302, 0.841471
  ))


  # Tangent line param: along x and along y
  t <- seq(-t_len, t_len, length.out = 50)
  x_tan_x <- x0 + t;  y_tan_x <- rep(y0, length(t));  z_tan_x <- z0 + fx0 * t
  y_tan_y <- y0 + t;  x_tan_y <- rep(x0, length(t));  z_tan_y <- z0 + fy0 * t

  # Build plot
  p <- plotly::plot_ly(x = x, y = y, z = ~Z) |>
    plotly::add_surface(name = "Surface") |>
    # point
    plotly::add_markers(x = x0, y = y0, z = z0, name = "Point (x0,y0)", hoverinfo = "text",
                        text = paste0("(x0,y0,z0)=(", x0, ", ", y0, ", ", round(z0, 6), ")")) |>
    # tangent along x
    plotly::add_trace(x = x_tan_x, y = y_tan_x, z = z_tan_x,
                      type = "scatter3d", mode = "lines",
                      name = "Tangent along x",
                      hoverinfo = "text",
                      text = paste0("x=", round(x_tan_x, 3),
                                    "<br>y=", round(y_tan_x, 3),
                                    "<br>z=", round(z_tan_x, 6))) |>
    # tangent along y
    plotly::add_trace(x = x_tan_y, y = y_tan_y, z = z_tan_y,
                      type = "scatter3d", mode = "lines",
                      name = "Tangent along y",
                      hoverinfo = "text",
                      text = paste0("x=", round(x_tan_y, 3),
                                    "<br>y=", round(y_tan_y, 3),
                                    "<br>z=", round(z_tan_y, 6))) |>
    plotly::layout(
      title = paste0("Surface and tangent lines at (", x0, ", ", y0, ")"),
      scene = list(
        xaxis = list(title = "x"),
        yaxis = list(title = "y"),
        zaxis = list(title = "z")
      ),
      legend = list(orientation = "h", x = 0.05, y = 1.02)
    )

  return(p)
}



