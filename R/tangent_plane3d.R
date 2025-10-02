#' Tangent plane and normal of a surface z = f(x,y) (numeric)
#'
#' @param F Function \code{F(x,y)} returning \code{z}.
#' @param point Vector \code{c(a,b)} where the tangent plane is computed.
#' @param h Step size for centered finite differences.
#' @param plot Logical; if \code{TRUE}, produces a 3D figure.
#' @param x_window,y_window Half-widths in \code{x} and \code{y} for the surface.
#' @param z_window Half-height for fixing the \code{z}-range around \code{f(a,b)}.
#' @param grid Number of divisions per axis for the surface mesh.
#' @param plane_window Half-width of the rectangle where the tangent plane \code{g} is drawn.
#' @param vec_N_factor Length of the normal vector segment (scale).
#' @param surface_opacity,plane_opacity Opacities of surface and plane (0–1).
#' @param colors List of colors/scales: \code{list(surface = "Viridis", plane = "Reds",
#'   xcut = "black", ycut = "black", point = "blue", normal = "red")}.
#'   \itemize{
#'     \item \code{surface}, \code{plane}: plotly colorscales (e.g. "Viridis", "Cividis", "Blues").
#'     \item \code{xcut}, \code{ycut}, \code{point}, \code{normal}: simple colors.
#'   }
#' @param show_surface_grid Logical; show mesh grid over the surface.
#' @param surface_grid_color Color of the surface mesh grid (e.g. \code{"rgba(0,0,0,0.35)"}).
#' @param surface_grid_width Width of the surface mesh grid.
#' @param show_axis_grid Logical; show grid on x,y,z axes.
#' @param axis_grid_color Color of the axis grid.
#' @param axis_grid_width Width of the axis grid.
#' @param scene Base 3D scene settings (axis/grid options are added).
#' @param bg Background colors: \code{list(paper="white", plot="white")}.
#'
#' @return A list with \code{fx}, \code{fy}, \code{f0}, \code{normal_unit},
#'   \code{normal_raw}, \code{plane_fun}, \code{plane_coeff}.
#' @export
tangent_plane3d <- function(
    F,
    point,
    h = 1e-4,
    plot = FALSE,
    x_window = 4, y_window = 4, z_window = 4,
    grid = 50,
    plane_window = 1,
    vec_N_factor = 1,
    surface_opacity = 0.85,
    plane_opacity   = 0.7,
    colors = list(surface = "Viridis", plane = "Reds",
                  xcut = "black", ycut = "black",
                  point = "blue", normal = "red"),
    show_surface_grid = FALSE,
    surface_grid_color = "rgba(0,0,0,0.35)",
    surface_grid_width = 1,
    show_axis_grid = TRUE,
    axis_grid_color = "#e5e5e5",
    axis_grid_width = 1,
    scene = list(aspectmode = "data",
#' @noRd
                 xaxis = list(title = "x"),
                 yaxis = list(title = "y"),
                 zaxis = list(title = "z")),
#' @noRd
#' @noRd
#' @noRd
#' @noRd
#' @noRd
#' @noRd
#' @noRd
    bg = list(paper = "white", plot = "white")
) {
  if (!is.numeric(point) || length(point) != 2L || any(!is.finite(point))) {
    stop("'point' must be numeric finite vector c(a,b).", call. = FALSE)
  }

  a <- point[1]; b <- point[2]

  # Partial derivatives (centered)
  d1x <- function(F, x, y, h) (F(x + h, y) - F(x - h, y)) / (2*h)
  d1y <- function(F, x, y, h) (F(x, y + h) - F(x, y - h)) / (2*h)

  f0 <- F(a, b)
  fx <- d1x(F, a, b, h)
  fy <- d1y(F, a, b, h)

  # Tangent plane g(x,y)
  g <- function(x, y) f0 + fx * (x - a) + fy * (y - b)

  # Normal (raw and unit)
  n_raw  <- c(-fx, -fy, 1)
  n_unit <- n_raw / sqrt(sum(n_raw^2))

  # Coefficients ax + by + cz + d = 0 (normal = n_raw, passes through (a,b,f0))
  dcoef <- fx * a + fy * b - f0
  plane_coeff <- c(a = -fx, b = -fy, c = 1, d = dcoef)

  if (isTRUE(plot)) {
    if (!requireNamespace("plotly", quietly = TRUE)) {
      warning("To plot you need 'plotly' installed.", call. = FALSE)
    } else {
      # Surface mesh
      xs <- seq(a - x_window, a + x_window, length.out = grid)
      ys <- seq(b - y_window, b + y_window, length.out = grid)
      Zf <- outer(ys, xs, function(yy, xx) F(xx, yy))

      # Plane mesh (smaller window)
      xg <- seq(a - plane_window, a + plane_window, length.out = max(20, grid %/% 2))
      yg <- seq(b - plane_window, b + plane_window, length.out = max(20, grid %/% 2))
      Zg <- outer(yg, xg, function(yy, xx) g(xx, yy))

      # Cuts
      xcut <- seq(a - plane_window, a + plane_window, length.out = 200)
      ycut <- seq(b - plane_window, b + plane_window, length.out = 200)
      z_xcut <- vapply(xcut, function(x) F(x, b), numeric(1))
      z_ycut <- vapply(ycut, function(y) F(a, y), numeric(1))

      # Point and normal
      p0 <- c(a, b, f0); p1 <- p0 + vec_N_factor * n_unit

      # Contours on surface f
      contours_arg <- NULL
      if (isTRUE(show_surface_grid)) {
        contours_arg <- list(
          x = list(show = TRUE, color = surface_grid_color, width = surface_grid_width),
          y = list(show = TRUE, color = surface_grid_color, width = surface_grid_width),
          z = list(show = FALSE)
        )
      }

      plt <- plotly::plot_ly()
      # Surface f
      plt <- plt |>
        plotly::add_surface(
          x = xs, y = ys, z = Zf,
          colorscale = colors$surface,
          showscale = FALSE, opacity = surface_opacity,
          contours = contours_arg
        )
      # Plane g
      plt <- plt |>
        plotly::add_surface(
          x = xg, y = yg, z = Zg,
          colorscale = colors$plane,
          showscale = FALSE, opacity = plane_opacity
        )
      # Cut y=b
      plt <- plt |>
        plotly::add_trace(
          x = xcut, y = rep(b, length(xcut)), z = z_xcut,
          type = "scatter3d", mode = "lines",
          line = list(color = colors$xcut, width = 2),
          showlegend = FALSE, hoverinfo = "none"
        )
      # Cut x=a
      plt <- plt |>
        plotly::add_trace(
          x = rep(a, length(ycut)), y = ycut, z = z_ycut,
          type = "scatter3d", mode = "lines",
          line = list(color = colors$ycut, width = 2),
          showlegend = FALSE, hoverinfo = "none"
        )
      # Point
      plt <- plt |>
        plotly::add_trace(
          x = p0[1], y = p0[2], z = p0[3],
          type = "scatter3d", mode = "markers",
          marker = list(color = colors$point, size = 4),
          showlegend = FALSE, hoverinfo = "none"
        )
      # Normal vector (segment)
      plt <- plt |>
        plotly::add_trace(
          x = c(p0[1], p1[1]), y = c(p0[2], p1[2]), z = c(p0[3], p1[3]),
          type = "scatter3d", mode = "lines",
          line = list(color = colors$normal, width = 3),
          showlegend = FALSE, hoverinfo = "none"
        )

      # z-range around f0
      zmin <- f0 - z_window; zmax <- f0 + z_window

      # Apply axis grid without overwriting other user settings
      scene_final <- scene
      if (is.null(scene_final$xaxis)) scene_final$xaxis <- list()
      if (is.null(scene_final$yaxis)) scene_final$yaxis <- list()
      if (is.null(scene_final$zaxis)) scene_final$zaxis <- list()
      scene_final$xaxis <- modifyList(scene_final$xaxis, list(showgrid = show_axis_grid, gridcolor = axis_grid_color, gridwidth = axis_grid_width))
      scene_final$yaxis <- modifyList(scene_final$yaxis, list(showgrid = show_axis_grid, gridcolor = axis_grid_color, gridwidth = axis_grid_width))
      scene_final$zaxis <- modifyList(scene_final$zaxis, list(showgrid = show_axis_grid, gridcolor = axis_grid_color, gridwidth = axis_grid_width, range = c(zmin, zmax)))

      plt <- plt |>
        plotly::layout(
          title = "Tangent plane and normal at (a,b)",
          scene = scene_final,
          paper_bgcolor = bg$paper,
          plot_bgcolor  = bg$plot
        )

      print(plt)
    }
  }

  list(
    fx = fx, fy = fy, f0 = f0,
    normal_unit = n_unit,
    normal_raw  = n_raw,
    plane_fun   = g,
    plane_coeff = plane_coeff
  )
}
