tangent_plane3d <- function(
    f,
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
                 xaxis = list(title = "x"),
                 yaxis = list(title = "y"),
                 zaxis = list(title = "z")),
    bg = list(paper = "white", plot = "white")
) {
  # Validación básica
  if (!is.numeric(point) || length(point) != 2L || any(!is.finite(point))) {
    stop("'point' must be a finite numeric vector c(x0, y0).", call. = FALSE)
  }
  f <- match.fun(f); force(f); force(h)
  a <- point[1]; b <- point[2]

  # Capturas para evitar defaults autoreferenciales
  ff <- f; hh <- h

  # Derivadas parciales (diferencias centradas)
  d1x <- function(x, y) (ff(x + hh, y) - ff(x - hh, y)) / (2 * hh)
  d1y <- function(x, y) (ff(x, y + hh) - ff(x, y - hh)) / (2 * hh)

  f0 <- ff(a, b)
  fx <- d1x(a, b)
  fy <- d1y(a, b)

  # Plano tangente g(x,y) y normales
  g <- function(x, y) f0 + fx * (x - a) + fy * (y - b)
  n_raw  <- c(-fx, -fy, 1)
  n_unit <- n_raw / sqrt(sum(n_raw^2))

  # Coeficientes ax + by + cz + d = 0 (normal = n_raw, pasa por (a,b,f0))
  dcoef <- fx * a + fy * b - f0
  plane_coeff <- c(a = -fx, b = -fy, c = 1, d = dcoef)

  fig <- NULL
  if (isTRUE(plot)) {
    if (!requireNamespace("plotly", quietly = TRUE)) {
      warning("To plot you need 'plotly' installed.", call. = FALSE)
    } else {
      # Mallas para superficie y plano
      xs <- seq(a - x_window, a + x_window, length.out = grid)
      ys <- seq(b - y_window, b + y_window, length.out = grid)
      Zf <- outer(ys, xs, function(yy, xx) ff(xx, yy))

      xg <- seq(a - plane_window, a + plane_window, length.out = max(20, grid %/% 2))
      yg <- seq(b - plane_window, b + plane_window, length.out = max(20, grid %/% 2))
      Zg <- outer(yg, xg, function(yy, xx) g(xx, yy))

      # Cortes por x=a y y=b
      xcut <- seq(a - plane_window, a + plane_window, length.out = 200)
      ycut <- seq(b - plane_window, b + plane_window, length.out = 200)
      z_xcut <- vapply(xcut, function(x) ff(x, b), numeric(1))
      z_ycut <- vapply(ycut, function(y) ff(a, y), numeric(1))

      # Punto y vector normal
      p0 <- c(a, b, f0); p1 <- p0 + vec_N_factor * n_unit

      contours_arg <- NULL
      if (isTRUE(show_surface_grid)) {
        contours_arg <- list(
          x = list(show = TRUE, color = surface_grid_color, width = surface_grid_width),
          y = list(show = TRUE, color = surface_grid_color, width = surface_grid_width),
          z = list(show = FALSE)
        )
      }

      plt <- plotly::plot_ly()
      plt <- plt |>
        plotly::add_surface(
          x = xs, y = ys, z = Zf,
          colorscale = colors$surface,
          showscale = FALSE, opacity = surface_opacity,
          contours = contours_arg
        ) |>
        plotly::add_surface(
          x = xg, y = yg, z = Zg,
          colorscale = colors$plane,
          showscale = FALSE, opacity = plane_opacity
        ) |>
        plotly::add_trace(
          x = xcut, y = rep(b, length(xcut)), z = z_xcut,
          type = "scatter3d", mode = "lines",
          line = list(color = colors$xcut, width = 2),
          showlegend = FALSE, hoverinfo = "none"
        ) |>
        plotly::add_trace(
          x = rep(a, length(ycut)), y = ycut, z = z_ycut,
          type = "scatter3d", mode = "lines",
          line = list(color = colors$ycut, width = 2),
          showlegend = FALSE, hoverinfo = "none"
        ) |>
        plotly::add_trace( # punto
          x = p0[1], y = p0[2], z = p0[3],
          type = "scatter3d", mode = "markers",
          marker = list(color = colors$point, size = 4),
          showlegend = FALSE, hoverinfo = "none"
        ) |>
        plotly::add_trace( # normal
          x = c(p0[1], p1[1]), y = c(p0[2], p1[2]), z = c(p0[3], p1[3]),
          type = "scatter3d", mode = "lines",
          line = list(color = colors$normal, width = 3),
          showlegend = FALSE, hoverinfo = "none"
        )

      # Rango z alrededor de f0
      zmin <- f0 - z_window; zmax <- f0 + z_window

      # Ajustes de escena sin pisar otros
      scene_final <- scene
      if (is.null(scene_final$xaxis)) scene_final$xaxis <- list()
      if (is.null(scene_final$yaxis)) scene_final$yaxis <- list()
      if (is.null(scene_final$zaxis)) scene_final$zaxis <- list()
      scene_final$xaxis <- utils::modifyList(scene_final$xaxis, list(showgrid = show_axis_grid, gridcolor = axis_grid_color, gridwidth = axis_grid_width))
      scene_final$yaxis <- utils::modifyList(scene_final$yaxis, list(showgrid = show_axis_grid, gridcolor = axis_grid_color, gridwidth = axis_grid_width))
      scene_final$zaxis <- utils::modifyList(scene_final$zaxis, list(showgrid = show_axis_grid, gridcolor = axis_grid_color, gridwidth = axis_grid_width, range = c(zmin, zmax)))

      plt <- plt |>
        plotly::layout(
          title = "Tangent plane and normal at (x0, y0)",
          scene = scene_final,
          paper_bgcolor = bg$paper,
          plot_bgcolor  = bg$plot
        )

      fig <- plt
      print(plt)
    }
  }

  list(
    fx = fx, fy = fy, f0 = f0,
    normal_unit = n_unit,
    normal_raw  = n_raw,
    plane_fun   = g,
    plane_coeff = plane_coeff,
    fig = fig
  )
}

