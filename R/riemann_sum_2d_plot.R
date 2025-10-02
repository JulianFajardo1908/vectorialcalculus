#' 2D Riemann sums (upper, lower, midpoint) with a 3D plot
#'
#' Visualize 2D Riemann sums for a scalar field \eqn{f(x,y)} over a rectangular
#' domain \eqn{[x_{\min},x_{\max}]\times[y_{\min},y_{\max}]}. The function
#' computes **upper**, **lower**, and **midpoint** sums and renders a 3D figure
#' showing step tiles for the selected methods. Optionally overlays the true
#' surface \eqn{z=f(x,y)} and a base grid on the \eqn{xy}-plane.
#'
#' Upper/lower tiles use corner samples (min/max of the four corners), and the
#' midpoint tiles sample at the cell center.
#'
#' @param f \code{function(x,y)} returning a numeric scalar \eqn{f(x,y)}.
#' @param xlim,ylim Numeric length-2 vectors \code{c(min,max)} for the domain.
#' @param nx,ny Positive integers: number of subintervals along \eqn{x} and \eqn{y}.
#' @param methods Character vector with any of \code{c("lower","upper","mid")}.
#'        Controls which Riemann tiles are drawn. (All estimates are returned.)
#' @param show_surface Logical; if \code{TRUE} overlays the true surface \eqn{z=f(x,y)}.
#' @param surface_res Integer vector \code{c(nx_s, ny_s)} for the surface mesh.
#' @param surface_colorscale Plotly colorscale for the true surface (e.g. \code{"Viridis"}).
#' @param surface_opacity Opacity of the true surface (0–1).
#' @param base_plane Logical; if \code{TRUE} draws a faint base plane at \code{z0}.
#' @param z0 Numeric; height of the base plane (typically 0).
#' @param base_opacity Opacity of the base plane (0–1).
#' @param base_color Color of the base plane.
#' @param show_base_grid Logical; if \code{TRUE} draws partition grid lines on the base plane.
#' @param grid_color Color of the base grid lines.
#' @param grid_width Line width of the base grid.
#' @param tile_opacity Opacity of Riemann tiles (0–1).
#' @param colors Named list for tile colors (hex or rgba), e.g.:
#'   \code{list(lower="#a1d99b", upper="#fc9272", mid="#9ecae1")}.
#' @param edge_color Edge color for vertical edges of tiles.
#' @param edge_width Line width for tile edges.
#' @param warn_heavy Logical; if \code{TRUE}, warns when \code{nx*ny} is large (many traces).
#'
#' @return A list with:
#' \itemize{
#'   \item \code{lower_sum}, \code{upper_sum}, \code{mid_sum}: numeric estimates,
#'   \item \code{dx}, \code{dy}: partition widths,
#'   \item \code{breaks}: list with \code{x_breaks}, \code{y_breaks},
#'   \item \code{figure}: the \pkg{plotly} object (or \code{NULL} if \pkg{plotly} is missing).
#' }
#'
#' @examples
#' \dontshow{if (interactive()) \{}
#' # A smooth bump:
#' f <- function(x,y) exp(-(x^2 + y^2)) * (1 + 0.3 * cos(3*x) * sin(2*y))
#' out <- riemann_sum_2d_plot(
#'   f, xlim = c(-2,2), ylim = c(-2,2),
#'   nx = 8, ny = 7, methods = c("lower","mid","upper"),
#'   show_surface = TRUE, surface_res = c(80,80),
#'   surface_colorscale = "YlGnBu", surface_opacity = 0.45
#' )
#' out$lower_sum; out$mid_sum; out$upper_sum
#' \dontshow{\}}
#'
#' @export
riemann_sum_2d_plot <- function(
    f,
    xlim, ylim,
    nx = 8, ny = 8,
    methods = c("lower","upper","mid"),
    show_surface = TRUE,
    surface_res = c(60L, 60L),
    surface_colorscale = "Viridis",
    surface_opacity = 0.5,
    base_plane = TRUE,
    z0 = 0,
    base_opacity = 0.15,
    base_color = "lightgray",
    show_base_grid = TRUE,
    grid_color = "gray50",
    grid_width = 1,
    tile_opacity = 0.92,
    colors = list(lower = "#a1d99b", upper = "#fc9272", mid = "#9ecae1"),
    edge_color = "black",
    edge_width = 1.2,
    warn_heavy = TRUE
){
  # ---- checks
  stopifnot(is.function(f))
  if (length(xlim)!=2 || xlim[2] <= xlim[1] || length(ylim)!=2 || ylim[2] <= ylim[1])
    stop("'xlim' and 'ylim' must be c(min,max) with min < max.", call. = FALSE)
  if (nx < 1 || ny < 1) stop("'nx' and 'ny' must be >= 1.", call. = FALSE)

  methods <- intersect(tolower(methods), c("lower","upper","mid"))
  if (length(methods) == 0) stop("Choose at least one of methods = c('lower','upper','mid').", call. = FALSE)

  have_plotly <- requireNamespace("plotly", quietly = TRUE)
  if (!have_plotly) warning("Package 'plotly' not found; returning estimates only.", call. = FALSE)

  # ---- partitions
  x_breaks <- seq(xlim[1], xlim[2], length.out = nx + 1)
  y_breaks <- seq(ylim[1], ylim[2], length.out = ny + 1)
  dx <- diff(x_breaks)[1]
  dy <- diff(y_breaks)[1]

  # helper: cell corners & center
  corners_and_center <- function(i, j) {
    x0 <- x_breaks[i]; x1 <- x_breaks[i+1]
    y0 <- y_breaks[j]; y1 <- y_breaks[j+1]
    xm <- 0.5 * (x0 + x1); ym <- 0.5 * (y0 + y1)
    z_corners <- c(f(x0,y0), f(x0,y1), f(x1,y0), f(x1,y1))
    list(x0=x0, x1=x1, y0=y0, y1=y1, xm=xm, ym=ym, z_corners=z_corners, z_mid=f(xm,ym))
  }

  # compute per-cell heights for each method
  z_lower <- matrix(NA_real_, nrow = nx, ncol = ny)
  z_upper <- matrix(NA_real_, nrow = nx, ncol = ny)
  z_mid   <- matrix(NA_real_, nrow = nx, ncol = ny)

  for (i in 1:nx) for (j in 1:ny) {
    cc <- corners_and_center(i,j)
    z_lower[i,j] <- min(cc$z_corners)
    z_upper[i,j] <- max(cc$z_corners)
    z_mid[i,j]   <- cc$z_mid
  }

  lower_sum <- sum(z_lower) * dx * dy
  upper_sum <- sum(z_upper) * dx * dy
  mid_sum   <- sum(z_mid)   * dx * dy

  # ---- build plot
  fig <- NULL
  if (have_plotly) {
    p <- plotly::plot_ly()

    # base plane (optional)
    if (isTRUE(base_plane)) {
      xb <- x_breaks[c(1, nx+1)]; yb <- y_breaks[c(1, ny+1)]
      p <- plotly::add_surface(
        p,
        x = xb, y = yb, z = matrix(z0, nrow = 2, ncol = 2),
        opacity = base_opacity,
        showscale = FALSE,
        colorscale = list(list(0, base_color), list(1, base_color))
      )
    }

    # base grid (optional)
    if (isTRUE(show_base_grid)) {
      # lines parallel to y (fixed x)
      for (xi in x_breaks) {
        p <- plotly::add_trace(
          p, x = c(xi, xi), y = c(ylim[1], ylim[2]), z = c(z0, z0),
          type = "scatter3d", mode = "lines",
          line = list(color = grid_color, width = grid_width),
          hoverinfo = "none", showlegend = FALSE
        )
      }
      # lines parallel to x (fixed y)
      for (yi in y_breaks) {
        p <- plotly::add_trace(
          p, x = c(xlim[1], xlim[2]), y = c(yi, yi), z = c(z0, z0),
          type = "scatter3d", mode = "lines",
          line = list(color = grid_color, width = grid_width),
          hoverinfo = "none", showlegend = FALSE
        )
      }
    }

    # helper to add a single tile (rectangle plateau) and its vertical edges
    add_tile <- function(p, x0,x1,y0,y1,zval, col, name=NULL) {
      # top face as a tiny 2x2 surface (constant z)
      p <- plotly::add_surface(
        p,
        x = c(x0, x1), y = c(y0, y1), z = matrix(zval, nrow=2, ncol=2),
        showscale = FALSE,
        opacity   = tile_opacity,
        colorscale = list(list(0, col), list(1, col)),
        name = name, inherit = FALSE
      )
      # vertical edges to base plane
      XX <- c(x0, x0, NA, x0, x1, NA, x1, x1, NA, x1, x0, NA)
      YY <- c(y0, y0, NA, y0, y1, NA, y1, y1, NA, y1, y0, NA)
      ZZ <- rep(c(z0, zval, NA), 4)
      p <- plotly::add_trace(
        p, x = XX, y = YY, z = ZZ,
        type = "scatter3d", mode = "lines",
        line = list(color = edge_color, width = edge_width),
        hoverinfo = "none", showlegend = FALSE
      )
      p
    }

    n_tiles <- nx * ny
    if (warn_heavy && n_tiles > 1600) {
      warning(sprintf("You are about to draw %d tiles (nx*ny). Consider reducing nx,ny for speed.", n_tiles),
              call. = FALSE)
    }

    # draw tiles for requested methods (layer order: lower, mid, upper)
    draw_method <- function(method) {
      col <- colors[[method]]
      nm  <- switch(method, lower="lower sum", upper="upper sum", mid="midpoint sum")
      for (i in 1:nx) for (j in 1:ny) {
        x0 <- x_breaks[i]; x1 <- x_breaks[i+1]
        y0 <- y_breaks[j]; y1 <- y_breaks[j+1]
        zval <- switch(method,
                       lower = z_lower[i,j],
                       upper = z_upper[i,j],
                       mid   = z_mid[i,j]
        )
        p <<- add_tile(p, x0,x1,y0,y1,zval, col, name = nm)
      }
    }
    for (m in c("lower","mid","upper")) if (m %in% methods) draw_method(m)

    # optional true surface
    if (isTRUE(show_surface)) {
      xs <- seq(xlim[1], xlim[2], length.out = surface_res[1])
      ys <- seq(ylim[1], ylim[2], length.out = surface_res[2])
      Z  <- outer(ys, xs, function(yy, xx) f(xx, yy))  # z[y,x]
      p <- plotly::add_surface(
        p, x = xs, y = ys, z = Z,
        colorscale = surface_colorscale,
        opacity = surface_opacity,
        showscale = FALSE,
        name = "f(x,y)"
      )
    }

    p <- plotly::layout(
      p,
      title = "2D Riemann sums (upper, lower, midpoint)",
      scene = list(
        aspectmode = "data",
        xaxis = list(title = "x"),
        yaxis = list(title = "y"),
        zaxis = list(title = "z")
      )
    )
    fig <- p
    print(fig)
  }

  list(
    lower_sum = lower_sum,
    upper_sum = upper_sum,
    mid_sum   = mid_sum,
    dx = dx, dy = dy,
    breaks = list(x_breaks = x_breaks, y_breaks = y_breaks),
    figure = fig
  )
}
