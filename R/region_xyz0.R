#' Region on the plane z = z0 between y = H1(x) and y = H2(x)
#'
#' Translation of WxMaxima block \code{region_xyz0(...)}.
#' It draws the region \eqn{\{(x,y,z_0): x\in[a,b],\, \min(H_1(x),H_2(x)) \le y \le \max(H_1(x),H_2(x))\}}
#' and places vertical partitions every \eqn{D} (using \eqn{n=\lceil(b-a)/D\rceil}).
#' Internally it delegates to \code{region_xy()} with \code{as_3d=TRUE}.
#'
#' @param H1,H2 \code{function(x)} defining the lower/upper boundaries in \code{y}.
#' @param a,b Endpoints in \code{x} (with \code{b > a}).
#' @param z0 Height of the plane \eqn{z=z_0}.
#' @param D Target spacing between partitions (\eqn{n=\lceil (b-a)/D\rceil}).
#' @param plot \code{TRUE}/\code{FALSE}. If \code{TRUE}, draw with \pkg{plotly}.
#' @param n_curve Number of samples for \code{H1}, \code{H2}.
#' @param show_surface \code{TRUE}/\code{FALSE}. If \code{TRUE}, fills the region as a faint surface on \eqn{z=z_0}.
#' @param surface_colorscale Plotly colorscale for the fill (if \code{show_surface=TRUE}).
#' @param surface_opacity Fill opacity (0–1).
#' @param show_surface_grid \code{TRUE}/\code{FALSE} grid over the filled surface.
#' @param surface_grid_color,surface_grid_width Grid style for the filled surface.
#' @param boundary_line Line style for boundary curves (\code{list(color, width, dash)}).
#' @param partition_line Line style for the vertical partitions.
#' @param show_end_edges \code{TRUE}/\code{FALSE} to also include \eqn{x=a} and \eqn{x=b}.
#' @param scene 3D scene settings (default \code{aspectmode="data"}).
#' @param bg Background colors: \code{list(paper="white", plot="white")}.
#'
#' @return The same list returned by \code{region_xy()}:
#' \itemize{
#'   \item \code{x}, \code{y1}, \code{y2}, \code{y_low}, \code{y_high}
#'   \item \code{x_part} (partition positions)
#'   \item \code{fig} (plotly object if \code{plot=TRUE}, otherwise \code{NULL})
#' }
#'
#' @examples
#' H1 <- function(x) sin(x) - 0.3
#' H2 <- function(x) 0.8 + 0.2*cos(2*x)
#' region_xyz0(H1, H2, a = 0, b = 2*pi, z0 = 1, D = 0.5,
#'             plot = TRUE, show_surface = TRUE,
#'             surface_colorscale = "Blues", surface_opacity = 0.30)
#' @export
region_xyz0 <- function(
    H1, H2,
    a, b,
    z0,
    D,
    plot = TRUE,
    n_curve = 800,
    show_surface = FALSE,
    surface_colorscale = "Blues",
    surface_opacity    = 0.30,
    show_surface_grid  = TRUE,
    surface_grid_color = "rgba(60,80,200,0.25)",
    surface_grid_width = 1,
    boundary_line  = list(color = "blue", width = 2),
    partition_line = list(color = "blue", width = 1),
    show_end_edges = TRUE,
    scene = list(
      aspectmode = "data",
      xaxis = list(title = "x"),
      yaxis = list(title = "y"),
#' @noRd
#' @noRd
#' @noRd
      zaxis = list(title = "z")
    ),
    bg = list(paper = "white", plot = "white")
) {
  # Minimal validations (region_xy performs the underlying checks too)
  if (b <= a) stop("'b' must be > 'a'.", call. = FALSE)
  if (D <= 0)  stop("'D' must be > 0.",  call. = FALSE)

  # Delegate to region_xy() in 3D mode over plane z = z0
  xy_region(
    H1 = H1, H2 = H2,
    a = a, b = b, D = D,
    plot = plot,
    n_curve = n_curve,
    fill = isTRUE(show_surface),          # fill as a surface
    fillcolor = "rgba(49,130,189,0.25)",  # unused in 3D but kept for compatibility
    boundary_line  = boundary_line,
    partition_line = partition_line,
    show_end_edges = show_end_edges,
    axis_equal = FALSE,       # not applicable in 3D
    as_3d = TRUE,
    plane_z = z0,
    n_u = 30,
    surface_colorscale = surface_colorscale,
    surface_opacity    = surface_opacity,
    show_surface_grid  = show_surface_grid,
    surface_grid_color = surface_grid_color,
    surface_grid_width = surface_grid_width,
    scene = scene,
    bg = bg
  )
}
