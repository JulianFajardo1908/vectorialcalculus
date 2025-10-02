#' Region in the xy–plane between two curves y = H1(x), y = H2(x) (2D or 3D)
#'
#' For \eqn{x \in [a,b]}, this plots the curves \eqn{y = H_1(x)} and \eqn{y = H_2(x)},
#' the vertical partitions every \eqn{D} (using \eqn{n=\lceil(b-a)/D\rceil}),
#' and optionally fills the region. With \code{as_3d=TRUE}, it draws the scene
#' in 3D over the plane \eqn{z=\texttt{plane\_z}} (useful to combine with other surfaces).
#'
#' @param H1,H2 \code{function(x)}.
#' @param a,b Endpoints of the interval in \eqn{x}.
#' @param D Desired step between vertical partitions (we use \eqn{n=\lceil (b-a)/D\rceil}).
#' @param plot \code{TRUE}/\code{FALSE}. If \code{TRUE}, draw with \pkg{plotly}.
#' @param n_curve Integer, number of samples for the curves.
#' @param fill \code{TRUE}/\code{FALSE}. If \code{TRUE}, fill the region.
#' @param fillcolor (2D) Fill color (rgba or simple color).
#' @param boundary_line Line style for the boundary curves.
#' @param partition_line Line style for the vertical partitions.
#' @param show_end_edges \code{TRUE}/\code{FALSE} to also include \eqn{x=a} and \eqn{x=b}.
#' @param axis_equal (2D) \code{TRUE}/\code{FALSE} to use the same scale in x and y.
#' @param as_3d \code{TRUE}/\code{FALSE} to draw in 3D over \code{plane_z}.
#' @param plane_z Height of the \eqn{z}-plane (default 0).
#' @param n_u (3D) “Transverse” resolution for the filled surface.
#' @param surface_colorscale (3D) Plotly colorscale for the fill (e.g., "Blues").
#' @param surface_opacity (3D) Opacity of the fill (0–1).
#' @param show_surface_grid (3D) \code{TRUE}/\code{FALSE} grid over the filled surface.
#' @param surface_grid_color,surface_grid_width (3D) Aesthetics of that grid.
#' @param scene (3D) Scene settings (default \code{aspectmode="data"}).
#' @param bg (3D) Background colors \code{list(paper="white", plot="white")}.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{x}, \code{y1}, \code{y2}, \code{y_low}, \code{y_high}
#'   \item \code{x_part} (partition positions)
#'   \item \code{fig} (plotly object if \code{plot=TRUE}, otherwise \code{NULL})
#' }
#'
#' @examples
#' H1 <- function(x) sin(x) - 0.3
#' H2 <- function(x) 0.8 + 0.2*cos(2*x)
#'
#' # Classic 2D:
#' xy_region(H1, H2, a = 0, b = 2*pi, D = 0.5,
#'           plot = TRUE, fill = TRUE)
#'
#' # 3D over z=0 (very light filled surface):
#' xy_region(H1, H2, a = 0, b = 2*pi, D = 0.5,
#'           plot = TRUE, as_3d = TRUE, fill = TRUE,
#'           surface_colorscale = "Blues", surface_opacity = 0.30,
#'           show_surface_grid = TRUE)
#'
#' @export
xy_region <- function(
    H1, H2,
    a, b,
    D,
    plot = TRUE,
    n_curve = 800,
    fill = FALSE,
    fillcolor = "rgba(49,130,189,0.25)",
    boundary_line  = list(color = "blue", width = 2),
    partition_line = list(color = "blue", width = 1),
    show_end_edges = TRUE,
    axis_equal = TRUE,
    as_3d = FALSE,
    plane_z = 0,
    n_u = 30,
    surface_colorscale = "Blues",
    surface_opacity    = 0.30,
    show_surface_grid  = TRUE,
    surface_grid_color = "rgba(60,80,200,0.25)",
    surface_grid_width = 1,
    scene = list(
      aspectmode = "data",
      xaxis = list(title = "x"),
#' @noRd
#' @noRd
#' @noRd
      yaxis = list(title = "y"),
      zaxis = list(title = "z")
    ),
    bg = list(paper = "white", plot = "white")
) {
  if (b <= a) stop("'b' must be > 'a'.", call. = FALSE)
  if (D <= 0)  stop("'D' must be > 0.", call. = FALSE)

  # Number of partitions
  n <- ceiling((b - a) / D)
  if (n < 1) n <- 1

  # Curve mesh
  xs  <- seq(a, b, length.out = n_curve)
  y1  <- H1(xs)
  y2  <- H2(xs)
  ylo <- pmin(y1, y2)
  yhi <- pmax(y1, y2)

  # Partition positions (includes a and b)
  x_part <- a + (b - a) * (0:n) / n

  fig <- NULL
  if (isTRUE(plot)) {
    if (!requireNamespace("plotly", quietly = TRUE)) {
      warning("To plot you need 'plotly' installed.", call. = FALSE)
    } else {
      if (!isTRUE(as_3d)) {
        # ---------- 2D MODE ----------
        fig <- plotly::plot_ly()

        if (isTRUE(fill)) {
          fig <- fig |>
            plotly::add_trace(
              x = c(xs, rev(xs)),
              y = c(ylo, rev(yhi)),
              type = "scatter",
              mode = "lines",
              fill = "toself",
              fillcolor = fillcolor,
              line = list(color = "rgba(0,0,0,0)"),
              hoverinfo = "none",
              showlegend = FALSE
            )
        }

        fig <- fig |>
          plotly::add_trace(
            x = xs, y = y1, type = "scatter", mode = "lines",
            line = boundary_line, hoverinfo = "none", showlegend = FALSE
          ) |>
          plotly::add_trace(
            x = xs, y = y2, type = "scatter", mode = "lines",
            line = boundary_line, hoverinfo = "none", showlegend = FALSE
          )

        k_seq <- if (isTRUE(show_end_edges)) 0:n else 1:(n-1)
        if (length(k_seq) > 0) {
          for (k in k_seq) {
            xk <- a + (b - a) * k / n
            yk1 <- H1(xk); yk2 <- H2(xk)
            yklo <- min(yk1, yk2); ykhi <- max(yk1, yk2)
            fig <- fig |>
              plotly::add_trace(
                x = c(xk, xk),
                y = c(yklo, ykhi),
                type = "scatter", mode = "lines",
                line = partition_line, hoverinfo = "none", showlegend = FALSE
              )
          }
        }

        layout_args <- list(
          xaxis = list(title = "x"),
          yaxis = list(title = "y")
        )
        if (isTRUE(axis_equal)) {
          layout_args$yaxis$scaleanchor <- "x"
          layout_args$yaxis$scaleratio  <- 1
        }
        fig <- fig |>
          plotly::layout(title = "Region in the xy-plane",
                         xaxis = layout_args$xaxis,
                         yaxis = layout_args$yaxis)

      } else {
        # ---------- 3D MODE ----------
        fig <- plotly::plot_ly()

        # Fill as a surface S(u,j) with u in [0,1]: y(u) = (1-u)*ylo + u*yhi
        if (isTRUE(fill)) {
          U <- seq(0, 1, length.out = n_u)
          Xmat <- matrix(rep(xs, each = n_u), nrow = n_u)
          # Ymat[i,j] = (1-U[i])*ylo[j] + U[i]*yhi[j]
          Ymat <- outer(U, seq_along(xs), function(u, j) (1 - u) * ylo[j] + u * yhi[j])
          Zmat <- matrix(plane_z, nrow = n_u, ncol = length(xs))

          contours_arg <- if (isTRUE(show_surface_grid)) list(
            x = list(show = TRUE, color = surface_grid_color, width = surface_grid_width),
            y = list(show = TRUE, color = surface_grid_color, width = surface_grid_width),
            z = list(show = FALSE)
          ) else NULL

          fig <- fig |>
            plotly::add_surface(
              x = Xmat, y = Ymat, z = Zmat,
              colorscale = surface_colorscale,
              showscale  = FALSE,
              opacity    = surface_opacity,
              contours   = contours_arg
            )
        }

        # Boundary curves over z=plane_z
        fig <- fig |>
          plotly::add_trace(
            x = xs, y = y1, z = rep(plane_z, length(xs)),
            type = "scatter3d", mode = "lines",
            line = boundary_line, hoverinfo = "none", showlegend = FALSE
          ) |>
          plotly::add_trace(
            x = xs, y = y2, z = rep(plane_z, length(xs)),
            type = "scatter3d", mode = "lines",
            line = boundary_line, hoverinfo = "none", showlegend = FALSE
          )

        # Vertical partitions over z=plane_z
        k_seq <- if (isTRUE(show_end_edges)) 0:n else 1:(n-1)
        if (length(k_seq) > 0) {
          for (k in k_seq) {
            xk <- a + (b - a) * k / n
            yk1 <- H1(xk); yk2 <- H2(xk)
            yklo <- min(yk1, yk2); ykhi <- max(yk1, yk2)
            fig <- fig |>
              plotly::add_trace(
                x = c(xk, xk),
                y = c(yklo, ykhi),
                z = c(plane_z, plane_z),
                type = "scatter3d", mode = "lines",
                line = partition_line, hoverinfo = "none", showlegend = FALSE
              )
          }
        }

        # 3D scene with true proportions (aspectmode="data")
        fig <- fig |>
          plotly::layout(
            title = sprintf("Region on plane z = %g", plane_z),
            scene = scene,
            paper_bgcolor = bg$paper,
            plot_bgcolor  = bg$plot
          )
      }

      print(fig)
    }
  }

  list(
    x = xs, y1 = y1, y2 = y2,
    y_low = ylo, y_high = yhi,
    x_part = x_part,
    fig = fig
  )
}
