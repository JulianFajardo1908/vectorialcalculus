#' Planar region between two curves y = H1(x) and y = H2(x)
#'
#' @description
#' Builds the region \eqn{\{(x,y): a \le x \le b,\ H_1(x) \le y \le H_2(x)\}} and (optionally) plots it.
#'
#' @param H1,H2 functions \code{H_i(x)}; lower/upper y-boundaries.
#' @param a,b numeric; x-interval endpoints with \eqn{a \le b}.
#' @param D integer > 0; number of x-subdivisions (default \code{200}).
#' @param plot logical; if \code{TRUE}, draw the region with \pkg{plotly}.
#' @param show_surface_grid logical; show grid/contours on the plotted surface.
#' @param surface_grid_color character; color for surface grid lines.
#' @param surface_grid_width numeric; width for surface grid lines.
#' @param scene list; Plotly 3D scene options.
#' @param bg list; background colors (\code{paper}, \code{plot}).
#' @param n_curve Integer; puntos para muestrear cada curva.
#' @param fill Logical; rellenar la región (2D) si \code{TRUE}.
#' @param fillcolor Character; color de relleno (2D).
#' @param boundary_line List; estilo de la frontera \code{y = H1, H2}.
#' @param partition_line List; estilo de las líneas de partición.
#' @param show_end_edges Logical; mostrar bordes en los extremos \code{x=a,b}.
#' @param axis_equal Logical; usar aspecto igual en 2D.
#' @param as_3d Logical; dibujar en 3D (superficie) si \code{TRUE}.
#' @param plane_z Numeric; altura \code{z} de la “lámina” si \code{as_3d = TRUE}.
#' @param n_u Integer; subdivisiones internas para discretización.
#' @param surface_colorscale Character; colorscale (3D).
#' @param surface_opacity Numeric in \eqn{[0,1]}; opacidad (3D).
#'
#' @return A list with the grid data (and a Plotly object if \code{plot = TRUE}).
#'
#' @examples
#' H1 <- function(x) 0; H2 <- function(x) 1 - x
#' xy_region(H1, H2, a = 0, b = 1, D = 200, plot = FALSE)
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
