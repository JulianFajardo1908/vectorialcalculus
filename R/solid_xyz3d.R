#' Solid (x,y,z)
#'
#' Builds and (optionally) renders with \pkg{plotly} the solid
#' \eqn{\{(x,y,z): x \in [a,b], y \in [H_1(x),H_2(x)], z \in [G_1(x,y),G_2(x,y)]\}}
#' using a curvilinear-prism parametrization. It supports:
#' \itemize{
#'   \item Display modes: \code{mode = "faces" | "wireframe" | "both"}.
#'   \item Numerical volume: \code{compute_volume = TRUE} with \code{vol_method = "adaptive"|"grid"}.
#'   \item Internal slices on planes \eqn{x=x_0}, \eqn{y=y_0}, \eqn{z=z_0}.
#'   \item Flexible colors: flat colors (R names, \code{"#RRGGBB"}, \code{"rgba(...)"});
#'         color vectors (gradients); or Plotly color scales (e.g. \code{"Viridis"}).
#'         You can pass 1 or 6 scales for the six faces.
#' }
#'
#' @param H1,H2 \code{function(x)} giving the lower/upper bounds in \eqn{y}.
#' @param G1,G2 \code{function(x,y)} giving the lower/upper bounds in \eqn{z}.
#' @param a,b Interval endpoints in \eqn{x} with \code{b > a}.
#' @param plot Logical. If \code{TRUE}, draw with \pkg{plotly}.
#' @param n_x,n_u,n_v Mesh resolution for the six faces (in \eqn{x}, \eqn{u}, \eqn{v}).
#' @param mode \code{"faces"}, \code{"wireframe"} or \code{"both"}.
#' @param show_faces Logical vector of length 6 indicating which faces to show,
#'   in the order \code{c("x=a","x=b","y=H1","y=H2","z=G1","z=G2")}. Length 1 is also accepted.
#' @param colorscales Color scales for faces. Accepts:
#'   \itemize{
#'     \item A single Plotly scale name (e.g. \code{"Blues"}, \code{"Viridis"}) applied to all faces.
#'     \item A single flat color (R name, \code{"#RRGGBB"}, \code{"rgba(...)"}) -> flat scale.
#'     \item A color vector \code{c("white","#2a9d8f",...)} -> evenly spaced gradient for all faces.
#'     \item A vector/list of 6 scales (one per face) in any of the formats above.
#'   }
#' @param opacities Opacities for faces (length 1 or 6).
#' @param show_surface_grid Logical. Draws a grid over surfaces.
#' @param surface_grid_color,surface_grid_width Grid style for surfaces.
#' @param show_edges Logical. Draw edges of each face.
#' @param edge_line Edge style (\code{list(color, width, dash)}).
#' @param wire_step Integer >= 1: draw every \code{wire_step}-th mesh line in wireframe mode.
#' @param wire_line Wireframe line style (\code{list(color, width, dash)}).
#' @param scene 3D scene options (default \code{aspectmode="data"}).
#' @param bg Background colors: \code{list(paper="white", plot="white")}.
#' @param compute_volume Logical. If \code{TRUE}, approximate volume.
#' @param vol_method Volume integration method: \code{"adaptive"} (nested \code{stats::integrate})
#'   or \code{"grid"} (trapezoidal rule on a regular grid).
#' @param nx_vol,ny_vol Grid sizes for \code{vol_method="grid"}.
#' @param slice List of slices: \code{slice = list(x = NULL, y = NULL, z = NULL)}. Each entry
#'   can be a number or numeric vector.
#' @param slice_mode Slice rendering mode: \code{"surface"}, \code{"wireframe"} or \code{"both"}.
#' @param slice_nx,slice_nu,slice_nv Mesh resolution for slices.
#' @param slice_colorscales Color scales for slices: \code{list(x=, y=, z=)} in the same formats as \code{colorscales}.
#' @param slice_opacity Opacity for slices (0-1).
#' @param slice_show_grid,slice_grid_color,slice_grid_width Grid options for slices.
#' @param slice_wire_step,slice_wire_line Wireframe step and style for slices.
#'
#' @return
#' A list with:
#' \itemize{
#'   \item \code{x_seq}, \code{u_seq}, \code{v_seq}: the parameter sequences used,
#'   \item \code{fig}: \pkg{plotly} object if \code{plot=TRUE}, otherwise \code{NULL},
#'   \item \code{volume}: \code{NULL} or a list with \code{estimate} and metadata if \code{compute_volume=TRUE}.
#' }
#'
#' @examples
#' # Note: examples avoid plotting for CRAN checks
#' H1 <- function(x) -1 - x
#' H2 <- function(x)  1 - x^2
#' G1 <- function(x,y) y
#' G2 <- function(x,y) y + 1
#' s <- solid_xyz3d(H1,H2,G1,G2, a=-1, b=1, plot=FALSE, compute_volume=TRUE, vol_method="grid",
#'                  nx_vol=50, ny_vol=50)
#' s$volume$estimate
#'
#' @export
solid_xyz3d <- function(
    H1, H2, G1, G2,
    a, b,
    plot = TRUE,
    n_x = 120, n_u = 60, n_v = 60,
    mode = c("faces","wireframe","both"),
    show_faces = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
    colorscales = c("Blues","Blues","Greens","Greens","Reds","Reds"),
    opacities   = 0.35,
    show_surface_grid  = TRUE,
    surface_grid_color = "rgba(60,80,200,0.25)",
    surface_grid_width = 1,
    show_edges = TRUE,
    edge_line  = list(color = "black", width = 2),
    wire_step  = 6,
    wire_line  = list(color = "black", width = 1),
    scene = list(
      aspectmode = "data",
      xaxis = list(title = "x"),
      yaxis = list(title = "y"),
      zaxis = list(title = "z")
    ),
    bg = list(paper = "white", plot = "white"),
    # Volume
    compute_volume = FALSE,
    vol_method = c("adaptive","grid"),
    nx_vol = 300, ny_vol = 300,
    # Slices
    slice = list(x = NULL, y = NULL, z = NULL),
    slice_mode = c("surface","wireframe","both"),
    slice_nx = 200, slice_nu = 120, slice_nv = 120,
    slice_colorscales = list(x = "Oranges", y = "Purples", z = "Greens"),
    slice_opacity = 0.55,
    slice_show_grid  = TRUE,
    slice_grid_color = "rgba(80,80,80,0.25)",
    slice_grid_width = 1,
    slice_wire_step = 8,
    slice_wire_line = list(color = "black", width = 2, dash = "dot")
)
{
  # ---------- helpers ----------
  assert_fun <- function(f, name) if (!is.function(f)) stop(sprintf("'%s' must be a function.", name), call. = FALSE)
  assert_num <- function(x, name) if (!is.numeric(x) || length(x)!=1L || !is.finite(x)) stop(sprintf("'%s' must be a finite numeric scalar.", name), call. = FALSE)
  assert_pos_int <- function(x, name) {
    if (!is.numeric(x) || length(x)!=1L || !is.finite(x) ||
        x < 1 || abs(x - round(x)) > .Machine$double.eps^0.5)
      stop(sprintf("'%s' must be a positive integer.", name), call. = FALSE)
  }
  as_vec6_exact <- function(x, name) {
    if (length(x) == 1L) rep(x, 6)
    else if (length(x) == 6L) x
    else stop(sprintf("'%s' must have length 1 or 6.", name), call. = FALSE)
  }

  is_color <- function(x) {
    if (!is.character(x) || length(x) != 1) return(FALSE)
    if (grepl("^rgba?\\(", x, ignore.case = TRUE)) return(TRUE)
    if (grepl("^#([0-9A-Fa-f]{6}|[0-9A-Fa-f]{8})$", x)) return(TRUE)
    ok <- TRUE; tryCatch(grDevices::col2rgb(x), error = function(...) ok <<- FALSE); ok
  }
  as_rgba <- function(col, alpha = NULL) {
    if (grepl("^rgba?\\(", col, ignore.case = TRUE)) return(col)
    rgb <- grDevices::col2rgb(col) # 0..255 integers
    a <- if (is.null(alpha)) 1 else max(0, min(1, alpha))
    sprintf("rgba(%d,%d,%d,%g)", rgb[1], rgb[2], rgb[3], a)
  }
  is_plotly_scale_list <- function(x) {
    is.list(x) && length(x) >= 2 && all(vapply(x, function(it)
      is.list(it) && length(it) == 2 && is.numeric(it[[1]]) && is.character(it[[2]]) && length(it[[2]])==1L, logical(1)))
  }
  as_colorscale <- function(x, alpha = NULL) {
    if (is.null(x)) return(NULL)
    if (is_plotly_scale_list(x)) return(x)
    if (is.character(x) && length(x) == 1) {
      if (is_color(x)) {
        ccol <- as_rgba(x, alpha); return(list(list(0, ccol), list(1, ccol)))
      } else return(x) # Plotly named scale
    }
    if (is.character(x) && length(x) > 1) {
      cols <- vapply(x, as_rgba, character(1), alpha = alpha)
      pos  <- seq(0, 1, length.out = length(cols))
      return(Map(function(p,c) list(p, c), pos, cols))
    }
    stop("Unrecognized 'colorscale' format.", call. = FALSE)
  }

  # ---------- validations ----------
  assert_fun(H1,"H1"); assert_fun(H2,"H2"); assert_fun(G1,"G1"); assert_fun(G2,"G2")
  assert_num(a,"a"); assert_num(b,"b"); if (b <= a) stop("'b' must be > 'a'.", call. = FALSE)
  for (nm in c("n_x","n_u","n_v","wire_step","nx_vol","ny_vol","slice_nx","slice_nu","slice_nv","slice_wire_step")) {
    assert_pos_int(get(nm), nm)
  }
  mode <- match.arg(mode); vol_method <- match.arg(vol_method); slice_mode <- match.arg(slice_mode)

  show_faces <- as_vec6_exact(show_faces, "show_faces")
  opacities  <- as_vec6_exact(opacities,  "opacities")

  # ---------- colors (6 faces + slices) ----------
  if (length(colorscales) == 6L) {
    colorscales <- lapply(colorscales, as_colorscale)
  } else {
    cs_global <- as_colorscale(colorscales)
    colorscales <- rep(list(cs_global), 6L)
  }
  slice_cs_x <- as_colorscale(slice_colorscales$x)
  slice_cs_y <- as_colorscale(slice_colorscales$y)
  slice_cs_z <- as_colorscale(slice_colorscales$z)

  # ---------- parametrization ----------
  y_blend <- function(x, u) (1 - u) * H1(x) + u * H2(x)
  z_blend <- function(x, y, v) (1 - v) * G1(x, y) + v * G2(x, y)

  x_seq <- seq(a, b, length.out = n_x)
  u_seq <- seq(0, 1, length.out = n_u)
  v_seq <- seq(0, 1, length.out = n_v)

  # clamp steps to mesh sizes (defensivo)
  wire_step       <- max(1L, min(wire_step, max(n_u, n_v)))
  slice_wire_step <- max(1L, min(slice_wire_step, max(slice_nu, slice_nv)))

  # ---------- drawing helpers ----------
  add_edges <- function(fig, X, Y, Z, line = edge_line) {
    i1 <- 1; i2 <- nrow(X); j1 <- 1; j2 <- ncol(X)
    fig |>
      plotly::add_trace(x=X[i1,], y=Y[i1,], z=Z[i1,], type="scatter3d", mode="lines",
                        line=line, hoverinfo="none", showlegend=FALSE) |>
      plotly::add_trace(x=X[i2,], y=Y[i2,], z=Z[i2,], type="scatter3d", mode="lines",
                        line=line, hoverinfo="none", showlegend=FALSE) |>
      plotly::add_trace(x=X[,j1], y=Y[,j1], z=Z[,j1], type="scatter3d", mode="lines",
                        line=line, hoverinfo="none", showlegend=FALSE) |>
      plotly::add_trace(x=X[,j2], y=Y[,j2], z=Z[,j2], type="scatter3d", mode="lines",
                        line=line, hoverinfo="none", showlegend=FALSE)
  }
  add_wireframe <- function(fig, X, Y, Z, step, line = wire_line) {
    nr <- nrow(X); nc <- ncol(X)
    r_idx <- seq(1, nr, by = max(1L, step))
    c_idx <- seq(1, nc, by = max(1L, step))
    for (i in r_idx) {
      fig <- fig |>
        plotly::add_trace(x=X[i,], y=Y[i,], z=Z[i,], type="scatter3d", mode="lines",
                          line=line, hoverinfo="none", showlegend=FALSE)
    }
    for (j in c_idx) {
      fig <- fig |>
        plotly::add_trace(x=X[,j], y=Y[,j], z=Z[,j], type="scatter3d", mode="lines",
                          line=line, hoverinfo="none", showlegend=FALSE)
    }
    fig
  }
  add_surface <- function(fig, X, Y, Z, cs, op, show_grid, grid_col, grid_w) {
    contours_arg <- if (isTRUE(show_grid)) list(
      x = list(show = TRUE, color = grid_col, width = grid_w),
      y = list(show = TRUE, color = grid_col, width = grid_w),
      z = list(show = FALSE)
    ) else NULL
    fig |>
      plotly::add_surface(
        x = X, y = Y, z = Z,
        colorscale = cs, showscale = FALSE,
        opacity = op, contours = contours_arg
      )
  }

  build_face_x <- function(xfix) {
    U <- seq(0, 1, length.out = n_u)
    V <- seq(0, 1, length.out = n_v)
    X <- matrix(xfix, nrow = n_u, ncol = n_v)
    y_line <- y_blend(xfix, U)
    Y <- matrix(rep(y_line, times = n_v), nrow = n_u)
    Z <- matrix(NA_real_, n_u, n_v)
    for (j in seq_len(n_v)) Z[, j] <- z_blend(xfix, y_line, V[j])
    list(X=X, Y=Y, Z=Z)
  }
  build_face_y <- function(y_fun) {
    V <- seq(0, 1, length.out = n_v)
    X <- matrix(rep(x_seq, each = n_v), nrow = n_v)
    y_line <- y_fun(x_seq)
    Y <- matrix(rep(y_line, each = n_v), nrow = n_v)
    Z <- matrix(NA_real_, n_v, n_x)
    for (j in seq_len(n_x)) Z[, j] <- z_blend(x_seq[j], y_line[j], V)
    list(X=X, Y=Y, Z=Z)
  }
  build_face_z <- function(z_fun) {
    U <- seq(0, 1, length.out = n_u)
    X <- matrix(rep(x_seq, each = n_u), nrow = n_u)
    Y <- matrix(NA_real_, n_u, n_x)
    Z <- matrix(NA_real_, n_u, n_x)
    for (j in seq_len(n_x)) {
      yline <- y_blend(x_seq[j], U)
      Y[, j] <- yline
      Z[, j] <- z_fun(x_seq[j], yline)
    }
    list(X=X, Y=Y, Z=Z)
  }

  # ---------- drawing ----------
  fig <- NULL
  if (isTRUE(plot)) {
    if (!requireNamespace("plotly", quietly = TRUE)) {
      warning("You need 'plotly' installed to draw.", call. = FALSE)
    } else {
      fig <- plotly::plot_ly()

      if (isTRUE(show_faces[1])) { F <- build_face_x(a)
      if (mode != "wireframe") fig <- add_surface(fig, F$X,F$Y,F$Z, colorscales[[1]], opacities[1],
                                                  show_surface_grid, surface_grid_color, surface_grid_width)
      if (mode != "faces")     fig <- add_wireframe(fig, F$X,F$Y,F$Z, wire_step)
      if (isTRUE(show_edges))  fig <- add_edges(fig,     F$X,F$Y,F$Z)
      }
      if (isTRUE(show_faces[2])) { F <- build_face_x(b)
      if (mode != "wireframe") fig <- add_surface(fig, F$X,F$Y,F$Z, colorscales[[2]], opacities[2],
                                                  show_surface_grid, surface_grid_color, surface_grid_width)
      if (mode != "faces")     fig <- add_wireframe(fig, F$X,F$Y,F$Z, wire_step)
      if (isTRUE(show_edges))  fig <- add_edges(fig,     F$X,F$Y,F$Z)
      }
      if (isTRUE(show_faces[3])) { F <- build_face_y(H1)
      if (mode != "wireframe") fig <- add_surface(fig, F$X,F$Y,F$Z, colorscales[[3]], opacities[3],
                                                  show_surface_grid, surface_grid_color, surface_grid_width)
      if (mode != "faces")     fig <- add_wireframe(fig, F$X,F$Y,F$Z, wire_step)
      if (isTRUE(show_edges))  fig <- add_edges(fig,     F$X,F$Y,F$Z)
      }
      if (isTRUE(show_faces[4])) { F <- build_face_y(H2)
      if (mode != "wireframe") fig <- add_surface(fig, F$X,F$Y,F$Z, colorscales[[4]], opacities[4],
                                                  show_surface_grid, surface_grid_color, surface_grid_width)
      if (mode != "faces")     fig <- add_wireframe(fig, F$X,F$Y,F$Z, wire_step)
      if (isTRUE(show_edges))  fig <- add_edges(fig,     F$X,F$Y,F$Z)
      }
      if (isTRUE(show_faces[5])) { F <- build_face_z(function(x,y) G1(x,y))
      if (mode != "wireframe") fig <- add_surface(fig, F$X,F$Y,F$Z, colorscales[[5]], opacities[5],
                                                  show_surface_grid, surface_grid_color, surface_grid_width)
      if (mode != "faces")     fig <- add_wireframe(fig, F$X,F$Y,F$Z, wire_step)
      if (isTRUE(show_edges))  fig <- add_edges(fig,     F$X,F$Y,F$Z)
      }
      if (isTRUE(show_faces[6])) { F <- build_face_z(function(x,y) G2(x,y))
      if (mode != "wireframe") fig <- add_surface(fig, F$X,F$Y,F$Z, colorscales[[6]], opacities[6],
                                                  show_surface_grid, surface_grid_color, surface_grid_width)
      if (mode != "faces")     fig <- add_wireframe(fig, F$X,F$Y,F$Z, wire_step)
      if (isTRUE(show_edges))  fig <- add_edges(fig,     F$X,F$Y,F$Z)
      }

      fig <- fig |>
        plotly::layout(
          title = sprintf("Solid (mode=%s)", mode),
          scene = scene, paper_bgcolor = bg$paper, plot_bgcolor = bg$plot
        )
      if (interactive()) print(fig)
    }
  }

  # ---------- volume ----------
  volume <- NULL
  if (isTRUE(compute_volume)) {
    # thickness always nonnegative
    f_gap  <- function(x, y) abs(G2(x, y) - G1(x, y))
    f_gapv <- Vectorize(f_gap, SIMPLIFY = TRUE)
    H1v <- Vectorize(H1); H2v <- Vectorize(H2)

    if (vol_method == "adaptive") {
      inner <- function(x) {
        yl <- pmin(H1v(x), H2v(x)); yh <- pmax(H1v(x), H2v(x))
        g <- function(y) f_gapv(x, y)
        sapply(seq_along(x), function(i) {
          res <- try(stats::integrate(g, lower = yl[i], upper = yh[i], rel.tol = 1e-6, stop.on.error = FALSE), silent = TRUE)
          if (inherits(res, "try-error") || !is.finite(res$value)) 0 else res$value
        })
      }
      res <- try(stats::integrate(function(x) inner(x), lower = a, upper = b, rel.tol = 1e-6, stop.on.error = FALSE), silent = TRUE)
      vol_val <- if (inherits(res, "try-error") || !is.finite(res$value)) NA_real_ else res$value
      volume  <- list(estimate = vol_val, method = "adaptive")
    } else {
      xs <- seq(a, b, length.out = nx_vol)
      dx <- (b - a) / (nx_vol - 1)
      vol_x <- numeric(nx_vol)
      for (i in seq_along(xs)) {
        xl <- xs[i]; yl <- H1(xl); yh <- H2(xl); if (yh < yl) { tmp <- yl; yl <- yh; yh <- tmp }
        ys <- seq(yl, yh, length.out = ny_vol)
        if (length(ys) < 2) { vol_x[i] <- 0; next }
        dy <- (yh - yl) / (ny_vol - 1)
        vals <- f_gapv(xl, ys)
        Iy <- dy * (sum(vals) - 0.5*vals[1] - 0.5*vals[length(vals)])
        vol_x[i] <- Iy
      }
      vol_val <- dx * (sum(vol_x) - 0.5*vol_x[1] - 0.5*vol_x[length(vol_x)])
      volume  <- list(estimate = vol_val, method = "grid", nx = nx_vol, ny = ny_vol)
    }
  }

  list(x_seq = x_seq, u_seq = u_seq, v_seq = v_seq, fig = fig, volume = volume)
}

