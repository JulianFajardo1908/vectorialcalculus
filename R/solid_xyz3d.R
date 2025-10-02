#' Solid \eqn{(x,y,z)} with \eqn{y \in [H_1(x), H_2(x)]} and \eqn{z \in [G_1(x,y), G_2(x,y)]}
#'
#' Builds and (optionally) renders with \pkg{plotly} the solid
#' \deqn{\{(x,y,z):\ x\in[a,b],\ y\in[H_1(x),H_2(x)],\ z\in[G_1(x,y),G_2(x,y)]\}}
#' using a curvilinear-prism parametrization. It supports:
#' \itemize{
#'   \item \strong{Display modes}: \code{mode = "faces" | "wireframe" | "both"}.
#'   \item \strong{Numerical volume}: \code{compute_volume = TRUE} with \code{vol_method = "adaptive"|"grid"}.
#'   \item \strong{Internal slices} on planes \eqn{x=x_0}, \eqn{y=y_0}, \eqn{z=z_0}.
#'   \item \strong{Flexible colors}: flat colors (R names, \code{"#RRGGBB"}, \code{"rgba(...)"});
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
#'   in the order \code{c("x=a","x=b","y=H1","y=H2","z=G1","z=G2")}.
#' @param colorscales Color scales for faces. Accepts:
#'   \itemize{
#'     \item A single Plotly scale name (e.g. \code{"Blues"}, \code{"Viridis"}) → applied to all faces.
#'     \item A single flat color (R name, \code{"#RRGGBB"}, \code{"rgba(...)"}) → flat scale.
#'     \item A color vector \code{c("white","#2a9d8f",...)} → evenly spaced gradient for all faces.
#'     \item A vector/list of 6 scales (one per face) in any of the formats above.
#'   }
#' @param opacities Opacities for faces (length 1 or 6).
#' @param show_surface_grid Logical. Draws a grid over surfaces.
#' @param surface_grid_color,surface_grid_width Grid style for surfaces.
#' @param show_edges Logical. Draw edges of each face.
#' @param edge_line Edge style (\code{list(color, width, dash)}).
#' @param wire_step Integer \eqn{\ge 1}: draw every \code{wire_step}-th mesh line in wireframe mode.
#' @param wire_line Wireframe line style (\code{list(color, width, dash)}).
#' @param scene 3D scene options (default \code{aspectmode="data"}).
#' @param bg Background colors: \code{list(paper="white", plot="white")}.
#' @param compute_volume Logical. If \code{TRUE}, estimates \eqn{\int_a^b \int_{H_1}^{H_2} (G_2-G_1)\,dy\,dx}.
#' @param vol_method Volume integration method: \code{"adaptive"} (nested \code{stats::integrate})
#'   or \code{"grid"} (trapezoidal rule on a regular grid).
#' @param nx_vol,ny_vol Grid sizes for \code{vol_method="grid"}.
#' @param slice List of slices: \code{slice = list(x = NULL, y = NULL, z = NULL)}. Each entry
#'   can be a number or numeric vector.
#' @param slice_mode Slice rendering mode: \code{"surface"}, \code{"wireframe"} or \code{"both"}.
#' @param slice_nx,slice_nu,slice_nv Mesh resolution for slices (\eqn{x}: \code{nu×nv}, \eqn{y}: \code{nx×nv}, \eqn{z}: \code{nx×nu}).
#' @param slice_colorscales Color scales for slices: \code{list(x=, y=, z=)} in the same formats as \code{colorscales}.
#' @param slice_opacity Opacity for slices (0–1).
#' @param slice_show_grid,slice_grid_color,slice_grid_width Grid options for slices.
#' @param slice_wire_step,slice_wire_line Wireframe step and style for slices.
#'
#' @details
#' \strong{About color scales}: Each \emph{colorscale} can be
#' \enumerate{
#'   \item A Plotly scale name (\code{"Blues"}, \code{"Viridis"}, …),
#'   \item A single color (R name, \code{"#RRGGBB"}, \code{"rgba(r,g,b,a)"}),
#'   \item A color vector (\code{c("white","#2a9d8f",...)}) → evenly spaced gradient,
#'   \item A ready Plotly list (pairs \code{list(list(pos,color), ...)}).
#' }
#' If \code{colorscales} has length 1 (or is a color vector), it is applied to all faces.
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
#' \dontshow{if (interactive()) \{}
#' H1 <- function(x) -1 - x
#' H2 <- function(x)  1 - x^2
#' G1 <- function(x,y) y
#' G2 <- function(x,y) y + 1
#'
#' # Six flat face colors
#' solid_xyz3d(
#'   H1,H2,G1,G2, a=-1, b=1, plot=TRUE, mode="faces",
#'   colorscales = c("#b3cde0","#ccebc5","#decbe4","#fed9a6","#ffffcc","#fbb4ae"),
#'   opacities   = 0.30
#' )
#'
#' # Global custom gradient + sparse wireframe
#' solid_xyz3d(
#'   H1,H2,G1,G2, a=-1, b=1, plot=TRUE, mode="both",
#'   colorscales = c("white","#2a9d8f"),
#'   opacities   = 0.25,
#'   wire_step   = 12,
#'   wire_line   = list(color="rgba(20,20,20,0.25)", width=1, dash="dot"),
#'   edge_line   = list(color="rgba(30,30,30,0.35)", width=1),
#'   slice = list(x = 0, z = 0.6),
#'   slice_mode = "surface",
#'   slice_colorscales = list(x = "#ffb703", y = "#9b5de5", z = c("white", "#2a9d8f")),
#'   slice_opacity = 0.6
#' )
#' \dontshow{\}}
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
) {
  # ---------- local helpers ----------
  assert_fun <- function(f, name) if (!is.function(f)) stop(sprintf("'%s' must be a function.", name), call. = FALSE)
  assert_num <- function(x, name) if (!is.numeric(x) || length(x)!=1L || !is.finite(x)) stop(sprintf("'%s' must be a finite numeric scalar.", name), call. = FALSE)
  as_vec6 <- function(x) if (length(x)==1L) rep(x, 6) else x

  is_color <- function(x) {
    if (!is.character(x) || length(x) != 1) return(FALSE)
    if (grepl("^rgba?\\(", x, ignore.case = TRUE)) return(TRUE)
    if (grepl("^#([0-9A-Fa-f]{6}|[0-9A-Fa-f]{8})$", x)) return(TRUE)
    ok <- TRUE; tryCatch(grDevices::col2rgb(x), error = function(...) ok <<- FALSE); ok
  }
  as_rgba <- function(col, alpha = NULL) {
    if (grepl("^rgba?\\(", col, ignore.case = TRUE)) return(col)
    rgb <- grDevices::col2rgb(col) / 255
    a <- if (is.null(alpha)) 1 else max(0, min(1, alpha))
    sprintf("rgba(%g,%g,%g,%g)", 255*rgb[1], 255*rgb[2], 255*rgb[3], a)
  }
  as_colorscale <- function(x, alpha = NULL) {
    if (is.null(x)) return(NULL)
    if (is.list(x) && length(x) >= 2 && is.numeric(x[[1]][[1]])) return(x)
    if (is.character(x) && length(x) == 1) {
      if (is_color(x)) {
        ccol <- as_rgba(x, alpha); return(list(list(0, ccol), list(1, ccol)))
      } else return(x)
    }
    if (is.character(x) && length(x) > 1) {
      cols <- vapply(x, as_rgba, character(1), alpha = alpha)
      pos  <- seq(0, 1, length.out = length(cols))
      return(lapply(seq_along(cols), function(i) list(pos[i], cols[i])))
    }
    stop("Unrecognized 'colorscale' format.", call. = FALSE)
  }

  # ---------- validations ----------
  assert_fun(H1,"H1"); assert_fun(H2,"H2"); assert_fun(G1,"G1"); assert_fun(G2,"G2")
  assert_num(a,"a"); assert_num(b,"b"); if (b <= a) stop("'b' must be > 'a'.", call. = FALSE)
  for (nm in c("n_x","n_u","n_v","wire_step","nx_vol","ny_vol","slice_nx","slice_nu","slice_nv","slice_wire_step")) {
    assert_num(get(nm), nm); if (get(nm) < 1) stop(nm," must be >= 1.", call. = FALSE)
  }
  mode <- match.arg(mode); vol_method <- match.arg(vol_method); slice_mode <- match.arg(slice_mode)

  # ---------- colors (6 faces + slices) ----------
  colorscales <- as_vec6(colorscales)
  if (length(colorscales) != 6) {
    cs_global <- as_colorscale(colorscales)
    colorscales <- rep(list(cs_global), 6)
  } else {
    colorscales <- lapply(colorscales, as_colorscale)
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

      # ---------- slices ----------
      has_x <- !is.null(slice$x) && length(slice$x) > 0
      has_y <- !is.null(slice$y) && length(slice$y) > 0
      has_z <- !is.null(slice$z) && length(slice$z) > 0

      if (has_x) {
        U <- seq(0,1,length.out = slice_nu)
        V <- seq(0,1,length.out = slice_nv)
        for (x0 in as.numeric(slice$x)) {
          y_line <- y_blend(x0, U)
          Xs <- matrix(x0, nrow = slice_nu, ncol = slice_nv)
          Ys <- matrix(rep(y_line, times = slice_nv), nrow = slice_nu)
          Zs <- matrix(NA_real_, slice_nu, slice_nv)
          for (j in seq_len(slice_nv)) Zs[, j] <- z_blend(x0, y_line, V[j])

          if (slice_mode != "wireframe") {
            fig <- add_surface(fig, Xs, Ys, Zs,
                               slice_cs_x, slice_opacity,
                               slice_show_grid, slice_grid_color, slice_grid_width)
          }
          if (slice_mode != "surface") {
            fig <- add_wireframe(fig, Xs, Ys, Zs, slice_wire_step, line = slice_wire_line)
            fig <- add_edges(fig,     Xs, Ys, Zs, line = slice_wire_line)
          }
        }
      }
      if (has_y) {
        V <- seq(0,1,length.out = slice_nv)
        xs <- seq(a, b, length.out = slice_nx)
        for (y0 in as.numeric(slice$y)) {
          Xs <- matrix(rep(xs, each = slice_nv), nrow = slice_nv)
          Ys <- matrix(y0, nrow = slice_nv, ncol = slice_nx)
          Zs <- matrix(NA_real_, slice_nv, slice_nx)
          for (j in seq_len(slice_nx)) Zs[, j] <- z_blend(xs[j], y0, V)

          if (slice_mode != "wireframe") {
            fig <- add_surface(fig, Xs, Ys, Zs,
                               slice_cs_y, slice_opacity,
                               slice_show_grid, slice_grid_color, slice_grid_width)
          }
          if (slice_mode != "surface") {
            fig <- add_wireframe(fig, Xs, Ys, Zs, slice_wire_step, line = slice_wire_line)
            fig <- add_edges(fig,     Xs, Ys, Zs, line = slice_wire_line)
          }
        }
      }
      if (has_z) {
        xs <- seq(a, b, length.out = slice_nx)
        U  <- seq(0, 1, length.out = slice_nu)
        for (z0 in as.numeric(slice$z)) {
          Xs <- matrix(rep(xs, each = slice_nu), nrow = slice_nu)
          Ys <- matrix(NA_real_, slice_nu, slice_nx)
          Zs <- matrix(z0,       slice_nu, slice_nx)
          for (j in seq_len(slice_nx)) {
            y_line <- y_blend(xs[j], U)
            g1 <- G1(xs[j], y_line); g2 <- G2(xs[j], y_line)
            inside <- z0 >= pmin(g1,g2) & z0 <= pmax(g1,g2)
            Ys[, j] <- ifelse(inside, y_line, NA_real_)
            Zs[, j] <- ifelse(inside, z0,     NA_real_)
          }

          if (slice_mode != "wireframe") {
            fig <- add_surface(fig, Xs, Ys, Zs,
                               slice_cs_z, slice_opacity,
                               slice_show_grid, slice_grid_color, slice_grid_width)
          }
          if (slice_mode != "surface") {
            fig <- add_wireframe(fig, Xs, Ys, Zs, slice_wire_step, line = slice_wire_line)
            fig <- add_edges(fig,     Xs, Ys, Zs, line = slice_wire_line)
          }
        }
      }

      fig <- fig |>
        plotly::layout(
          title = sprintf("Solid (mode=%s)%s", mode, if (has_x||has_y||has_z) " with slices" else ""),
          scene = scene, paper_bgcolor = bg$paper, plot_bgcolor = bg$plot
        )
      print(fig)
    }
  }

  # ---------- volume ----------
  volume <- NULL
  if (isTRUE(compute_volume)) {
    f_gap  <- function(x, y) G2(x, y) - G1(x, y)
    f_gapv <- Vectorize(f_gap, SIMPLIFY = TRUE)
    H1v <- Vectorize(H1); H2v <- Vectorize(H2)

    if (vol_method == "adaptive") {
      inner <- function(x) {
        yl <- pmin(H1v(x), H2v(x)); yh <- pmax(H1v(x), H2v(x))
        g <- function(y) f_gapv(x, y)
        sapply(seq_along(x), function(i) {
          stats::integrate(g, lower = yl[i], upper = yh[i], rel.tol = 1e-6)$value
        })
      }
      vol_val <- stats::integrate(function(x) inner(x), lower = a, upper = b, rel.tol = 1e-6)$value
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
