#' Cylindrical solid defined by radial and vertical bounds (with optional plot)
#'
#' Builds (and optionally plots with \pkg{plotly}) the solid
#' \deqn{\{(x,y,z):\ \theta \in [\theta_{\min}, \theta_{\max}],\ r \in [R_1(\theta), R_2(\theta)],\ z \in [Z_1(r,\theta), Z_2(r,\theta)]\}.}
#' Rendering is done by sampling a curvilinear grid in \eqn{(\theta, u, v)},
#' where \eqn{u,v \in [0,1]} linearly blend the inner/outer radius and lower/upper \eqn{z}, respectively.
#'
#' The volume is computed as
#' \deqn{\int_{\theta_{\min}}^{\theta_{\max}} \int_{R_1(\theta)}^{R_2(\theta)} \int_{Z_1(r,\theta)}^{Z_2(r,\theta)} r\,dz\,dr\,d\theta.}
#'
#' @param R1,R2 Functions \code{function(theta)} giving the inner/outer radius bounds.
#' @param Z1,Z2 Functions \code{function(r,theta)} giving the lower/upper \code{z} bounds.
#' @param th_min,th_max Angular limits (numeric scalars) with \code{th_max > th_min}.
#' @param plot Logical. If \code{TRUE}, plot with \pkg{plotly}.
#' @param n_theta,n_u,n_v Mesh resolution in \code{theta} (angle), \code{u} (radial blend), \code{v} (vertical blend).
#' @param mode \code{"faces"}, \code{"wireframe"} or \code{"both"}.
#' @param colorscale Plotly colorscale (name, single color, or vector of colors) for the surface.
#' @param opacity Surface opacity (0-1).
#' @param show_surface_grid Logical. Draw a grid on the surface.
#' @param surface_grid_color,surface_grid_width Grid aesthetics.
#' @param edge_line,wire_line Line styles for edges and wireframe.
#' @param scene,bg Plotly 3D scene and background colors.
#' @param compute_volume Logical. If \code{TRUE}, compute the volume.
#' @param vol_method \code{"adaptive"} (nested \code{stats::integrate}) or \code{"grid"} (trapezoidal on regular mesh).
#' @param ntheta_vol,nr_vol Mesh sizes for \code{vol_method = "grid"}.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{theta_seq}, \code{u_seq}, \code{v_seq}: the parameter sequences used,
#'   \item \code{fig}: a \pkg{plotly} object if \code{plot = TRUE}, otherwise \code{NULL},
#'   \item \code{volume}: \code{NULL} or a list with \code{estimate} and metadata when \code{compute_volume = TRUE}.
#' }
#'
#' @examples
#' \dontrun{
#' # Example: a quarter-twisted "cup": R in [0, 1+0.2*cos(theta)], z in [0, 1 + 0.5*r]
#' R1 <- function(theta) 0
#' R2 <- function(theta) 1 + 0.2*cos(theta)
#' Z1 <- function(r,theta) 0
#' Z2 <- function(r,theta) 1 + 0.5*r
#' solid_cylindrical3d(
#'   R1, R2, Z1, Z2, th_min = 0, th_max = pi/2,
#'   plot = TRUE, mode = "both",
#'   colorscale = c("white", "#2a9d8f"), opacity = 0.35, show_surface_grid = TRUE,
#'   compute_volume = TRUE, vol_method = "adaptive"
#' )$volume
#' }
#' @export
#' @importFrom stats integrate
solid_cylindrical3d <- function(
    R1, R2, Z1, Z2,
    th_min, th_max,
    plot = TRUE,
    n_theta = 160, n_u = 70, n_v = 70,
    mode = c("faces","wireframe","both"),
    colorscale = "Blues",
    opacity    = 0.35,
    show_surface_grid  = TRUE,
    surface_grid_color = "rgba(60,80,200,0.25)",
    surface_grid_width = 1,
    edge_line = list(color = "black", width = 2),
    wire_line = list(color = "rgba(0,0,0,0.35)", width = 1),
    scene = list(
      aspectmode = "data",
      xaxis = list(title = "x"),
      yaxis = list(title = "y"),
      zaxis = list(title = "z")
    ),
    bg = list(paper = "white", plot = "white"),
    compute_volume = FALSE,
    vol_method = c("adaptive","grid"),
    ntheta_vol = 400, nr_vol = 400
) {
  mode <- match.arg(mode); vol_method <- match.arg(vol_method)
  if (!is.function(R1) || !is.function(R2)) stop("'R1','R2' must be functions of theta.", call. = FALSE)
  if (!is.function(Z1) || !is.function(Z2)) stop("'Z1','Z2' must be functions of (r,theta).", call. = FALSE)
  if (!is.numeric(th_min) || !is.numeric(th_max) || th_max <= th_min)
    stop("'th_max' must be > 'th_min'.", call. = FALSE)

  # helpers (colorscale parsing)
  is_color <- function(x) {
    if (!is.character(x) || length(x) != 1) return(FALSE)
    if (grepl("^rgba?\\(", x, ignore.case = TRUE)) return(TRUE)
    if (grepl("^#([0-9A-Fa-f]{6}|[0-9A-Fa-f]{8})$", x)) return(TRUE)
    ok <- TRUE; tryCatch(grDevices::col2rgb(x), error = function(...) ok <<- FALSE); ok
  }
  as_rgba <- function(col, alpha = NULL) {
    if (grepl("^rgba?\\(", col, ignore.case = TRUE)) return(col)
    rgb <- grDevices::col2rgb(col) / 255
    a   <- if (is.null(alpha)) 1 else max(0, min(1, alpha))
    sprintf("rgba(%g,%g,%g,%g)", 255*rgb[1], 255*rgb[2], 255*rgb[3], a)
  }
  as_colorscale <- function(x, alpha = NULL) {
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

  # parameter grids
  theta_seq <- seq(th_min, th_max, length.out = n_theta)
  u_seq     <- seq(0, 1, length.out = n_u)   # radius blend
  v_seq     <- seq(0, 1, length.out = n_v)   # z blend

  # Build three pairs of faces to cover the surface:
  # 1) theta = const (two faces) -> vary u (r) and v (z)
  build_face_theta <- function(th) {
    rr_in  <- R1(th); rr_out <- R2(th)
    r_line <- (1 - u_seq) * rr_in + u_seq * rr_out
    # matrices n_u x n_v
    R <- matrix(rep(r_line, times = n_v), nrow = n_u)
    Z <- matrix(NA_real_, n_u, n_v)
    for (j in seq_len(n_v)) {
      z_low  <- Z1(r_line, th)
      z_high <- Z2(r_line, th)
      Z[, j] <- (1 - v_seq[j]) * z_low + v_seq[j] * z_high
    }
    X <- R * cos(th); Y <- R * sin(th)
    list(X=X, Y=Y, Z=Z)
  }

  # 2) r = R1(theta) and r = R2(theta): sweep theta,v, radius fixed
  build_face_r <- function(which = c("inner","outer")) {
    which <- match.arg(which)
    Rfix  <- if (which=="inner") vapply(theta_seq, R1, numeric(1)) else vapply(theta_seq, R2, numeric(1))
    R <- matrix(rep(Rfix, each = n_v), nrow = n_v)
    TH <- matrix(rep(theta_seq, each = n_v), nrow = n_v)
    Z  <- matrix(NA_real_, n_v, n_theta)
    for (j in seq_len(n_theta)) {
      z_low  <- Z1(Rfix[j], theta_seq[j])
      z_high <- Z2(Rfix[j], theta_seq[j])
      Z[, j] <- (1 - v_seq) * z_low + v_seq * z_high
    }
    X <- R * cos(TH); Y <- R * sin(TH)
    list(X=X, Y=Y, Z=Z)
  }

  # 3) z = Z1(r,theta) and z = Z2(r,theta): sweep theta,u, set z to bound
  build_face_z <- function(which = c("low","high")) {
    which <- match.arg(which)
    TH <- matrix(rep(theta_seq, each = n_u), nrow = n_u)
    rr_in  <- vapply(theta_seq, R1, numeric(1))
    rr_out <- vapply(theta_seq, R2, numeric(1))
    R <- matrix(NA_real_, n_u, n_theta)
    for (j in seq_len(n_theta)) R[, j] <- (1 - u_seq) * rr_in[j] + u_seq * rr_out[j]
    Z <- matrix(NA_real_, n_u, n_theta)
    for (j in seq_len(n_theta)) {
      z_low  <- Z1(R[,j], theta_seq[j])
      z_high <- Z2(R[,j], theta_seq[j])
      Z[, j] <- if (which=="low") z_low else z_high
    }
    X <- R * cos(TH); Y <- R * sin(TH)
    list(X=X, Y=Y, Z=Z)
  }

  add_edges <- function(fig, F, line = edge_line) {
    X <- F$X; Y <- F$Y; Z <- F$Z
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
  add_wireframe <- function(fig, F, line = wire_line, step = 6) {
    X <- F$X; Y <- F$Y; Z <- F$Z
    nr <- nrow(X); nc <- ncol(X)
    r_idx <- seq(1, nr, by = max(1L, step))
    c_idx <- seq(1, nc, by = max(1L, step))
    for (i in r_idx) fig <- fig |>
      plotly::add_trace(x=X[i,], y=Y[i,], z=Z[i,], type="scatter3d", mode="lines",
                        line=line, hoverinfo="none", showlegend=FALSE)
    for (j in c_idx) fig <- fig |>
      plotly::add_trace(x=X[,j], y=Y[,j], z=Z[,j], type="scatter3d", mode="lines",
                        line=line, hoverinfo="none", showlegend=FALSE)
    fig
  }
  add_surface <- function(fig, F) {
    contours_arg <- if (isTRUE(show_surface_grid)) list(
      x = list(show = TRUE, color = surface_grid_color, width = surface_grid_width),
      y = list(show = TRUE, color = surface_grid_color, width = surface_grid_width),
      z = list(show = FALSE)
    ) else NULL
    fig |>
      plotly::add_surface(
        x = F$X, y = F$Y, z = F$Z,
        colorscale = as_colorscale(colorscale),
        showscale  = FALSE,
        opacity    = opacity,
        contours   = contours_arg
      )
  }

  fig <- NULL
  if (isTRUE(plot)) {
    if (!requireNamespace("plotly", quietly = TRUE)) {
      warning("You need 'plotly' installed to plot.", call. = FALSE)
    } else {
      fig <- plotly::plot_ly()

      # theta faces
      Fth_min <- build_face_theta(th_min)
      Fth_max <- build_face_theta(th_max)
      for (Fth in list(Fth_min, Fth_max)) {
        if (mode != "wireframe") fig <- add_surface(fig, Fth)
        if (mode != "faces")     fig <- add_wireframe(fig, Fth)
        fig <- add_edges(fig, Fth)
      }

      # inner/outer radius faces
      Fin  <- build_face_r("inner")
      Fout <- build_face_r("outer")
      for (Fr in list(Fin, Fout)) {
        if (mode != "wireframe") fig <- add_surface(fig, Fr)
        if (mode != "faces")     fig <- add_wireframe(fig, Fr)
        fig <- add_edges(fig, Fr)
      }

      # z=Z1, z=Z2 faces
      Fz1 <- build_face_z("low")
      Fz2 <- build_face_z("high")
      for (Fz in list(Fz1, Fz2)) {
        if (mode != "wireframe") fig <- add_surface(fig, Fz)
        if (mode != "faces")     fig <- add_wireframe(fig, Fz)
        fig <- add_edges(fig, Fz)
      }

      fig <- fig |>
        plotly::layout(
          title = "Cylindrical solid",
          scene = scene, paper_bgcolor = bg$paper, plot_bgcolor = bg$plot
        )
      print(fig)
    }
  }

  # ----- Volume -------------------------------------------------------
  volume <- NULL
  if (isTRUE(compute_volume)) {
    R1v <- Vectorize(R1); R2v <- Vectorize(R2)
    Z1v <- Vectorize(Z1); Z2v <- Vectorize(Z2)

    if (vol_method == "adaptive") {
      inner_r <- function(theta) {
        rL <- R1v(theta); rU <- R2v(theta)
        if (rU < rL) { tmp <- rL; rL <- rU; rU <- tmp }
        if (!is.finite(rL) || !is.finite(rU) || rU <= rL) return(0)
        g <- function(r) (Z2v(r, theta) - Z1v(r, theta)) * r
        stats::integrate(g, lower = rL, upper = rU, rel.tol = 1e-6)$value
      }
      val <- stats::integrate(function(th) inner_r(th), lower = th_min, upper = th_max, rel.tol = 1e-6)$value
      volume <- list(estimate = val, method = "adaptive")
    } else {
      ths <- seq(th_min, th_max, length.out = ntheta_vol)
      dth <- (th_max - th_min) / (ntheta_vol - 1)
      accum_th <- numeric(ntheta_vol)
      for (i in seq_along(ths)) {
        th <- ths[i]
        rL <- R1v(th); rU <- R2v(th)
        if (rU < rL) { tmp <- rL; rL <- rU; rU <- tmp }
        if (!is.finite(rL) || !is.finite(rU) || rU <= rL) { accum_th[i] <- 0; next }
        rs <- seq(rL, rU, length.out = nr_vol)
        if (length(rs) < 2) { accum_th[i] <- 0; next }
        dr <- (rU - rL) / (nr_vol - 1)
        vals <- (Z2v(rs, th) - Z1v(rs, th)) * rs
        Ir <- dr * (sum(vals) - 0.5*vals[1] - 0.5*vals[length(vals)])
        accum_th[i] <- Ir
      }
      val <- dth * (sum(accum_th) - 0.5*accum_th[1] - 0.5*accum_th[length(accum_th)])
      volume <- list(estimate = val, method = "grid", ntheta = ntheta_vol, nr = nr_vol)
    }
  }

  list(theta_seq = theta_seq, u_seq = u_seq, v_seq = v_seq, fig = fig, volume = volume)
}

