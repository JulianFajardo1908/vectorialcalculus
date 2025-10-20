#' Vector field + streamline in 3D (single combined figure)
#'
#' Plots the 3D vector field \eqn{\mathbf F(x,y,z)} over a curvilinear volume and
#' overlays the \strong{streamline} \eqn{\mathbf r'(t)=\mathbf F(\mathbf r(t),t)}
#' with \eqn{\mathbf r(0)=\mathbf p}, integrated by fixed-step \strong{RK4} in a
#' \emph{single} \pkg{plotly} figure.
#'
#' The volume is defined by linear blends:
#' \eqn{
#' \begin{aligned}
#' x &= (1-i)\,a + i\,b,\quad i\in[0,1],\\
#' y &= (1-j)\,H_1(x) + j\,H_2(x),\quad j\in[0,1],\\
#' z &= (1-k)\,G_1(x,y) + k\,G_2(x,y),\quad k\in[0,1].
#' \end{aligned}
#' }
#'
#'
#' @param F \code{function(x,y,z)} or \code{function(x,y,z,t)} returning \code{c(Fx,Fy,Fz)}.
#'   If it includes \code{t}, the time dependence is used in the integration.
#' @param H1,H2 \code{function(x)}: lower/upper limits in \eqn{y}.
#' @param G1,G2 \code{function(x,y)}: lower/upper limits in \eqn{z}.
#' @param a,b Domain limits in \eqn{x} (\code{b > a}).
#' @param NX,NY,NZ Integers \eqn{\ge 1}: sampling density of the field in \eqn{i,j,k}
#' @param p Streamline initial point \code{c(x0,y0,z0)}.
#' @param T,step Final time and RK4 step. If \code{T<0}, integrates backward.
#'
#' @param arrows Arrow mode: \code{"none"}, \code{"line"}, \code{"cone"} or \code{"both"}.
#' @param arrow_scale Global arrow length scale (fraction of the largest span of the bounding box).
#' @param normalize_bias Saturation bias for the field magnitude:
#'   uses \eqn{\tilde F = F / \sqrt{\|F\|^2 + \texttt{bias}}}. Default \code{1}.
#' @param normal_color,normal_width Color and width of the arrow \emph{shaft} (line).
#' @param arrow_color,arrow_opacity,arrow_size Color, opacity and relative size of the head
#'   (cone if available; otherwise chevron).
#'
#' @param traj_color,traj_width Streamline color and width.
#' @param traj_markers,traj_marker_size Whether to show markers on the streamline and their size.
#'
#' @param scene List of 3D scene settings (e.g., \code{aspectmode="data"}, axis titles).
#' @param bg Background colors: list with \code{paper} and \code{plot}.
#'
#' @details
#' - Arrow length scales with a \emph{saturated} version of \eqn{\|\mathbf F\|}:
#'   normalize with \code{normalize_bias}, then \eqn{L \propto \|\tilde{\mathbf F}\|}.
#' - If the plot looks cluttered, reduce \code{NX,NY,NZ} or \code{arrow_scale}.
#' - If the streamline does not stand out, increase \code{traj_width} or change \code{traj_color}.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{field_points}: \code{data.frame} with \code{x,y,z,mag} at arrow bases,
#'   \item \code{field_segments}: \code{data.frame} with \code{x0,y0,z0,x1,y1,z1} for arrow shafts,
#'   \item \code{traj}: \code{data.frame} with the trajectory (\code{t,x,y,z,speed}),
#'   \item \code{fig}: \pkg{plotly} object with field + streamline.
#' }
#'
#' @examples
#' \dontshow{if (interactive()) \{}
#' H1 <- function(x) -1; H2 <- function(x) 1
#' G1 <- function(x,y) -0.5; G2 <- function(x,y) 0.5
#' F  <- function(x,y,z) c(-y, x, 0.6)
#' streamline_and_field3d(
#'   F, H1, H2, G1, G2,
#'   a = -2, b = 2, NX = 10, NY = 8, NZ = 5,
#'   p = c(1,0,0), T = 12, step = 0.02,
#'   arrows = "both", arrow_scale = 0.12, normalize_bias = 1,
#'   normal_color = "rgba(0,0,0,0.55)", normal_width = 2,
#'   arrow_color  = "#1d3557", arrow_opacity = 0.95, arrow_size = 0.4,
#'   traj_color   = "#e63946", traj_width = 5, traj_markers = TRUE
#' )
#' \dontshow{\}}
#'
#' @export
streamline_and_field3d <- function(
    F, H1, H2, G1, G2,
    a, b, NX = 8, NY = 6, NZ = 6,
    p, T, step,
    # field style
    arrows = c("both","line","cone","none"),
    arrow_scale = 0.08,
    normalize_bias = 1,
    normal_color = "rgba(0,0,0,0.55)",
    normal_width = 2,
    arrow_color  = "#1d3557",
    arrow_opacity = 0.95,
    arrow_size    = 0.35,
    # streamline style
    traj_color = "#e63946",
    traj_width = 5,
    traj_markers = TRUE,
    traj_marker_size = 2,
    # scene / backgrounds
    scene = list(
      aspectmode = "data",
      xaxis = list(title = "x"),
      yaxis = list(title = "y"),
      zaxis = list(title = "z")
    ),
    bg = list(paper = "white", plot = "white")
) {
  arrows <- match.arg(arrows)
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("Package 'plotly' is required for plotting.", call. = FALSE)
  }
  # --- basic checks
  if (!is.function(F))  stop("'F' must be function(x,y,z) or function(x,y,z,t).", call. = FALSE)
  if (!is.function(H1) || !is.function(H2)) stop("'H1' and 'H2' must be functions of x.", call. = FALSE)
  if (!is.function(G1) || !is.function(G2)) stop("'G1' and 'G2' must be functions of (x,y).", call. = FALSE)
  if (!is.numeric(a) || !is.numeric(b) || length(a)!=1L || length(b)!=1L || !is.finite(a) || !is.finite(b) || b <= a)
    stop("'a' and 'b' must be scalars with b > a.", call. = FALSE)
  for (nm in c("NX","NY","NZ")) {
    v <- get(nm); if (!is.numeric(v) || length(v)!=1L || !is.finite(v) || v < 1)
      stop(sprintf("'%s' must be integer >= 1.", nm), call. = FALSE)
  }
  if (!is.numeric(p) || length(p) != 3L || any(!is.finite(p)))
    stop("'p' must be c(x0,y0,z0).", call. = FALSE)
  if (!is.numeric(T) || length(T) != 1L || !is.finite(T))
    stop("'T' must be a numeric scalar.", call. = FALSE)
  if (!is.numeric(step) || length(step) != 1L || !is.finite(step) || step <= 0)
    stop("'step' must be > 0.", call. = FALSE)

  # --- helpers
  has_cone <- "add_cone" %in% getNamespaceExports("plotly")
  perp_unit <- function(n) {
    if (abs(n[1]) <= abs(n[2]) && abs(n[1]) <= abs(n[3])) v <- c(0, -n[3], n[2])
    else if (abs(n[2]) <= abs(n[1]) && abs(n[2]) <= abs(n[3])) v <- c(-n[3], 0, n[1])
    else v <- c(-n[2], n[1], 0)
    v / sqrt(sum(v*v))
  }

  # --- grids & curvilinear mapping
  Ii <- seq(0, 1, length.out = NX + 1)
  Jj <- seq(0, 1, length.out = NY + 1)
  Kk <- seq(0, 1, length.out = NZ + 1)

  n_tot <- length(Ii) * length(Jj) * length(Kk)
  P0 <- matrix(NA_real_, n_tot, 3)  # bases
  V  <- matrix(NA_real_, n_tot, 3)  # vectors F
  idx <- 0L

  for (i in Ii) {
    x <- (1 - i) * a + i * b
    y1 <- H1(x); y2 <- H2(x)
    for (j in Jj) {
      y <- (1 - j) * y1 + j * y2
      g1 <- G1(x, y); g2 <- G2(x, y)
      for (k in Kk) {
        z <- (1 - k) * g1 + k * g2
        idx <- idx + 1L
        P0[idx, ] <- c(x, y, z)
        vec <- F(x, y, z)
        if (!is.numeric(vec) || length(vec)!=3L) stop("F(x,y,z) must return c(Fx,Fy,Fz).", call. = FALSE)
        V[idx, ] <- vec
      }
    }
  }

  # --- bounding box & saturated arrow lengths
  rx <- range(P0[,1], finite = TRUE); ry <- range(P0[,2], finite = TRUE); rz <- range(P0[,3], finite = TRUE)
  span <- max(diff(rx), diff(ry), diff(rz)); if (!is.finite(span) || span <= 0) span <- 1

  vmag <- sqrt(rowSums(V*V))
  Vsat <- V / sqrt(vmag^2 + normalize_bias)
  lsat <- sqrt(rowSums(Vsat*Vsat))
  diru <- Vsat
  nzr  <- lsat > 0
  diru[nzr, ] <- Vsat[nzr, , drop = FALSE] / lsat[nzr]
  diru[!nzr, ] <- 0

  L   <- arrow_scale * span * lsat
  P1  <- P0 + diru * L

  field_points   <- data.frame(x = P0[,1], y = P0[,2], z = P0[,3], mag = vmag)
  field_segments <- data.frame(x0 = P0[,1], y0 = P0[,2], z0 = P0[,3],
                               x1 = P1[,1], y1 = P1[,2], z1 = P1[,3])

  # --- streamline (internal RK4)
  F_eval <- function(x, y, z, t) {
    if (length(formals(F)) >= 4L) F(x, y, z, t) else F(x, y, z)
  }

  if (abs(T) < .Machine$double.eps) {
    traj <- data.frame(t = 0, x = p[1], y = p[2], z = p[3],
                       speed = sqrt(sum(F_eval(p[1],p[2],p[3],0)^2)))
  } else {
    n_steps <- ceiling(abs(T) / step)
    signT   <- if (T >= 0) 1 else -1
    dt_def  <- signT * step
    t <- numeric(n_steps + 1L); t[1] <- 0
    Y <- matrix(NA_real_, n_steps + 1L, 3); Y[1, ] <- p

    for (ii in 1:n_steps) {
      dt <- dt_def
      if ((signT > 0 && (t[ii] + dt) > T) || (signT < 0 && (t[ii] + dt) < T)) dt <- T - t[ii]
      x <- Y[ii, 1]; y <- Y[ii, 2]; z <- Y[ii, 3]; ti <- t[ii]
      k1 <- F_eval(x, y, z, ti)
      k2 <- F_eval(x + 0.5*dt*k1[1], y + 0.5*dt*k1[2], z + 0.5*dt*k1[3], ti + 0.5*dt)
      k3 <- F_eval(x + 0.5*dt*k2[1], y + 0.5*dt*k2[2], z + 0.5*dt*k2[3], ti + 0.5*dt)
      k4 <- F_eval(x + dt*k3[1],    y + dt*k3[2],    z + dt*k3[3],    ti + dt)
      incr <- (k1 + 2*k2 + 2*k3 + k4) / 6
      Y[ii+1, ] <- c(x, y, z) + dt * incr
      t[ii+1]   <- t[ii] + dt
    }

    speed <- vapply(seq_along(t), function(jj) {
      fj <- F_eval(Y[jj,1], Y[jj,2], Y[jj,3], t[jj])
      sqrt(sum(fj * fj))
    }, numeric(1))

    traj <- data.frame(t = t, x = Y[,1], y = Y[,2], z = Y[,3], speed = speed)
  }

  # --- figure
  plt <- plotly::plot_ly()

  # field shafts
  if (arrows %in% c("both","line","cone")) {
    xs <- as.numeric(rbind(field_segments$x0, field_segments$x1, NA))
    ys <- as.numeric(rbind(field_segments$y0, field_segments$y1, NA))
    zs <- as.numeric(rbind(field_segments$z0, field_segments$z1, NA))
    plt <- plotly::add_trace(
      plt, x = xs, y = ys, z = zs,
      type = "scatter3d", mode = "lines",
      line = list(color = normal_color, width = normal_width),
      hoverinfo = "none", showlegend = FALSE
    )
  }

  # field heads
  if (arrows %in% c("both","cone")) {
    if (isTRUE(has_cone)) {
      cone_len <- max(1e-8, arrow_scale * span * arrow_size)
      cs_cone  <- list(list(0, arrow_color), list(1, arrow_color))
      plt <- plotly::add_trace(type = "cone",
        plt,
        x = field_segments$x1, y = field_segments$y1, z = field_segments$z1,
        u = diru[,1], v = diru[,2], w = diru[,3],
        anchor   = "tip",
        sizemode = "absolute",
        sizeref  = cone_len,
        colorscale = cs_cone, showscale = FALSE,
        opacity    = arrow_opacity
      )
    } else {
      head_len <- max(1e-8, arrow_scale * span * arrow_size)
      head_w   <- head_len * 0.5
      hx <- hy <- hz <- numeric(0)
      for (i in seq_len(nrow(field_segments))) {
        tip  <- c(field_segments$x1[i], field_segments$y1[i], field_segments$z1[i])
        n    <- diru[i, ]
        base <- tip - head_len * n
        pperp <- perp_unit(n)
        left  <- base + head_w * pperp
        right <- base - head_w * pperp
        hx <- c(hx, tip[1], left[1], NA, tip[1], right[1], NA)
        hy <- c(hy, tip[2], left[2], NA, tip[2], right[2], NA)
        hz <- c(hz, tip[3], left[3], NA, tip[3], right[3], NA)
      }
      plt <- plotly::add_trace(
        plt, x = hx, y = hy, z = hz,
        type = "scatter3d", mode = "lines",
        line = list(color = arrow_color, width = max(1, normal_width)),
        hoverinfo = "none", showlegend = FALSE
      )
    }
  }

  # streamline
  plt <- plotly::add_trace(
    plt,
    x = traj$x, y = traj$y, z = traj$z,
    type = "scatter3d",
    mode = if (isTRUE(traj_markers)) "lines+markers" else "lines",
    line = list(color = traj_color, width = traj_width),
    marker = list(size = traj_marker_size, color = traj_color),
    showlegend = FALSE, name = "streamline"
  )

  # start/end markers
  plt <- plotly::add_markers(
    plt, x = traj$x[1], y = traj$y[1], z = traj$z[1],
    marker = list(size = traj_marker_size + 2, color = "forestgreen"),
    showlegend = FALSE, hoverinfo = "none"
  )
  plt <- plotly::add_markers(
    plt, x = tail(traj$x,1), y = tail(traj$y,1), z = tail(traj$z,1),
    marker = list(size = traj_marker_size + 2, color = "firebrick"),
    showlegend = FALSE, hoverinfo = "none"
  )

  plt <- plotly::layout(
    plt, title = "Vector field + streamline (combined)",
    scene = scene, paper_bgcolor = bg$paper, plot_bgcolor = bg$plot
  )

  print(plt)
  list(field_points = field_points, field_segments = field_segments, traj = traj, fig = plt)
}
