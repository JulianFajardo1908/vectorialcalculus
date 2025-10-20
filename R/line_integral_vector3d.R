#' 3D line integral of a vector field with an auto-fitted field domain
#'
#' @description
#' Computes the line integral
#' \deqn{\int_{a}^{b} \mathbf{F}(\mathbf{r}(t)) \cdot \mathbf{r}'(t)\, dt}
#' and optionally plots the curve with arrows of \eqn{\mathbf{F}}.
#' The arrow grid is fitted to the curve's bounding box unless you pass explicit ranges.
#'
#' @param F function(x,y,z) -> c(Fx,Fy,Fz)  (or function(x,y,z,t))
#' @param r function(t) -> c(x,y,z)
#' @param a,b parameter limits (b > a)
#' @param plot logical; if TRUE, draw with plotly
#' @param n_curve samples for the curve
#' @param n_field integer >=2; arrows per axis (produces n_field^3 arrows)
#' @param field_ranges optional list(x = c(xmin,xmax), y = c(ymin,ymax), z = c(zmin,zmax)).
#'        If NULL, they are inferred from the curve.
#' @param pad fraction of range added on each side when auto-ranging (e.g. 0.15 = 15%)
#' @param arrows "line" | "cone" | "both" | "none"
#' @param arrow_scale relative arrow length (fraction of max span)
#' @param normalize_bias saturation bias for |F| (bigger -> shorter arrows)
#' @param scene,bg plotly scene & background
#'
#' @return list(value=work, fig=plotly_or_NULL)
#' @export
line_integral3d_work <- function(
    F, r, a, b,
    plot = TRUE,
    n_curve = 600,
    n_field = 7,
    field_ranges = NULL,
    pad = 0.15,
    arrows = c("both","line","cone","none"),
    arrow_scale = 0.10,
    normalize_bias = 1,
    scene = list(aspectmode="data",
                 xaxis=list(title="x"),
                 yaxis=list(title="y"),
                 zaxis=list(title="z")),
    bg = list(paper="white", plot="white")
){
  arrows <- match.arg(arrows)

  if (!is.function(F)) stop("'F' must be function(x,y,z) or (x,y,z,t).")
  if (!is.function(r)) stop("'r' must be function(t)->c(x,y,z).")
  stopifnot(is.numeric(a), is.numeric(b), b > a)
  stopifnot(is.numeric(n_curve), n_curve >= 10, is.numeric(n_field), n_field >= 2)

  # --- helper to call F with/without time ------------------------------
  F_eval <- function(x,y,z,t){
    if (length(formals(F)) >= 4L) F(x,y,z,t) else F(x,y,z)
  }

  # --- curve sampling ---------------------------------------------------
  t_curve <- seq(a, b, length.out = n_curve)
  Rxyz    <- t(vapply(t_curve, r, numeric(3)))
  cx <- Rxyz[,1]; cy <- Rxyz[,2]; cz <- Rxyz[,3]

  # numerical r'(t) (central differences along sampled grid)
  dt  <- (b - a) / (n_curve - 1)
  Rp  <- rbind(Rxyz[2,], Rxyz[3:n_curve,] - Rxyz[1:(n_curve-2),], Rxyz[n_curve-1,]) / (c(dt, rep(2*dt, n_curve-2), dt))
  # integrand F(r(t))·r'(t)
  Fr  <- t(mapply(function(x,y,z,t) F_eval(x,y,z,t), cx, cy, cz, t_curve))
  dot <- rowSums(Fr * Rp)
  # composite Simpson (ensure even subintervals)
  if ((n_curve - 1) %% 2 == 1) {
    t_curve <- t_curve[-length(t_curve)]
    dot     <- dot[-length(dot)]
  }
  hS <- (tail(t_curve,1) - t_curve[1]) / (length(t_curve)-1)
  work <- hS * (dot[1] + dot[length(dot)] +
                  4*sum(dot[seq(2, length(dot)-1, by=2)]) +
                  2*sum(dot[seq(3, length(dot)-2, by=2)])) / 3

  # --- arrow domain (auto from curve bbox unless provided) -------------
  if (is.null(field_ranges)) {
    rx <- range(cx); ry <- range(cy); rz <- range(cz)
    # avoid zero span (flat curves); give a small cushion
    padv <- function(rg) {
      span <- diff(rg); if (span <= 0) span <- 1
      c(rg[1] - pad*span, rg[2] + pad*span)
    }
    xlim <- padv(rx); ylim <- padv(ry); zlim <- padv(rz)
  } else {
    stopifnot(is.list(field_ranges), !is.null(field_ranges$x),
              !is.null(field_ranges$y), !is.null(field_ranges$z))
    xlim <- field_ranges$x; ylim <- field_ranges$y; zlim <- field_ranges$z
  }

  # grid for arrows
  xs <- seq(xlim[1], xlim[2], length.out = n_field)
  ys <- seq(ylim[1], ylim[2], length.out = n_field)
  zs <- seq(zlim[1], zlim[2], length.out = n_field)
  G   <- as.matrix(expand.grid(xs, ys, zs))
  colnames(G) <- c("x","y","z")

  V <- t(mapply(function(x,y,z) F_eval(x,y,z,0), G[,"x"], G[,"y"], G[,"z"]))
  vmag <- sqrt(rowSums(V*V))
  Vsat <- V / sqrt(vmag^2 + normalize_bias)             # saturation
  lsat <- sqrt(rowSums(Vsat*Vsat))
  nz   <- lsat > 0
  diru <- Vsat; diru[nz,] <- Vsat[nz,] / lsat[nz]; diru[!nz,] <- 0
  span <- max(diff(xlim), diff(ylim), diff(zlim)); if (!is.finite(span) || span<=0) span <- 1
  L    <- arrow_scale * span * lsat
  tips <- G + diru * L

  fig <- NULL
  if (isTRUE(plot)) {
    if (!requireNamespace("plotly", quietly = TRUE)) {
      warning("Install 'plotly' to plot.")
    } else {
      plt <- plotly::plot_ly()

      # field body (lines)
      if (arrows %in% c("both","line","cone")) {
        xsL <- as.numeric(rbind(G[,"x"], tips[,"x"], NA))
        ysL <- as.numeric(rbind(G[,"y"], tips[,"y"], NA))
        zsL <- as.numeric(rbind(G[,"z"], tips[,"z"], NA))
        plt <- plotly::add_trace(
          plt, x = xsL, y = ysL, z = zsL,
          type = "scatter3d", mode = "lines",
          line = list(color = "rgba(60,60,60,0.55)", width = 2),
          showlegend = FALSE, hoverinfo = "none"
        )
      }
      # field heads (cones if available)
      if (arrows %in% c("both","cone")) {
        if ("add_cone" %in% getNamespaceExports("plotly")) {
          plt <- plotly::add_trace(type = "cone",
            plt,
            x = tips[,"x"], y = tips[,"y"], z = tips[,"z"],
            u = diru[,1], v = diru[,2], w = diru[,3],
            anchor   = "tip",
            sizemode = "absolute",
            sizeref  = max(1e-8, arrow_scale * span * 0.35),
            colorscale = list(list(0,"#444"), list(1,"#444")),
            showscale  = FALSE, opacity = 0.9
          )
        }
      }

      # curve colored by instantaneous power F·r'
      pwr_col <- grDevices::colorRampPalette(c("#2c7bb6","#abd9e9","#ffffbf","#fdae61","#d7191c"))(length(dot))
      plt <- plotly::add_trace(
        plt, x = cx, y = cy, z = cz,
        type = "scatter3d", mode = "lines+markers",
        line = list(width = 6, color = dot, colorscale = "RdBu", reversescale = TRUE),
        marker = list(size = 2, color = dot, colorscale = "RdBu", reversescale = TRUE),
        showlegend = FALSE
      )
      # start/end markers
      plt <- plotly::add_markers(plt, x = cx[1], y = cy[1], z = cz[1],
                                 marker=list(size=4, color="forestgreen"),
                                 showlegend = FALSE, hoverinfo="none")
      plt <- plotly::add_markers(plt, x = tail(cx,1), y = tail(cy,1), z = tail(cz,1),
                                 marker=list(size=4, color="firebrick"),
                                 showlegend = FALSE, hoverinfo="none")

      ttl <- sprintf("3D line integral (work) \u2248 %.6g", work)
      plt <- plotly::layout(plt, title = ttl, scene = scene,
                            paper_bgcolor = bg$paper, plot_bgcolor = bg$plot)
      fig <- plt
      print(fig)
    }
  }

  list(value = work, fig = fig,
       field_box = list(x=xlim, y=ylim, z=zlim))
}


