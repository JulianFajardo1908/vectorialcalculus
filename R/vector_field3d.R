#' 3D Vector field in a curvilinear prism \eqn{(x,y,z)}
#'
#' Plots the vector field \eqn{\mathbf{F}(x,y,z)} over the volume
#' \deqn{\{(x,y,z):\ x \in [a,b],\ y \in [(1-u)H_1(x)+uH_2(x)],\ z \in [(1-w)G_1(x,y)+wG_2(x,y)]\}}
#' sampled on a regular grid \eqn{(i,j,k)} of \code{(NX+1)\times(NY+1)\times(NZ+1)} points
#' using the linear blends \eqn{x=(1-i)a+i b}, \eqn{y=(1-j)H_1(x)+j H_2(x)},
#' \eqn{z=(1-k)G_1(x,y)+k G_2(x,y)}.
#'
#' Arrows can be drawn as lines only, cones only (\pkg{plotly} \code{add_cone()}), both,
#' or disabled. If your \pkg{plotly} build does not support \code{add_cone()}, a \emph{chevron}
#' head (two “V” line segments) is used as a fallback.
#'
#' Arrow length is scaled with a saturated norm \eqn{\|\mathbf{F}\|}:
#' \deqn{\tilde{\mathbf{F}} = \mathbf{F} / \sqrt{\|\mathbf{F}\|^2 + \texttt{normalize\_bias}}}
#' and then \eqn{L \propto \|\tilde{\mathbf{F}}\|}. Thus small vectors are short
#' and very large vectors saturate.
#'
#' @param F A \code{function(x,y,z)} returning a numeric vector \code{c(Fx,Fy,Fz)}.
#' @param H1,H2 \code{function(x)} lower/upper bounds in \eqn{y}.
#' @param G1,G2 \code{function(x,y)} lower/upper bounds in \eqn{z}.
#' @param a,b Endpoints in \eqn{x} (\code{b > a}).
#' @param NX,NY,NZ Integers \eqn{\ge 1}: density in \eqn{x}, \eqn{y} (u), \eqn{z} (w).
#'   Sampling uses \code{seq(0,1,length.out = N+1)} for each parameter.
#' @param plot Logical. If \code{TRUE}, draw with \pkg{plotly}.
#' @param arrows Arrow mode: \code{"both"}, \code{"line"}, \code{"cone"} or \code{"none"}.
#' @param arrow_scale Global length scale (fraction of the largest span of the bounding box).
#'   \code{0.08} usually looks good.
#' @param normalize_bias Saturation bias (Maxima-style): uses \eqn{\sqrt{\|F\|^2 + \texttt{bias}}}.
#'   Larger values saturate earlier. Default \code{1}.
#' @param normal_color,normal_width Color and width of the arrow \emph{shaft} (line).
#' @param arrow_color Color of the head (cone/chevron).
#' @param arrow_opacity Opacity for cones (if \code{add_cone()} is available).
#' @param arrow_size Relative head size w.r.t. \code{arrow_scale} (default \code{0.35}).
#' @param scene,bg 3D scene options and background colors (\code{list(paper=, plot=)}).
#'
#' @return A list with
#' \itemize{
#'   \item \code{points}: data.frame with positions \code{x,y,z} and magnitude \code{mag},
#'   \item \code{segments}: data.frame with \code{x0,y0,z0,x1,y1,z1} (arrow shafts),
#'   \item \code{fig}: a \pkg{plotly} object if \code{plot=TRUE}, otherwise \code{NULL}.
#' }
#'
#' @examples
#' \dontshow{if (interactive()) \{}
#' # Domain
#' H1 <- function(x) -1 - x
#' H2 <- function(x)  1 - x^2
#' G1 <- function(x,y) y
#' G2 <- function(x,y) y + 1
#'
#' # Swirling upward field
#' F <- function(x,y,z) c(-y, x, 1)
#'
#' vector_field3d(
#'   F, H1=H1, H2=H2, G1=G1, G2=G2, a=-1, b=1, NX=8, NY=6, NZ=6,
#'   plot=TRUE, arrows="both",
#'   arrow_scale = 0.08, normalize_bias = 1,
#'   normal_color = "rgba(0,0,0,0.6)", normal_width = 2,
#'   arrow_color  = "#1d3557", arrow_opacity = 0.95, arrow_size = 0.35
#' )
#' \dontshow{\}}
#'
#' @export
vector_field3d <- function(
    F, H1, H2, G1, G2,
    a, b, NX = 8, NY = 6, NZ = 6,
    plot = TRUE,
    arrows = c("both","line","cone","none"),
    arrow_scale = 0.08,
    normalize_bias = 1,
    normal_color = "black",
    normal_width = 1.5,
    arrow_color  = "black",
    arrow_opacity = 0.9,
    arrow_size    = 0.35,
    scene = list(
      aspectmode = "data",
      xaxis = list(title = "x"),
      yaxis = list(title = "y"),
      zaxis = list(title = "z")
    ),
    bg = list(paper = "white", plot = "white")
) {
  arrows <- match.arg(arrows)

  # ---- Basic validations
  if (!is.function(F))  stop("'F' must be function(x,y,z) -> c(Fx,Fy,Fz).", call. = FALSE)
  if (!is.function(H1) || !is.function(H2)) stop("'H1' and 'H2' must be functions of x.", call. = FALSE)
  if (!is.function(G1) || !is.function(G2)) stop("'G1' and 'G2' must be functions of (x,y).", call. = FALSE)
  if (!is.numeric(a) || !is.numeric(b) || length(a)!=1L || length(b)!=1L || !is.finite(a) || !is.finite(b) || b <= a)
    stop("'a' and 'b' must be scalars with b > a.", call. = FALSE)
  for (nm in c("NX","NY","NZ")) {
    v <- get(nm); if (!is.numeric(v) || length(v)!=1L || !is.finite(v) || v < 1)
      stop(sprintf("'%s' must be an integer >= 1.", nm), call. = FALSE)
  }

  # ---- Color helpers for cones
  is_color <- function(x) {
    if (!is.character(x) || length(x) != 1) return(FALSE)
    if (grepl("^rgba?\\(", x, ignore.case = TRUE)) return(TRUE)
    if (grepl("^#([0-9A-Fa-f]{6}|[0-9A-Fa-f]{8})$", x)) return(TRUE)
    ok <- TRUE; tryCatch(grDevices::col2rgb(x), error = function(...) ok <<- FALSE); ok
  }
  as_rgba <- function(col, alpha = NULL) {
    if (grepl("^rgba?\\(", col, ignore.case = TRUE)) return(col)
    rgb <- grDevices::col2rgb(col) / 255
    a_  <- if (is.null(alpha)) 1 else max(0, min(1, alpha))
    sprintf("rgba(%g,%g,%g,%g)", 255*rgb[1], 255*rgb[2], 255*rgb[3], a_)
  }
  as_colorscale <- function(x, alpha = NULL) {
    if (is_color(x)) return(list(list(0, as_rgba(x, alpha)), list(1, as_rgba(x, alpha))))
    if (is.character(x) && length(x) > 1) {
      cols <- vapply(x, as_rgba, character(1), alpha = alpha)
      pos  <- seq(0, 1, length.out = length(cols))
      return(lapply(seq_along(cols), function(i) list(pos[i], cols[i])))
    }
    x
  }

  # ---- Grids in [0,1] for i,j,k and mapping to the volume
  Ii <- seq(0, 1, length.out = NX + 1)
  Jj <- seq(0, 1, length.out = NY + 1)
  Kk <- seq(0, 1, length.out = NZ + 1)

  n_tot <- length(Ii) * length(Jj) * length(Kk)
  P0 <- matrix(NA_real_, n_tot, 3)  # base points
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

  # ---- Bounding box and length scale
  rx <- range(P0[,1], finite = TRUE); ry <- range(P0[,2], finite = TRUE); rz <- range(P0[,3], finite = TRUE)
  span <- max(diff(rx), diff(ry), diff(rz))
  span <- ifelse(is.finite(span) && span > 0, span, 1)

  # ---- Saturated lengths and heads
  vmag <- sqrt(rowSums(V*V))
  Vsat <- V / sqrt(vmag^2 + normalize_bias)  # saturation
  lsat <- sqrt(rowSums(Vsat*Vsat))           # factor in [0,1)
  diru <- Vsat
  nzr  <- lsat > 0
  diru[nzr, ] <- Vsat[nzr, , drop = FALSE] / lsat[nzr]
  diru[!nzr, ] <- 0

  L   <- arrow_scale * span * lsat
  P1  <- P0 + diru * L

  # ---- Output data (useful if plot=FALSE)
  points_df <- data.frame(x = P0[,1], y = P0[,2], z = P0[,3], mag = vmag)
  seg_df    <- data.frame(x0 = P0[,1], y0 = P0[,2], z0 = P0[,3],
                          x1 = P1[,1], y1 = P1[,2], z1 = P1[,3])

  fig <- NULL
  if (isTRUE(plot)) {
    if (!requireNamespace("plotly", quietly = TRUE)) {
      warning("To plot you need 'plotly' installed.", call. = FALSE)
    } else {
      plt <- plotly::plot_ly()

      # Shafts (lines)
      if (arrows %in% c("both","line")) {
        xs <- as.numeric(rbind(seg_df$x0, seg_df$x1, NA))
        ys <- as.numeric(rbind(seg_df$y0, seg_df$y1, NA))
        zs <- as.numeric(rbind(seg_df$z0, seg_df$z1, NA))
        plt <- plt |>
          plotly::add_trace(
            x = xs, y = ys, z = zs,
            type = "scatter3d", mode = "lines",
            line = list(color = normal_color, width = normal_width),
            hoverinfo = "none", showlegend = FALSE
          )
      }

      # Heads: cones if available; otherwise chevrons
      if (arrows %in% c("both","cone")) {
        has_cone <- "add_cone" %in% getNamespaceExports("plotly")
        if (isTRUE(has_cone)) {
          cone_len <- max(1e-8, arrow_scale * span * arrow_size)
          cs_cone  <- as_colorscale(arrow_color, alpha = arrow_opacity)
          plt <- plt |>
            plotly::add_cone(
              x = P1[,1], y = P1[,2], z = P1[,3],     # at the tip
              u = diru[,1], v = diru[,2], w = diru[,3],
              anchor   = "tip",
              sizemode = "absolute",
              sizeref  = cone_len,
              colorscale = cs_cone,
              showscale  = FALSE,
              opacity    = arrow_opacity
            )
        } else {
          # Fallback: chevron (two short lines) at the tip
          head_len <- max(1e-8, arrow_scale * span * arrow_size)
          head_w   <- head_len * 0.5

          perp_unit <- function(n) {
            if (abs(n[1]) <= abs(n[2]) && abs(n[1]) <= abs(n[3])) v <- c(0, -n[3], n[2])
            else if (abs(n[2]) <= abs(n[1]) && abs(n[2]) <= abs(n[3])) v <- c(-n[3], 0, n[1])
            else v <- c(-n[2], n[1], 0)
            v / sqrt(sum(v*v))
          }

          hx <- hy <- hz <- numeric(0)
          for (i in seq_len(nrow(P0))) {
            if (lsat[i] <= 0) next
            tip <- P1[i, ]; n <- diru[i, ]
            p   <- perp_unit(n)
            base  <- tip - head_len * n
            left  <- base + head_w * p
            right <- base - head_w * p
            hx <- c(hx, tip[1], left[1], NA, tip[1], right[1], NA)
            hy <- c(hy, tip[2], left[2], NA, tip[2], right[2], NA)
            hz <- c(hz, tip[3], left[3], NA, tip[3], right[3], NA)
          }
          plt <- plt |>
            plotly::add_trace(
              x = hx, y = hy, z = hz,
              type = "scatter3d", mode = "lines",
              line = list(color = arrow_color, width = max(1, normal_width)),
              hoverinfo = "none", showlegend = FALSE
            )
        }
      }

      plt <- plt |>
        plotly::layout(
          title = "3D vector field",
          scene = scene,
          paper_bgcolor = bg$paper,
          plot_bgcolor  = bg$plot
        )

      fig <- plt
      print(fig)
    }
  }

  list(points = points_df, segments = seg_df, fig = fig)
}
