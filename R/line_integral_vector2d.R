#' Line integral of a vector field along a parametric curve (with optional plot)
#'
#' Computes \eqn{\int_a^b F(r(t)) \cdot r'(t)\, dt} for
#' \eqn{F:\mathbb{R}^n \to \mathbb{R}^n} and \eqn{r:[a,b]\to \mathbb{R}^n}
#' with \eqn{n\in\{2,3\}}. Uses a central finite difference for \eqn{r'(t)}.
#' Optionally returns a plotly visualization coloring the curve by the integrand.
#'
#' @param F Vector field; a function \code{function(x)} returning a numeric
#'   vector of the same length as \code{x} (2 or 3).
#' @param r Parametric curve; a function \code{function(t)} returning a numeric
#'   vector of length 2 or 3.
#' @param a Numeric, left endpoint of the parameter interval.
#' @param b Numeric, right endpoint of the parameter interval (must satisfy \eqn{b>a}).
#' @param method Integration method: \code{"adaptive"} (uses \code{stats::integrate})
#'   or \code{"simpson"} (composite Simpson).
#' @param n_simpson Integer, number of subintervals for Simpson (will be forced even).
#' @param h_t Step size for the finite-difference derivative of \eqn{r}; if \code{NULL},
#'   an automatic step is chosen.
#' @param plot Logical. If \code{TRUE}, draws the curve (2D with \eqn{z=0} or 3D in place)
#'   colored by the integrand.
#' @param n_plot Integer, number of samples for the plotted curve.
#' @param scene List with plotly \code{scene} configuration (used if \code{plot=TRUE}).
#' @param bg List with plotly background configuration (used if \code{plot=TRUE}).
#'
#' @return A \code{list} with components:
#' \describe{
#'   \item{\code{value}}{Numeric scalar, the value of the line integral.}
#'   \item{\code{fig}}{A plotly object if \code{plot=TRUE}, otherwise \code{NULL}.}
#' }
#'
#' @examples
#' \dontrun{
#' P <- function(x) c(x[1], x[2])                # example: F(x,y) = (x, y)
#' r <- function(t) c(cos(t), sin(t))            # unit circle
#' # line_integral_vector2d(F = P, r = r, a = 0, b = 2*pi, method = "simpson")
#' }
#'
#' @export
#' @importFrom stats integrate
line_integral_vector <- function(
    F, r, a, b,
    method = c("adaptive","simpson"),
    n_simpson = 1000,
    h_t = NULL,
    plot = TRUE,
    n_plot = 400,
    scene = list(
      aspectmode = "data",
      xaxis = list(title = "x"),
      yaxis = list(title = "y"),
      zaxis = list(title = "z")
    ),
    bg = list(paper = "white", plot = "white")
) {
  method <- match.arg(method)
  # --- checks
  if (!is.function(F)) stop("'F' must be function(x)-> numeric vector.", call. = FALSE)
  if (!is.function(r)) stop("'r' must be function(t)-> numeric vector.", call. = FALSE)
  if (!is.numeric(a) || !is.numeric(b) || length(a)!=1L || length(b)!=1L || !is.finite(a) || !is.finite(b) || b <= a)
    stop("'a' and 'b' must be scalars with b > a.", call. = FALSE)

  # dimension check once
  r0 <- r((a+b)/2)
  if (!is.numeric(r0) || !(length(r0) %in% c(2,3))) stop("'r(t)' must return length-2 or length-3 numeric.")
  n <- length(r0)

  # --- derivative r'(t) by central differences on the fly
  if (is.null(h_t)) h_t <- (b - a) * 1e-4 + 1e-6
  rprime <- function(t) {
    t1 <- max(a, min(b, t - h_t))
    t2 <- max(a, min(b, t + h_t))
    (r(t2) - r(t1)) / (t2 - t1)
  }

  # --- integrand: F(r(t)) · r'(t)
  integrand <- function(t) {
    rt <- r(t)
    Ft <- F(rt)
    if (!is.numeric(Ft) || length(Ft) != n) stop("F(x) must return vector of same length as x.")
    vt <- rprime(t)
    sum(Ft * vt)
  }
  integrand_vec <- Vectorize(integrand, "t")

  # --- integral
  value <- if (method == "adaptive") {
    stats::integrate(integrand_vec, lower = a, upper = b, rel.tol = 1e-6)$value
  } else {
    m <- as.integer(n_simpson); if (m %% 2L == 1L) m <- m + 1L
    tt <- seq(a, b, length.out = m + 1L)
    yy <- integrand_vec(tt)
    hS <- (b - a) / m
    hS * (yy[1] + yy[length(yy)] + 4 * sum(yy[seq(2, m, by = 2)]) + 2 * sum(yy[seq(3, m-1, by = 2)])) / 3
  }

  # --- plot (optional)
  fig <- NULL
  if (isTRUE(plot)) {
    if (!requireNamespace("plotly", quietly = TRUE)) {
      warning("Plot requested but 'plotly' is not installed.")
    } else {
      ts <- seq(a, b, length.out = n_plot)
      R  <- t(vapply(ts, r, numeric(n)))
      I  <- integrand_vec(ts)

      if (n == 2) {
        # draw curve lifted by integrand as color (z=0 baseline)
        fig <- plotly::plot_ly() |>
          plotly::add_trace(
            x = R[,1], y = R[,2], type = "scatter", mode = "lines",
            line = list(color = "black", width = 3),
            hoverinfo = "none", showlegend = FALSE
          ) |>
          plotly::add_trace(
            x = R[,1], y = R[,2], z = I,
            type = "scatter3d", mode = "lines",
            line = list(width = 6, color = I, colorscale = "RdBu", reversescale = TRUE),
            showlegend = FALSE, name = "integrand"
          ) |>
          plotly::layout(
            title = sprintf("Line integral (value \u2248 %.6g)", value),
            scene = scene, paper_bgcolor = bg$paper, plot_bgcolor = bg$plot
          )
      } else {
        # 3D curve colored by integrand
        fig <- plotly::plot_ly() |>
          plotly::add_trace(
            x = R[,1], y = R[,2], z = R[,3],
            type = "scatter3d", mode = "lines",
            line = list(width = 6, color = I, colorscale = "RdBu", reversescale = TRUE),
            showlegend = FALSE
          ) |>
          plotly::layout(
            title = sprintf("Line integral (value \u2248 %.6g)", value),
            scene = scene, paper_bgcolor = bg$paper, plot_bgcolor = bg$plot
          )
      }
      print(fig)
    }
  }

  list(value = value, fig = fig)
}#' Vector line integral in R^2 with field + trajectory visualization
#'
#' Computes \eqn{W=\int_a^b F(r(t)) \cdot r'(t)\,dt} and draws:
#' - A quiver (arrows) of the vector field \eqn{F(x,y)} on a grid around the path.
#' - The trajectory \eqn{r(t)} colored by the power density \eqn{p(t)=F(r(t))\cdot r'(t)}.
#'
#' @param F function(x,y) -> numeric(2) returning c(Fx,Fy).
#' @param r function(t)   -> numeric(2) returning c(x,y).
#' @param a,b scalars with b > a (parameter limits).
#' @param plot logical; if TRUE, draws an interactive plotly figure.
#' @param n_curve integer; samples of the path for plotting.
#' @param grid_n integer; grid resolution per axis for the field.
#' @param padding numeric; margin added around path bounding box for field grid.
#' @param h numeric or NULL; step for finite differences of r'(t). If NULL, chosen automatically.
#' @param method "adaptive" (stats::integrate) or "simpson".
#' @param n_simpson integer (even); subintervals if method="simpson".
#' @param arrow_scale overall arrow length as fraction of plot span (e.g., 0.08).
#' @param normalize_bias saturation bias in arrow scaling: F / sqrt(|F|^2 + bias). Default 1.
#' @param field_color color (or rgba) for arrow bodies.
#' @param field_width numeric line width for arrows.
#' @param traj_palette color scale for power along the trajectory (Plotly name or vector of colors).
#' @param traj_width line width for trajectory.
#' @param show_markers logical; if TRUE, shows markers along trajectory.
#' @param scene,bg plotly layout helpers.
#'
#' @return list(value = numeric, samples = data.frame(t,x,y,px,py,speed,power), fig = plotly or NULL)
#' @export
line_integral_vector2d <- function(
    F, r, a, b,
    plot = TRUE,
    n_curve = 600,
    grid_n  = 15,
    padding = 0.15,
    h = NULL,
    method = c("adaptive","simpson"),
    n_simpson = 1000,
    arrow_scale = 0.08,
    normalize_bias = 1,
    field_color = "rgba(0,0,0,0.55)",
    field_width = 1.8,
    traj_palette = "RdBu",
    traj_width = 5,
    show_markers = FALSE,
    scene = list(aspectmode="data",
                 xaxis=list(title="x"),
                 yaxis=list(title="y"),
                 zaxis=list(title="z")),
    bg = list(paper="white", plot="white")
){
  method <- match.arg(method)
  if (!is.function(F)) stop("'F' must be function(x,y)->c(Fx,Fy).", call. = FALSE)
  if (!is.function(r)) stop("'r' must be function(t)->c(x,y).", call. = FALSE)
  if (!is.numeric(a) || !is.numeric(b) || length(a)!=1L || length(b)!=1L || !is.finite(a) || !is.finite(b) || b <= a)
    stop("'a' and 'b' must be finite scalars with b > a.", call. = FALSE)

  # --- path samples (for plotting and domain box)
  tt   <- seq(a, b, length.out = n_curve)
  Rxy  <- t(vapply(tt, r, numeric(2)))
  xs   <- Rxy[,1]; ys <- Rxy[,2]

  # --- r'(t) via centered finite differences (adaptive near ends)
  if (is.null(h)) h <- (b - a) * 1e-4 + 1e-6
  rprime <- function(t) {
    t1 <- max(a, min(b, t - h))
    t2 <- max(a, min(b, t + h))
    (r(t2) - r(t1)) / (t2 - t1)
  }

  # --- power density p(t) = F(r(t)) · r'(t)
  Fr   <- t(vapply(seq_along(tt), function(i) F(xs[i], ys[i]), numeric(2)))
  rp   <- t(vapply(tt, rprime, numeric(2)))
  speed <- sqrt(rowSums(rp^2))
  power <- rowSums(Fr * rp)

  # --- integrand and integral value
  integrand <- function(t) {
    p <- r(t); v <- rprime(t); fval <- F(p[1], p[2])
    sum(fval * v)
  }
  integrand_vec <- Vectorize(integrand, "t")
  value <- if (method == "adaptive") {
    stats::integrate(integrand_vec, lower = a, upper = b, rel.tol = 1e-6)$value
  } else {
    n <- as.integer(n_simpson); if (n %% 2L == 1L) n <- n + 1L
    tS <- seq(a, b, length.out = n + 1L)
    yS <- integrand_vec(tS)
    hS <- (b - a) / n
    hS * (yS[1] + yS[length(yS)] + 4 * sum(yS[seq(2, n, by = 2)]) + 2 * sum(yS[seq(3, n-1, by = 2)])) / 3
  }

  # --- field grid around path
  rx <- range(xs); ry <- range(ys)
  pad_x <- padding * diff(rx); pad_y <- padding * diff(ry)
  if (!is.finite(pad_x) || pad_x == 0) pad_x <- 1
  if (!is.finite(pad_y) || pad_y == 0) pad_y <- 1
  xr <- c(rx[1] - pad_x, rx[2] + pad_x)
  yr <- c(ry[1] - pad_y, ry[2] + pad_y)

  gx <- seq(xr[1], xr[2], length.out = grid_n)
  gy <- seq(yr[1], yr[2], length.out = grid_n)

  grid_pts <- expand.grid(x = gx, y = gy)
  Fgrid    <- t(vapply(seq_len(nrow(grid_pts)), function(i) F(grid_pts$x[i], grid_pts$y[i]), numeric(2)))
  mag      <- sqrt(rowSums(Fgrid^2))

  # saturate & scale arrows
  Fsat <- Fgrid / sqrt(mag^2 + normalize_bias)
  mlen <- sqrt(rowSums(Fsat^2)) # in [0,1)
  span <- max(diff(xr), diff(yr)); if (!is.finite(span) || span <= 0) span <- 1
  L    <- arrow_scale * span * mlen
  diru <- Fsat
  nz   <- mlen > 0
  diru[nz, ] <- Fsat[nz, , drop = FALSE] / mlen[nz]
  diru[!nz,] <- 0
  X0 <- grid_pts$x; Y0 <- grid_pts$y
  X1 <- X0 + diru[,1] * L
  Y1 <- Y0 + diru[,2] * L

  # --- build plot
  fig <- NULL
  if (isTRUE(plot)) {
    if (!requireNamespace("plotly", quietly = TRUE)) {
      warning("To plot you need the 'plotly' package installed.", call. = FALSE)
    } else {
      plt <- plotly::plot_ly()

      # field as line segments (with NA-separators)
      xs_ar <- as.numeric(rbind(X0, X1, NA))
      ys_ar <- as.numeric(rbind(Y0, Y1, NA))
      plt <- plt |>
        plotly::add_trace(
          x = xs_ar, y = ys_ar, z = 0,
          type = "scatter3d", mode = "lines",
          line = list(color = field_color, width = field_width),
          hoverinfo = "none", showlegend = FALSE
        )

      # trajectory colored by power
      # map power to colors by supplying 'color' + colorscale
      plt <- plt |>
        plotly::add_trace(
          x = xs, y = ys, z = rep(0, length(xs)),
          type = "scatter3d",
          mode = if (isTRUE(show_markers)) "lines+markers" else "lines",
          line = list(width = traj_width, color = power, colorscale = traj_palette),
          marker = list(size = 3, color = power, colorscale = traj_palette, showscale = TRUE),
          showlegend = FALSE,
          name = "trajectory",
          hovertemplate = paste(
            "x=%{x:.3f}, y=%{y:.3f}<br>",
            "speed=", sprintf("%.3f", mean(speed)), "<br>",
            "power=%{marker.color:.3f}<extra></extra>"
          )
        ) |>
        plotly::layout(
          title = sprintf("Vector line integral: Work \u2248 %.6g", value),
          scene = scene,
          paper_bgcolor = bg$paper,
          plot_bgcolor  = bg$plot
        )

      fig <- plt
      print(fig)
    }
  }

  samples <- data.frame(t = tt, x = xs, y = ys,
                        px = rp[,1], py = rp[,2],
                        speed = speed, power = power)

  list(value = value, samples = samples, fig = fig)
}
