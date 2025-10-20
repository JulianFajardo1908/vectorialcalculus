#' 1D Riemann sums (upper, lower, midpoint) with a 2D plot
#'
#' Computes and visualizes 1D Riemann sums for \eqn{y=f(x)} on \eqn{[x_{\min}, x_{\max}]}.
#' It draws the chosen rectangles (upper, lower, and/or midpoint) and optionally
#' overlays the true curve \eqn{y=f(x)} and a baseline \eqn{y=0}.
#'
#' \strong{Notes:}
#' \itemize{
#'   \item \emph{Upper/lower} tiles use the two endpoints of each subinterval:
#'         \code{upper = max(f(x_i), f(x_{i+1}))}, \code{lower = min(...)}.
#'   \item \emph{Midpoint} tiles use \code{f((x_i + x_{i+1})/2)}.
#'   \item Rectangles are anchored to \code{baseline} (default 0) and extend above/below it.
#' }
#'
#' @param f \code{function(x)} returning a numeric vector \eqn{f(x)}.
#' @param xlim Numeric length-2 vector \code{c(xmin, xmax)} with \code{xmax > xmin}.
#' @param n Integer \eqn{\ge 1}: number of subintervals.
#' @param methods Character vector containing any of \code{c("lower","upper","mid")};
#'   controls which rectangles are drawn. All estimates are always returned.
#' @param show_curve Logical; if \code{TRUE} overlays the true curve.
#' @param curve_res Integer; number of points for the curve polyline.
#' @param colors Named list for rectangle colors: \code{list(lower=, upper=, mid=)}.
#' @param alpha Fill opacity for rectangles in \eqn{[0,1]}.
#' @param edge_color Rectangle border color; \code{edge_width} its width.
#' @param curve_color,curve_width Style for the true curve.
#' @param show_baseline Logical; draw a horizontal baseline at \code{baseline}.
#' @param baseline Baseline y-value (default 0); \code{baseline_color}, \code{baseline_width} style it.
#' @param warn_heavy Warn if \code{n} is large (many polygons).
#' @param edge_width Numeric width for bar edges.
#' @param baseline_color Color for the baseline (y=0).
#' @param baseline_width Numeric width for the baseline.

#'
#' @return A list with:
#' \itemize{
#'   \item \code{lower_sum}, \code{upper_sum}, \code{mid_sum}: numeric estimates,
#'   \item \code{dx}: subinterval width,
#'   \item \code{x_breaks}: partition points,
#'   \item \code{figure}: the \pkg{plotly} object (or \code{NULL} if \pkg{plotly} is missing).
#' }
#'
#' @examples
#' \dontshow{if (interactive()) \{}
#' f <- function(x) sin(2*x) + 0.5*x
#' out <- riemann_sum_1d_plot(
#'   f, xlim = c(0, pi), n = 12,
#'   methods = c("lower","mid","upper"),
#'   show_curve = TRUE
#' )
#' out$lower_sum; out$mid_sum; out$upper_sum
#' \dontshow{\}}
#'
#' @export
riemann_sum_1d_plot <- function(
    f,
    xlim,
    n = 12L,
    methods = c("lower","upper","mid"),
    show_curve = TRUE,
    curve_res = 400L,
    colors = list(lower = "#a1d99b", upper = "#fc9272", mid = "#9ecae1"),
    alpha = 0.8,
    edge_color = "black",
    edge_width = 1.2,
    curve_color = "black",
    curve_width = 2,
    show_baseline = TRUE,
    baseline = 0,
    baseline_color = "gray50",
    baseline_width = 1,
    warn_heavy = TRUE
){
  stopifnot(is.function(f))
  if (length(xlim) != 2 || xlim[2] <= xlim[1]) stop("'xlim' must be c(min, max) with max > min.", call. = FALSE)
  if (n < 1) stop("'n' must be >= 1.", call. = FALSE)

  methods <- intersect(tolower(methods), c("lower","upper","mid"))
  if (length(methods) == 0) stop("Choose at least one of methods = c('lower','upper','mid').", call. = FALSE)

  have_plotly <- requireNamespace("plotly", quietly = TRUE)
  if (!have_plotly) warning("Package 'plotly' not found; returning estimates only.", call. = FALSE)
  if (warn_heavy && n > 800) warning("High n: many polygons to draw; consider reducing 'n'.", call. = FALSE)

  # partition
  x_breaks <- seq(xlim[1], xlim[2], length.out = n + 1L)
  dx <- diff(x_breaks)[1]

  # endpoint & midpoint samples
  f_left  <- f(x_breaks[-(n+1)])
  f_right <- f(x_breaks[-1])
  f_mid   <- f(0.5 * (x_breaks[-(n+1)] + x_breaks[-1]))

  # per-interval heights
  h_lower <- pmin(f_left, f_right)
  h_upper <- pmax(f_left, f_right)
  h_mid   <- f_mid

  # sums
  lower_sum <- sum(h_lower) * dx
  upper_sum <- sum(h_upper) * dx
  mid_sum   <- sum(h_mid)   * dx

  # --- Plot
  fig <- NULL
  if (have_plotly) {
    add_alpha <- function(col, a) {
      if (grepl("^rgba\\(", col)) return(col)
      rgb <- grDevices::col2rgb(col)
      sprintf("rgba(%d,%d,%d,%.3f)", rgb[1], rgb[2], rgb[3], max(0, min(1, a)))
    }

    p <- plotly::plot_ly()

    # baseline
    if (isTRUE(show_baseline)) {
      p <- plotly::add_trace(
        p, x = xlim, y = c(baseline, baseline),
        type = "scatter", mode = "lines",
        line = list(color = baseline_color, width = baseline_width),
        hoverinfo = "none", showlegend = FALSE
      )
    }

    # true curve
    if (isTRUE(show_curve)) {
      xs <- seq(xlim[1], xlim[2], length.out = curve_res)
      ys <- f(xs)
      p <- plotly::add_trace(
        p, x = xs, y = ys,
        type = "scatter", mode = "lines",
        line = list(color = curve_color, width = curve_width),
        name = "f(x)"
      )
    }

    # helper to add all rectangles of one method as a single polygons trace
    add_rects <- function(p, heights, method_name, color_hex) {
      xpoly <- ypoly <- numeric(0)
      for (i in seq_len(n)) {
        x0 <- x_breaks[i]; x1 <- x_breaks[i+1]; h <- heights[i]
        # polygon: (x0,baseline) -> (x1,baseline) -> (x1,h) -> (x0,h) -> close
        xpoly <- c(xpoly, x0, x1, x1, x0, NA_real_)
        ypoly <- c(ypoly, baseline, baseline, h,  h,  NA_real_)
      }
      plotly::add_trace(
        p,
        x = xpoly, y = ypoly,
        type = "scatter", mode = "lines",
        fill = "toself", fillcolor = add_alpha(color_hex, alpha),
        line = list(color = edge_color, width = edge_width),
        name = switch(method_name,
                      lower = sprintf("Lower (sum = %.6g)", lower_sum),
                      upper = sprintf("Upper (sum = %.6g)", upper_sum),
                      mid   = sprintf("Midpoint (sum = %.6g)", mid_sum)
        ),
        hoverinfo = "skip",
        inherit = FALSE
      )
    }

    # draw in a logical order (lower, mid, upper)
    if ("lower" %in% methods) p <- add_rects(p, h_lower, "lower", colors$lower %||% "#a1d99b")
    if ("mid"   %in% methods) p <- add_rects(p, h_mid,   "mid",   colors$mid   %||% "#9ecae1")
    if ("upper" %in% methods) p <- add_rects(p, h_upper, "upper", colors$upper %||% "#fc9272")

    p <- plotly::layout(
      p,
      title = "1D Riemann sums: upper, lower, midpoint",
      xaxis = list(title = "x"),
      yaxis = list(title = "y"),
      showlegend = TRUE
    )
    fig <- p
    print(fig)
  }

  list(
    lower_sum = lower_sum,
    upper_sum = upper_sum,
    mid_sum   = mid_sum,
    dx = dx,
    x_breaks = x_breaks,
    figure = fig
  )
}

# small helper for defaulting
`%||%` <- function(a, b) if (!is.null(a)) a else b
