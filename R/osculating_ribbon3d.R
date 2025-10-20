#' Osculating ribbon along a 3D curve
#'
#' Constructs and (optionally) plots a narrow ribbon that follows a spatial
#' parametric curve \eqn{\mathbf{r}(t) = (x(t), y(t), z(t))}, based on its Frenet
#' frame \eqn{(\mathbf{T}, \mathbf{N}, \mathbf{B})}. The ribbon extends a small
#' distance along the normal direction, forming a “band” that visualizes
#' curvature and torsion.
#'
#' @param X,Y,Z Functions of \code{t} returning \eqn{x(t)}, \eqn{y(t)}, \eqn{z(t)}.
#' @param a,b Numeric scalars; interval endpoints for \eqn{t \in [a,b]}.
#' @param h Numeric step for finite differences (default \code{1e-4}).
#' @param plot Logical; if \code{TRUE}, draws the ribbon using \pkg{plotly}.
#' @param n_t Integer; number of samples along the base curve.
#' @param n_u Integer; subdivisions across the ribbon width.
#' @param u_max Numeric; half-width of the ribbon (in normal direction units).
#' @param colorscale Character; Plotly colorscale for the ribbon (default \code{"Blues"}).
#' @param opacity Numeric in \eqn{[0,1]}; ribbon opacity.
#' @param show_curve Logical; draw the base curve \eqn{\mathbf{r}(t)}.
#' @param show_centers Logical; draw the centerline joining ribbon midpoints.
#' @param curve_line List; style for the base curve (e.g., \code{list(color="black", width=2)}).
#' @param centers_line List; style for the centerline (e.g., \code{list(color="red", width=2)}).
#' @param show_surface_grid Logical; show contour grid on the ribbon surface.
#' @param surface_grid_color Character; color of the grid lines.
#' @param surface_grid_width Numeric; width of the grid lines.
#' @param show_axis_grid Logical; show axis background gridlines.
#' @param scene List; Plotly 3D scene configuration.
#' @param bg List; background colors (\code{paper}, \code{plot}).
#' @param lighting List; Plotly lighting options for the surface.
#' @param tol Numeric tolerance used for normalizations and zero checks.
#'
#' @return
#' \describe{
#'   \item{\code{data}}{List or tibble with \code{t}, base curve coordinates,
#'     and frame vectors \code{T, N, B}.}
#'   \item{\code{plot}}{Plotly object if \code{plot = TRUE}, else \code{NULL}.}
#' }
#'
#' @examples
#' X <- function(t) cos(t)
#' Y <- function(t) sin(t)
#' Z <- function(t) 0.2 * t
#' # osculating_ribbon3d(X, Y, Z, a = 0, b = 4*pi, plot = TRUE)
#'
#' @export
osculating_ribbon3d <- function(
    X, Y, Z,
    a, b,
    h = 1e-4,
    plot = FALSE,
    n_t = 400,
    n_u = 25,
    u_max = 1,
    colorscale = "Blues",
    opacity = 0.35,
    show_curve = TRUE,
    show_centers = TRUE,
    curve_line   = list(color = "black", width = 2),
    centers_line = list(color = "red",   width = 2),
    show_surface_grid = TRUE,
    surface_grid_color = "rgba(60,80,200,0.25)",
    surface_grid_width = 1,
    show_axis_grid = FALSE,
    scene = list(
      aspectmode = "data",
      xaxis = list(title = "x(t)"),
      yaxis = list(title = "y(t)"),
      zaxis = list(title = "z(t)")
    ),
    bg = list(paper = "white", plot = "white"),
    lighting = list(ambient = 1, diffuse = 0.15, specular = 0, roughness = 1, fresnel = 0),
    tol = 1e-10
) {
  if (b < a) stop("'b' must be >= 'a'.", call. = FALSE)
  if (!(u_max > 0 && u_max <= 1)) stop("'u_max' must satisfy 0 < u_max <= 1.", call. = FALSE)

  # Derivatives and utilities
  d1 <- function(f, t, h) (f(t + h) - f(t - h)) / (2*h)
  d2 <- function(f, t, h) (f(t + h) - 2*f(t) + f(t - h)) / (h*h)
  r  <- function(t) c(X(t), Y(t), Z(t))
  dot   <- function(a,b) sum(a*b)
  norm  <- function(a) sqrt(dot(a,a))
  cross <- function(a,b) c(a[2]*b[3]-a[3]*b[2],
                           a[3]*b[1]-a[1]*b[3],
                           a[1]*b[2]-a[2]*b[1])

  t_seq <- seq(a, b, length.out = n_t)
  u_seq <- seq(0, u_max, length.out = n_u)   # fraction of the radius

  # Preallocate matrices (n_u x n_t)
  Xmat <- matrix(NA_real_, n_u, n_t)
  Ymat <- matrix(NA_real_, n_u, n_t)
  Zmat <- matrix(NA_real_, n_u, n_t)

  # Curve and centers (for optional traces)
  curve_df   <- data.frame(t = t_seq, x = NA_real_, y = NA_real_, z = NA_real_)
  centers_df <- data.frame(t = t_seq, cx = NA_real_, cy = NA_real_, cz = NA_real_)

  for (j in seq_along(t_seq)) {
    t0 <- t_seq[j]
    r0 <- r(t0)
    r1 <- c(d1(X, t0, h), d1(Y, t0, h), d1(Z, t0, h))
    r2 <- c(d2(X, t0, h), d2(Y, t0, h), d2(Z, t0, h))

    n1 <- norm(r1)
    if (n1 < tol) next
    T <- r1 / n1
    c12 <- cross(r1, r2)
    n12 <- norm(c12)
    if (n12 < tol) next

    B <- c12 / n12
    N <- cross(B, T); N <- N / norm(N)
    kappa <- n12 / (n1^3)
    if (kappa <= tol) next

    R <- 1 / kappa
    C <- r0 + N * R

    # Save r(t) and C(t)
    curve_df$x[j] <- r0[1]; curve_df$y[j] <- r0[2]; curve_df$z[j] <- r0[3]
    centers_df$cx[j] <- C[1]; centers_df$cy[j] <- C[2]; centers_df$cz[j] <- C[3]

    # Column j of the ribbon: u in [0, u_max] (fraction of radius)
    Xmat[, j] <- r0[1] + u_seq * R * N[1]
    Ymat[, j] <- r0[2] + u_seq * R * N[2]
    Zmat[, j] <- r0[3] + u_seq * R * N[3]
  }

  # Plotting
  if (isTRUE(plot)) {
    if (!requireNamespace("plotly", quietly = TRUE)) {
      warning("To plot you need 'plotly' installed.", call. = FALSE)
    } else {
      contours_arg <- if (isTRUE(show_surface_grid)) list(
        x = list(show = TRUE, color = surface_grid_color, width = surface_grid_width),
        y = list(show = TRUE, color = surface_grid_color, width = surface_grid_width),
        z = list(show = FALSE)
      ) else NULL

      plt <- plotly::plot_ly() |>
        plotly::add_surface(
          x = Xmat, y = Ymat, z = Zmat,
          colorscale = colorscale,
          showscale  = FALSE,
          opacity    = opacity,
          lighting   = lighting,
          contours   = contours_arg
        )

      if (isTRUE(show_curve) && any(is.finite(curve_df$x))) {
        plt <- plt |>
          plotly::add_trace(
            data = curve_df, x = ~x, y = ~y, z = ~z,
            type = "scatter3d", mode = "lines",
            line = curve_line, hoverinfo = "none", showlegend = FALSE
          )
      }

      if (isTRUE(show_centers) && any(is.finite(centers_df$cx))) {
        plt <- plt |>
          plotly::add_trace(
            x = centers_df$cx, y = centers_df$cy, z = centers_df$cz,
            type = "scatter3d", mode = "lines",
            line = centers_line, hoverinfo = "none", showlegend = FALSE
          )
      }

      scene_final <- scene
      for (ax in c("xaxis","yaxis","zaxis")) {
        if (is.null(scene_final[[ax]])) scene_final[[ax]] <- list()
        scene_final[[ax]]$showgrid <- isTRUE(show_axis_grid)
      }

      plt <- plt |>
        plotly::layout(
          title = sprintf("Osculating ribbon (u_max = %.2f)", u_max),
          scene = scene_final,
          paper_bgcolor = bg$paper,
          plot_bgcolor  = bg$plot
        )

      print(plt)
    }
  }

  list(
    t_seq = t_seq, u_seq = u_seq,
    Xmat = Xmat, Ymat = Ymat, Zmat = Zmat,
    curve   = curve_df,
    centers = centers_df
  )
}
