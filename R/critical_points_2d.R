#' Critical points in 2D via gradient and Hessian
#'
#' Finds stationary points of \eqn{f(x,y)} over a given domain using grid starts
#' and local refinement, and classifies them by Hessian eigenvalues.
#'
#' @param f Function \code{f(x,y)} returning numeric scalar.
#' @param xlim,ylim Numeric length-2 vectors, domain limits.
#' @param start_n Integer; number of regular starting points per axis.
#' @param n_rand Integer; extra random starts.
#' @param h Numeric step for finite differences.
#' @param tol_grad,tol_merge,tol_eig Numeric tolerances.
#' @param maxit Integer; max iterations per start.
#' @param optim_method Character; e.g. \code{"BFGS"} or \code{"CG"}.
#' @param plot Logical; draw surface and points.
#' @param grid_plot Integer; grid resolution for plotting.
#' @param surface_colorscale Character; Plotly colorscale.
#' @param surface_opacity Numeric in \eqn{[0,1]}.
#' @param cp_size Numeric point size (suggested \eqn{> 0}).
#' @param cp_colors Named vector for point colors.
#' @param cp_size Numeric; point size.
#' @param scene List; Plotly 3D scene options.
#'
#' @return Tibble/data.frame with columns (x, y, f, type, ...).
#' @export
critical_points_2d <- function(
    f,
    xlim, ylim,
    start_n = c(7L, 7L),
    n_rand = 40L,
    h = NULL,
    tol_grad = 1e-6,
    tol_merge = 1e-3,
    tol_eig = 1e-6,
    maxit = 200,
    optim_method = c("BFGS","Nelder-Mead"),
    plot = TRUE,
    grid_plot = c(60L, 60L),
    surface_colorscale = "YlGnBu",
    surface_opacity = 0.85,
    cp_colors = list(minimum = "#2ca02c", maximum = "#d62728", saddle = "#1f77b4", flat = "#ff7f0e"),
    cp_size = 6,
    scene = list(
      aspectmode = "data",
      xaxis = list(title = "x"),
      yaxis = list(title = "y"),
      zaxis = list(title = "f(x,y)")
    )
){
  x <- y <- grad_norm <- NULL
  optim_method <- match.arg(optim_method)
  stopifnot(is.function(f))
  if (length(xlim)!=2 || length(ylim)!=2 || xlim[2] <= xlim[1] || ylim[2] <= ylim[1]) {
    stop("'xlim' and 'ylim' must be c(min,max) with min < max.", call. = FALSE)
  }

  # ---- steps per coordinate
  .steps <- function(x, y) {
    if (is.null(h)) c(1e-4*(1+abs(x)), 1e-4*(1+abs(y)))
    else if (length(h)==1L) rep(as.numeric(h), 2L)
    else as.numeric(h[1:2])
  }
  # central gradient (order 2)
  .grad <- function(x, y) {
    st <- .steps(x,y); hx <- st[1]; hy <- st[2]
    fxp <- f(x+hx, y); fxm <- f(x-hx, y)
    fyp <- f(x, y+hy); fym <- f(x, y-hy)
    c((fxp - fxm)/(2*hx), (fyp - fym)/(2*hy))
  }
  # central Hessian (order 2)
  .hess <- function(x, y) {
    st <- .steps(x,y); hx <- st[1]; hy <- st[2]
    fv  <- f(x,y)
    fxx <- (f(x+hx,y) - 2*fv + f(x-hx,y)) / (hx^2)
    fyy <- (f(x,y+hy) - 2*fv + f(x,y-hy)) / (hy^2)
    fxy <- (f(x+hx,y+hy) - f(x+hx,y-hy) - f(x-hx,y+hy) + f(x-hx,y-hy)) / (4*hx*hy)
    matrix(c(fxx, fxy, fxy, fyy), 2, 2)
  }
  # ||grad||^2 objective
  .gnorm2 <- function(p) {
    g <- .grad(p[1], p[2]); sum(g*g)
  }
  # classification from Hessian eigenvalues
  .classify <- function(H) {
    ev <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
    pos <- ev >  tol_eig
    neg <- ev < -tol_eig
    if (all(pos)) "minimum"
    else if (all(neg)) "maximum"
    else if (any(pos) && any(neg)) "saddle"
    else "flat"
  }
  # merge nearby points
  .merge_pts <- function(P, r) {
    if (nrow(P) <= 1) return(P)
    keep <- rep(TRUE, nrow(P))
    for (i in seq_len(nrow(P))) {
      if (!keep[i]) next
      di <- sqrt((P$x[i] - P$x)^2 + (P$y[i] - P$y)^2)
      dup <- which(di < r)
      dup <- dup[dup > i]
      keep[dup] <- FALSE
    }
    P[keep, , drop = FALSE]
  }
  # clamp to domain
  .clip <- function(p) {
    c(
      min(max(p[1], xlim[1]), xlim[2]),
      min(max(p[2], ylim[1]), ylim[2])
    )
  }
  # safe optim with fallback
  .safe_optim <- function(par0) {
    par0 <- .clip(par0)
    fn   <- function(p) .gnorm2(p)
    ctrl <- list(maxit = maxit)
    res <- try(stats::optim(par0, fn, method = optim_method, control = ctrl), silent = TRUE)
    if (inherits(res, "try-error")) {
      res <- stats::optim(par0, fn, method = "Nelder-Mead", control = ctrl)
    }
    res
  }

  # ---- multi-start (grid + random)
  xs <- seq(xlim[1], xlim[2], length.out = start_n[1])
  ys <- seq(ylim[1], ylim[2], length.out = start_n[2])
  grid_starts <- as.data.frame(expand.grid(x = xs, y = ys))
  rand_starts <- data.frame(
    x = stats::runif(n_rand, min = xlim[1], max = xlim[2]),
    y = stats::runif(n_rand, min = ylim[1], max = ylim[2])
  )
  starts <- rbind(grid_starts, rand_starts)

  # ---- minimize g from each start
  sols <- lapply(seq_len(nrow(starts)), function(i) .safe_optim(as.numeric(starts[i, ])))
  pars <- do.call(rbind, lapply(sols, function(s) s$par))
  vals <- vapply(sols, function(s) s$value, numeric(1))
  df   <- data.frame(x = pars[,1], y = pars[,2], g = vals)

  # filter by grad norm and domain
  df$grad_norm <- sqrt(df$g)
  df <- subset(df, is.finite(grad_norm) & grad_norm <= tol_grad &
                 x >= xlim[1] & x <= xlim[2] & y >= ylim[1] & y <= ylim[2])

  # order and merge duplicates
  if (nrow(df) > 0) {
    df <- df[order(df$grad_norm), , drop = FALSE]
    df <- .merge_pts(df[, c("x","y","grad_norm")], r = tol_merge)
  }

  # evaluate f and classify
  if (nrow(df) > 0) {
    df$z <- mapply(f, df$x, df$y)
    cls <- character(nrow(df))
    for (i in seq_len(nrow(df))) {
      H  <- .hess(df$x[i], df$y[i])
      cls[i] <- .classify(H)
    }
    df$class <- cls
    df <- df[, c("x","y","z","grad_norm","class")]
  } else {
    df <- data.frame(x = numeric(0), y = numeric(0), z = numeric(0),
                     grad_norm = numeric(0), class = character(0))
  }

  # ---- plot
  fig <- NULL
  if (isTRUE(plot)) {
    if (!requireNamespace("plotly", quietly = TRUE)) {
      warning("Package 'plotly' is required for plotting.", call. = FALSE)
    } else {
      gx <- seq(xlim[1], xlim[2], length.out = grid_plot[1])
      gy <- seq(ylim[1], ylim[2], length.out = grid_plot[2])
      # plotly::add_surface expects z[y,x]
      Z <- outer(gy, gx, function(yy, xx) f(xx, yy))

      fig <- plotly::plot_ly() |>
        plotly::add_surface(
          x = gx, y = gy, z = Z,
          colorscale = surface_colorscale,
          opacity = surface_opacity,
          showscale = FALSE
        ) |>
        plotly::layout(
          title = "Surface f(x,y) with critical points",
          scene = scene
        )

      add_cls <- function(p, cls_name, sym) {
        sub <- df[df$class == cls_name, , drop = FALSE]
        if (nrow(sub) == 0) return(p)
        plotly::add_markers(
          p,
          x = sub$x, y = sub$y, z = sub$z,
          type = "scatter3d", mode = "markers",
          marker = list(
            size  = cp_size,
            color = if (!is.null(cp_colors[[cls_name]])) cp_colors[[cls_name]] else "black",
            line  = list(width = 1, color = "black")
          ),
          symbol = sym,
          name = switch(cls_name,
                        minimum = "minimum", maximum = "maximum",
                        saddle  = "saddle",  flat    = "flat / indeterminate"),
          hovertemplate = paste0(
            "<b>", toupper(cls_name), "</b><br>",
            "x=%{x:.4f}<br>y=%{y:.4f}<br>f=%{z:.6f}<br>",
            "||grad||=%{text:.2e}<extra></extra>"
          ),
          text = sub$grad_norm,
          showlegend = TRUE
        )
      }

      fig <- add_cls(fig, "minimum", "circle")
      fig <- add_cls(fig, "maximum", "diamond")
      fig <- add_cls(fig, "saddle",  "cross")
      fig <- add_cls(fig, "flat",    "square")

      print(fig)
    }
  }

  list(critical_points = df, fig = fig)
}
