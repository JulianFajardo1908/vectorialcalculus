#' Critical points of an nD scalar field (no plot, n >= 3)
#'
#' Finds critical points of a scalar field \eqn{f:\mathbb{R}^n \to \mathbb{R}} by
#' numerically solving \eqn{\nabla f(x)=0} via multi-start minimization of
#' \eqn{g(x)=\|\nabla f(x)\|^2}. Candidates closer than \code{tol_merge} are merged.
#' Each accepted point is classified using the Hessian eigenvalues.
#'
#' Central finite differences (order 2) are used for gradients and the Hessian.
#'
#' @param f Function \code{function(x)} that returns a numeric scalar; \code{x} is a numeric vector (length \eqn{n}).
#' @param bounds Domain bounds. Either:
#'   \itemize{
#'     \item an \eqn{n \times 2} numeric matrix (or data frame) with columns \code{lower, upper}, or
#'     \item a list of length \eqn{n}, each element \code{c(lower, upper)}.
#'   }
#' @param start_grid Integer vector of length \eqn{n} with the number of grid points per dimension
#'   for deterministic starts (default \code{rep(5, n)} internally). The total grid size is capped by \code{max_grid_starts}.
#' @param n_random Number of additional random starts sampled uniformly in the hyper-rectangle.
#' @param max_grid_starts Maximum number of grid starts to actually instantiate (to avoid explosion).
#'   If the full grid exceeds this, a random subset of that size is used. Default \code{2000}.
#' @param h Step(s) for finite differences. \code{NULL} (auto per component: \code{1e-4 * (1+abs(x_i))}),
#'   a scalar, or a numeric vector of length \eqn{n}.
#' @param tol_grad Threshold on \eqn{\|\nabla f(x)\|} to accept a critical point (default \code{1e-6}).
#' @param tol_merge Euclidean distance to merge nearby candidates (default \code{1e-3}).
#' @param tol_eig Eigenvalue tolerance for classification (positive/negative vs. near-zero).
#' @param maxit Max iterations per \code{stats::optim()} call.
#' @param optim_method Primary method for \code{stats::optim()} (\code{"BFGS"} or \code{"Nelder-Mead"}).
#'   If it errors, a fallback to \code{"Nelder-Mead"} is used.
#' @param seed Optional integer seed for reproducibility of random starts.
#' @param store_hessian Logical; if \code{TRUE}, returns the Hessian matrices for each critical point.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{critical_points}: a data.frame with columns \code{x1..xn}, \code{f}, \code{grad_norm}, \code{class}.
#'   \item \code{eigvals}: a list column of numeric eigenvalue vectors (one per critical point).
#'   \item \code{hessians}: (optional) a list column of Hessian matrices if \code{store_hessian=TRUE}.
#'   \item \code{starts_info}: list with the actual number of grid and random starts used.
#' }
#'
#' @examples
#' \dontshow{if (interactive()) \{}
#' # Example 1: unique minimum at (1,1,1)
#' f1 <- function(x) sum((x - 1)^2)
#' B  <- rbind(c(-2, 3), c(-2, 3), c(-2, 3))   # 3D bounds
#' res1 <- critical_points_nd(f1, bounds = B, start_grid = c(5,5,5), n_random = 50, seed = 1)
#' res1$critical_points
#'
#' # Example 2: saddle at the origin in 3D
#' f2 <- function(x) x[1]^2 + x[2]^2 - x[3]^2
#' B2 <- rbind(c(-1,1), c(-1,1), c(-1,1))
#' res2 <- critical_points_nd(f2, bounds = B2, start_grid = c(5,5,5), n_random = 30, seed = 123)
#' res2$critical_points
#'
#' # Example 3 (4D): multiple critical points
#' f3 <- function(x) sum(x^4 - 2*x^2)  # per-dim extrema at {-1, 0, 1}
#' B3 <- do.call(rbind, replicate(4, c(-2,2), simplify = FALSE))
#' res3 <- critical_points_nd(f3, bounds = B3, start_grid = rep(4,4), n_random = 200, seed = 42)
#' head(res3$critical_points)
#' \dontshow{\}}
#'
#' @export
critical_points_nd <- function(
    f,
    bounds,
    start_grid = NULL,
    n_random = 50L,
    max_grid_starts = 2000L,
    h = NULL,
    tol_grad = 1e-6,
    tol_merge = 1e-3,
    tol_eig = 1e-6,
    maxit = 200,
    optim_method = c("BFGS","Nelder-Mead"),
    seed = NULL,
    store_hessian = FALSE
){
  stopifnot(is.function(f))
  optim_method <- match.arg(optim_method)

  # --- normalize bounds to n x 2 matrix
  if (is.list(bounds) && !is.data.frame(bounds)) {
    B <- do.call(rbind, bounds)
  } else {
    B <- as.matrix(bounds)
  }
  if (!is.numeric(B) || ncol(B) != 2L) stop("'bounds' must be n x 2 (lower, upper).", call. = FALSE)
  if (any(!is.finite(B))) stop("'bounds' must be finite.", call. = FALSE)
  if (any(B[,2] <= B[,1])) stop("Each row of 'bounds' must satisfy upper > lower.", call. = FALSE)
  n <- nrow(B)

  # --- defaults for grid resolution
  if (is.null(start_grid)) start_grid <- rep(5L, n)
  if (length(start_grid) != n) stop("'start_grid' must have length n.", call. = FALSE)
  start_grid <- as.integer(start_grid); if (any(start_grid < 2L)) stop("'start_grid' entries must be >= 2.", call. = FALSE)

  # --- reproducibility
  if (!is.null(seed)) set.seed(as.integer(seed))

  # ---- finite-diff steps
  .steps <- function(x) {
    if (is.null(h)) {
      1e-4 * (1 + abs(x))
    } else if (length(h) == 1L) {
      rep(as.numeric(h), n)
    } else {
      hx <- as.numeric(h)
      if (length(hx) != n) stop("'h' must be length n if not scalar/NULL.", call. = FALSE)
      hx
    }
  }

  # ---- central gradient (order 2)
  .grad <- function(x) {
    hx <- .steps(x)
    g  <- numeric(n)
    for (i in seq_len(n)) {
      ei <- numeric(n); ei[i] <- hx[i]
      g[i] <- (f(x + ei) - f(x - ei)) / (2 * hx[i])
    }
    g
  }

  # ---- central Hessian (order 2)
  .hess <- function(x) {
    hx <- .steps(x)
    H  <- matrix(0, n, n)
    fx <- f(x)
    # diagonal
    for (i in seq_len(n)) {
      ei <- numeric(n); ei[i] <- hx[i]
      fpp <- f(x + ei); fmm <- f(x - ei)
      H[i,i] <- (fpp - 2*fx + fmm) / (hx[i]^2)
    }
    # off-diagonal
    if (n >= 2) {
      for (i in 1:(n-1)) for (j in (i+1):n) {
        ei <- numeric(n); ej <- numeric(n)
        ei[i] <- hx[i]; ej[j] <- hx[j]
        fpp <- f(x + ei + ej)
        fpm <- f(x + ei - ej)
        fmp <- f(x - ei + ej)
        fmm <- f(x - ei - ej)
        val <- (fpp - fpm - fmp + fmm) / (4 * hx[i] * hx[j])
        H[i,j] <- H[j,i] <- val
      }
    }
    H
  }

  # ---- objective g(x) = ||grad f(x)||^2
  .gnorm2 <- function(x) {
    g <- .grad(x); sum(g*g)
  }

  # ---- clamp to bounds
  .clip <- function(x) pmin(pmax(x, B[,1]), B[,2])

  # ---- safe optim with fallback
  .safe_optim <- function(x0) {
    x0 <- .clip(x0)
    ctrl <- list(maxit = maxit)
    fn   <- function(z) .gnorm2(z)
    res  <- try(stats::optim(x0, fn, method = optim_method, control = ctrl), silent = TRUE)
    if (inherits(res, "try-error")) {
      res <- stats::optim(x0, fn, method = "Nelder-Mead", control = ctrl)
    }
    res$par <- .clip(res$par)  # ensure inside bounds
    res
  }

  # ---- generate starts (grid + random)
  make_grid <- function() {
    seqs <- lapply(seq_len(n), function(i) seq(B[i,1], B[i,2], length.out = start_grid[i]))
    sizes <- vapply(seqs, length, integer(1))
    total <- prod(sizes)
    if (total <= max_grid_starts) {
      as.matrix(expand.grid(rev(seqs)))[, n:1, drop = FALSE]  # keep x1..xn order
    } else {
      # sample 'max_grid_starts' random grid tuples
      idx_mat <- sapply(sizes, function(m) sample.int(m, size = max_grid_starts, replace = TRUE))
      out <- matrix(NA_real_, nrow = max_grid_starts, ncol = n)
      for (i in seq_len(n)) out[, i] <- seqs[[i]][ idx_mat[, i] ]
      out
    }
  }

  grid_starts <- make_grid()
  rand_starts <- matrix(stats::runif(n_random * n), nrow = n_random, ncol = n)
  # scale random to bounds
  for (i in seq_len(n)) rand_starts[, i] <- B[i,1] + rand_starts[, i] * (B[i,2] - B[i,1])

  starts <- rbind(grid_starts, rand_starts)
  n_grid_used <- nrow(grid_starts)

  # ---- optimize from each start
  sols <- lapply(seq_len(nrow(starts)), function(i) .safe_optim(starts[i, ]))
  pars <- do.call(rbind, lapply(sols, `[[`, "par"))
  vals <- vapply(sols, `[[`, numeric(1), "value")

  df <- as.data.frame(pars)
  names(df) <- paste0("x", seq_len(n))
  df$g <- vals
  df$grad_norm <- sqrt(df$g)

  # filter by grad_norm and bounds (redundant but safe)
  inside <- apply(df[ , paste0("x", seq_len(n)), drop = FALSE], 1L,
                  function(xx) all(xx >= B[,1] & xx <= B[,2]))
  df <- df[ is.finite(df$grad_norm) & df$grad_norm <= tol_grad & inside, , drop = FALSE ]

  # order and merge nearby candidates
  merge_points <- function(D, r) {
    if (nrow(D) <= 1) return(D)
    coord <- as.matrix(D[ , paste0("x", seq_len(n)), drop = FALSE])
    keep  <- rep(TRUE, nrow(D))
    for (i in seq_len(nrow(D))) {
      if (!keep[i]) next
      di <- sqrt(rowSums( (coord[rep(i, nrow(D)), , drop = FALSE] - coord)^2 ))
      dup <- which(di < r); dup <- dup[dup > i]
      keep[dup] <- FALSE
    }
    D[keep, , drop = FALSE]
  }

  if (nrow(df) > 0) {
    df <- df[order(df$grad_norm), , drop = FALSE]
    df <- merge_points(df, tol_merge)
  }

  # evaluate f and classify by Hessian eigenvalues
  eig_list <- list(); H_list <- NULL
  if (store_hessian) H_list <- vector("list", length = nrow(df))

  if (nrow(df) > 0) {
    df$f <- apply(df[ , paste0("x", seq_len(n)), drop = FALSE], 1L, function(xx) f(as.numeric(xx)))
    classes <- character(nrow(df))
    for (i in seq_len(nrow(df))) {
      xi <- as.numeric(df[i, paste0("x", seq_len(n)), drop = TRUE])
      H  <- .hess(xi)
      ev <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
      pos <- ev >  tol_eig
      neg <- ev < -tol_eig
      cls <- if (all(pos)) "minimum" else if (all(neg)) "maximum" else if (any(pos) && any(neg)) "saddle" else "flat"
      classes[i] <- cls
      eig_list[[i]] <- ev
      if (store_hessian) H_list[[i]] <- H
    }
    df$class <- classes
    # reorder cols
    df <- df[ , c(paste0("x", seq_len(n)), "f", "grad_norm", "class"), drop = FALSE ]
  } else {
    df <- as.data.frame(matrix(numeric(0), nrow = 0, ncol = n))
    names(df) <- paste0("x", seq_len(n))
    df$f <- df$grad_norm <- numeric(0)
    df$class <- character(0)
  }

  out <- list(
    critical_points = df,
    eigvals = eig_list,
    starts_info = list(grid_used = n_grid_used, random_used = n_random)
  )
  if (store_hessian) out$hessians <- H_list
  out
}
