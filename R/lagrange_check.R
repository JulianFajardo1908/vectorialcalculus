#' Optimality check with Lagrange multipliers and bordered Hessian
#'
#' Evaluates, at a candidate point \eqn{x^\*}, the Lagrange conditions for a problem
#' with \strong{m} equality constraints \eqn{g_i(x)=0}, builds the \strong{bordered Hessian},
#' and classifies the candidate as \emph{minimum}, \emph{maximum}, or \emph{indeterminate/saddle}
#' using the signs of the \emph{leading principal minors}.
#'
#' Let \eqn{f:\mathbb{R}^n\to\mathbb{R}} and \eqn{g:\mathbb{R}^n\to\mathbb{R}^m}.
#' At \eqn{x^\*}, with Jacobian \eqn{J = \nabla g(x^\*)} (m×n) and
#' \eqn{H_f = \nabla^2 f(x^\*)}, \eqn{H_{g_i} = \nabla^2 g_i(x^\*)}, the Hessian of
#' the Lagrangian is \eqn{H_L = H_f - \sum_{i=1}^m \lambda_i H_{g_i}}, where
#' \eqn{\lambda} satisfies \eqn{\nabla f(x^\*) = J^\top \lambda}. The \strong{bordered Hessian} is
#' \deqn{B = \begin{bmatrix} 0_{m\times m} & J \\ J^\top & H_L \end{bmatrix}.}
#'
#' Criterion (leading principal minors \eqn{\Delta_k = \det B_{(m+k)}}, \eqn{k=1,\dots,n},
#' where \eqn{B_{(r)}} is the \eqn{r\times r} top-left submatrix):
#' \itemize{
#'   \item \strong{Minimum}: \eqn{(-1)^m \Delta_k > 0} for all \eqn{k=1,\dots,n}.
#'   \item \strong{Maximum}: the signs of \eqn{(-1)^m \Delta_k} \strong{alternate} with k:
#'         negative, positive, negative, \dots
#' }
#'
#' Everything is computed with \strong{second-order central finite differences}.
#'
#' @param f Objective, \code{function(x)} returning a scalar. \code{x} is numeric (length n).
#' @param g Equality constraints. Either:
#'   \itemize{
#'     \item \code{function(x)} returning a numeric vector of length m, or
#'     \item a \code{list} of scalar functions \code{function(x)} (one per constraint).
#'   }
#' @param x Candidate point to evaluate (numeric vector of length n).
#' @param h Finite-difference step(s). Scalar, length-n vector, or \code{NULL}
#'   (default \code{1e-4*(1+abs(x))} componentwise).
#' @param tol Tolerance for considering \eqn{g(x)\approx 0}, \eqn{rank(J)}, and minors \eqn{\Delta_k}.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{ok_stationarity}: norm of KKT residual \eqn{\|\nabla f - J^\top\lambda\|}.
#'   \item \code{ok_feasible}: \eqn{\max_i |g_i(x)|}.
#'   \item \code{lambda}: Lagrange multipliers \eqn{\lambda} (length m).
#'   \item \code{J}: Jacobian (m×n).
#'   \item \code{H_f}: Hessian of \eqn{f} (n×n).
#'   \item \code{H_g}: list of Hessians \eqn{H_{g_i}} (each n×n).
#'   \item \code{H_L}: Hessian of the Lagrangian (n×n).
#'   \item \code{B}: bordered Hessian ( (m+n)×(m+n) ).
#'   \item \code{minors}: \code{data.frame} with \eqn{\Delta_k}, its \code{sign} and \code{logAbsDet}.
#'   \item \code{clasificacion}: \code{"minimo"}, \code{"maximo"} or \code{"indeterminado"}.
#'   \item \code{notas}: notes about \eqn{rank(J)}, near-singularity, etc.
#' }
#'
#' @examples
#' ## 1) Minimum with 1 constraint:  f(x,y)=x^2+y^2,  g(x,y)=x+y-1=0  -> (0.5,0.5)
#' f1 <- function(x) x[1]^2 + x[2]^2
#' g1 <- function(x) x[1] + x[2] - 1
#' lagrange_check(f1, g1, x = c(0.5, 0.5))
#'
#' ## 2) Maximum with 1 constraint:  f(x,y)=-(x^2+y^2),  g(x,y)=x+y-1=0  -> (0.5,0.5)
#' f2 <- function(x) -(x[1]^2 + x[2]^2)
#' lagrange_check(f2, g1, x = c(0.5, 0.5))
#'
#' ## 3) Two constraints in R^3 (min-norm with two planes)
#' f3 <- function(x) sum(x^2)
#' g3 <- list(function(x) x[1] + x[2] + x[3] - 1,
#'            function(x) x[1] - x[3])
#' # Candidate solution: x1=x3, 2*x1 + x2 = 1 -> x=(1/3,1/3,1/3)
#' lagrange_check(f3, g3, x = c(1,1,1)/3)
#'
#' @export
lagrange_check <- function(f, g, x, h = NULL, tol = 1e-6) {
  # ---------- internal helpers (not exported)
  .as_list_g <- function(g) {
    if (is.null(g)) return(list())
    if (is.list(g)) return(g)
    if (is.function(g)) {
      out <- g(x)
      if (!is.numeric(out)) stop("'g(x)' must be numeric.", call. = FALSE)
      m <- length(out)
      return(lapply(seq_len(m), function(i) function(xx) g(xx)[i]))
    }
    stop("'g' must be function(x) or a list of function(x).", call. = FALSE)
  }
  .num_grad <- function(fun, x, h) {
    n <- length(x)
    if (is.null(h)) hi <- 1e-4 * (1 + abs(x)) else hi <- if (length(h)==1) rep(h,n) else as.numeric(h)
    g <- numeric(n)
    for (i in seq_len(n)) {
      ei <- rep(0, n); ei[i] <- hi[i]
      g[i] <- (fun(x + ei) - fun(x - ei)) / (2 * hi[i])
    }
    g
  }
  .num_hess <- function(fun, x, h) {
    n <- length(x)
    if (is.null(h)) hi <- 1e-4 * (1 + abs(x)) else hi <- if (length(h)==1) rep(h,n) else as.numeric(h)
    H <- matrix(0, n, n)
    fx <- fun(x)
    # diagonals
    for (i in seq_len(n)) {
      ei <- rep(0,n); ei[i] <- hi[i]
      fpp <- fun(x + ei); fmm <- fun(x - ei)
      H[i,i] <- (fpp - 2*fx + fmm) / (hi[i]^2)
    }
    # off-diagonals
    if (n >= 2) {
      for (i in 1:(n-1)) for (j in (i+1):n) {
        ei <- rep(0,n); ej <- rep(0,n)
        ei[i] <- hi[i]; ej[j] <- hi[j]
        fpp <- fun(x + ei + ej)
        fpm <- fun(x + ei - ej)
        fmp <- fun(x - ei + ej)
        fmm <- fun(x - ei - ej)
        val <- (fpp - fpm - fmp + fmm) / (4 * hi[i] * hi[j])
        H[i,j] <- H[j,i] <- val
      }
    }
    H
  }
  .det_sign_logabs <- function(M) {
    # LU-based sign and log|det| (more stable than raw det)
    d <- suppressWarnings(determinant(M, logarithm = TRUE))
    sign <- as.numeric(d$sign)      # 0 if singular
    logabs <- as.numeric(d$modulus) # log |det|
    list(sign = sign, logabs = logabs)
  }
  .classify <- function(m, deltas, tol) {
    # Classification from signs of (-1)^m * Δ_k
    s <- (-1)^m * deltas
    if (any(abs(deltas) < tol)) return("indeterminado")
    if (all(s > 0)) return("minimo")
    alt <- all( s[seq(1, length(s), 2)] < 0 ) && all( s[seq(2, length(s), 2)] > 0 )
    if (alt) return("maximo")
    "indeterminado"
  }

  # ---------- validations
  if (!is.function(f)) stop("'f' must be function(x)->scalar.", call. = FALSE)
  if (!is.numeric(x) || any(!is.finite(x))) stop("'x' must be a finite numeric vector.", call. = FALSE)
  x <- as.numeric(x); n <- length(x)

  glist <- .as_list_g(g); m <- length(glist)

  # ---------- gradient of f, Jacobian J, Hessians
  gvals <- if (m>0) vapply(glist, function(gi) gi(x), numeric(1)) else numeric(0)
  grad_f <- .num_grad(f, x, h)
  H_f    <- .num_hess(f, x, h)

  if (m > 0) {
    J <- matrix(NA_real_, nrow = m, ncol = n)
    H_g <- vector("list", m)
    for (i in seq_len(m)) {
      gi <- glist[[i]]
      J[i, ] <- .num_grad(gi, x, h)
      H_g[[i]] <- .num_hess(gi, x, h)
    }
  } else {
    J <- matrix(0, 0, n)
    H_g <- list()
  }

  # ---------- lambda from J J^T lambda = J grad_f (least squares)
  notas <- character()
  if (m > 0) {
    sv <- svd(J)
    rJ <- sum(sv$d > tol * max(1, sv$d[1]))
    if (rJ < m) notas <- c(notas, sprintf("Warning: rank(J)=%d < m=%d; regularity fails.", rJ, m))
    rhs <- as.numeric(J %*% grad_f)
    JJt <- J %*% t(J)
    if (abs(det(JJt)) < tol) {
      notas <- c(notas, "JJ^T nearly singular; using pseudoinverse (SVD).")
      sv2 <- svd(JJt)
      pinv <- sv2$v %*% (t(sv2$u) * ifelse(sv2$d > tol, 1/sv2$d, 0))
      lambda <- as.numeric(pinv %*% rhs)
    } else {
      lambda <- as.numeric(solve(JJt, rhs))
    }
  } else {
    lambda <- numeric(0)
  }

  # ---------- Hessian of the Lagrangian
  H_L <- H_f
  if (m > 0) for (i in seq_len(m)) H_L <- H_L - lambda[i] * H_g[[i]]

  # ---------- Bordered Hessian B = [[0, J], [J^T, H_L]]
  if (m > 0) {
    B <- rbind(
      cbind(matrix(0, m, m), J),
      cbind(t(J), H_L)
    )
  } else {
    B <- H_L
    notas <- c(notas, "No constraints: using H_L (Hessian of f) for the standard test.")
  }

  # ---------- Leading principal minors Δ_k, k=1..n (orders m+k)
  deltas <- numeric(n); signs <- integer(n); logs <- numeric(n)
  for (k in 1:n) {
    idx <- 1:(m + k)
    d <- .det_sign_logabs(B[idx, idx, drop = FALSE])
    signs[k] <- d$sign
    logs[k]  <- d$logabs
    deltas[k] <- if (d$sign == 0) 0 else d$sign * exp(d$logabs)
  }

  clas <- .classify(m, deltas, tol)

  list(
    ok_stationarity = if (m>0) sqrt(sum((grad_f - t(J) %*% lambda)^2)) else sqrt(sum(grad_f^2)),
    ok_feasible = if (m>0) max(abs(gvals)) else 0,
    lambda = lambda,
    J = J,
    H_f = H_f,
    H_g = H_g,
    H_L = H_L,
    B = B,
    minors = data.frame(
      k = 1:n,
      order = m + (1:n),
      Delta_k = deltas,
      sign = signs,
      logAbsDet = logs
    ),
    clasificacion = clas,
    notas = notas
  )
}
