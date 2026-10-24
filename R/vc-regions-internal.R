# Internal finite box and planar simple-region prototype.
# Do not export: representation and validation semantics are under review.

.vc_new_box <- function(limits, closed = TRUE) {
  if (!is.matrix(limits) || !is.numeric(limits) ||
      nrow(limits) < 1L || ncol(limits) != 2L ||
      anyNA(limits) || any(!is.finite(limits)) ||
      any(limits[, 1L] >= limits[, 2L])) {
    stop("'limits' must be a finite numeric n-by-2 matrix with lower < upper.",
         call. = FALSE)
  }
  if (is.logical(closed) && length(closed) == 1L && !is.na(closed)) {
    closed <- matrix(closed, nrow = nrow(limits), ncol = 2L)
  }
  if (!is.matrix(closed) || !is.logical(closed) ||
      !identical(dim(closed), dim(limits)) || anyNA(closed)) {
    stop("'closed' must be a nonmissing logical n-by-2 matrix or scalar.",
         call. = FALSE)
  }

  structure(
    list(limits = unname(limits),
         closed = unname(closed),
         dimension = as.integer(nrow(limits))),
    class = c("vc_box", "vc_region", "vc_domain")
  )
}

.vc_box_contains <- function(x, box) {
  if (!inherits(box, "vc_box")) {
    stop("'box' must be a vc_box.", call. = FALSE)
  }
  if (!is.numeric(x) || length(x) != box$dimension ||
      anyNA(x) || any(!is.finite(x))) {
    stop("'x' must be a finite numeric point of the box dimension.",
         call. = FALSE)
  }
  left <- ifelse(box$closed[, 1L],
                 x >= box$limits[, 1L], x > box$limits[, 1L])
  right <- ifelse(box$closed[, 2L],
                  x <= box$limits[, 2L], x < box$limits[, 2L])
  all(left & right)
}

.vc_box_within <- function(study, domain) {
  if (!inherits(study, "vc_box") || !inherits(domain, "vc_box") ||
      study$dimension != domain$dimension) {
    stop("Both boxes must have the same dimension.", call. = FALSE)
  }
  left <- study$limits[, 1L] > domain$limits[, 1L] |
    (study$limits[, 1L] == domain$limits[, 1L] &
      (!study$closed[, 1L] | domain$closed[, 1L]))
  right <- study$limits[, 2L] < domain$limits[, 2L] |
    (study$limits[, 2L] == domain$limits[, 2L] &
      (!study$closed[, 2L] | domain$closed[, 2L]))
  all(left & right)
}

.vc_region_contains <- function(region, at, ...) {
  UseMethod(".vc_region_contains")
}

.vc_region_contains.vc_box <- function(region, at, ...) {
  .vc_box_contains(at, region)
}

.vc_check_planar_limits <- function(limits, arg) {
  if (!is.numeric(limits) || length(limits) != 2L ||
      anyNA(limits) || any(!is.finite(limits)) ||
      limits[1L] >= limits[2L]) {
    stop(sprintf("'%s' must be a finite increasing numeric pair.", arg),
         call. = FALSE)
  }
  as.numeric(limits)
}

.vc_new_region_type1 <- function(xlim, lower, upper) {
  xlim <- .vc_check_planar_limits(xlim, "xlim")
  if (!is.function(lower) || !is.function(upper)) {
    stop("'lower' and 'upper' must be functions of x.", call. = FALSE)
  }
  structure(
    list(xlim = xlim, lower = lower, upper = upper,
         dimension = 2L, equivalence = "unverified"),
    class = c("vc_region_type1", "vc_region")
  )
}

.vc_new_region_type2 <- function(ylim, lower, upper) {
  ylim <- .vc_check_planar_limits(ylim, "ylim")
  if (!is.function(lower) || !is.function(upper)) {
    stop("'lower' and 'upper' must be functions of y.", call. = FALSE)
  }
  structure(
    list(ylim = ylim, lower = lower, upper = upper,
         dimension = 2L, equivalence = "unverified"),
    class = c("vc_region_type2", "vc_region")
  )
}

.vc_check_boundary_value <- function(f, x) {
  val <- f(x)
  if (!is.numeric(val) || length(val) != 1L ||
      !is.finite(val)) {
    stop("Each boundary function must return one finite numeric value.",
         call. = FALSE)
  }
  as.numeric(val)
}

.vc_region_contains.vc_region_type1 <- function(region, at, ...) {
  if (!is.numeric(at) || length(at) != 2L ||
      anyNA(at) || any(!is.finite(at))) {
    stop("'at' must be a finite numeric 2D point.", call. = FALSE)
  }
  x <- at[1L]
  if (x < region$xlim[1L] || x > region$xlim[2L]) {
    return(FALSE)
  }
  lo <- .vc_check_boundary_value(region$lower, x)
  hi <- .vc_check_boundary_value(region$upper, x)
  if (lo > hi) {
    stop("Type I region has lower(x) > upper(x) at the sampled x.",
         call. = FALSE)
  }
  at[2L] >= lo && at[2L] <= hi
}

.vc_region_contains.vc_region_type2 <- function(region, at, ...) {
  if (!is.numeric(at) || length(at) != 2L ||
      anyNA(at) || any(!is.finite(at))) {
    stop("'at' must be a finite numeric 2D point.", call. = FALSE)
  }
  y <- at[2L]
  if (y < region$ylim[1L] || y > region$ylim[2L]) {
    return(FALSE)
  }
  lo <- .vc_check_boundary_value(region$lower, y)
  hi <- .vc_check_boundary_value(region$upper, y)
  if (lo > hi) {
    stop("Type II region has lower(y) > upper(y) at the sampled y.",
         call. = FALSE)
  }
  at[1L] >= lo && at[1L] <= hi
}

.vc_new_region_type3 <- function(by_x, by_y) {
  if (!inherits(by_x, "vc_region_type1") ||
      !inherits(by_y, "vc_region_type2")) {
    stop("Type III needs one Type I and one Type II representation.",
         call. = FALSE)
  }
  # Arbitrary R boundary functions cannot generally be proven equivalent.
  structure(
    list(by_x = by_x, by_y = by_y,
         dimension = 2L, equivalence = "unverified"),
    class = c("vc_region_type3", "vc_region")
  )
}

.vc_region_contains.vc_region_type3 <- function(region, at, ...) {
  a <- .vc_region_contains(region$by_x, at)
  b <- .vc_region_contains(region$by_y, at)
  if (!identical(a, b)) {
    stop("Type III representations disagree at the queried point.",
         call. = FALSE)
  }
  a
}
