# Internal prototype of a one-dimensional mathematical interval.
# Do not export until the public S3 contract is approved.

.vc_new_interval <- function(lower, upper, closed = c(TRUE, TRUE)) {
  if (!is.numeric(lower) || length(lower) != 1L ||
      is.na(lower) || !is.finite(lower) ||
      !is.numeric(upper) || length(upper) != 1L ||
      is.na(upper) || !is.finite(upper) || lower >= upper) {
    stop("'lower' and 'upper' must be finite scalars with lower < upper.",
         call. = FALSE)
  }
  if (!is.logical(closed) || length(closed) != 2L ||
      anyNA(closed)) {
    stop("'closed' must be a logical vector of length two without NA.",
         call. = FALSE)
  }

  structure(
    list(lower = as.numeric(lower),
         upper = as.numeric(upper),
         closed = unname(closed)),
    class = "vc_interval"
  )
}

.vc_in_interval <- function(x, interval) {
  if (!inherits(interval, "vc_interval")) {
    stop("'interval' must be a vc_interval.", call. = FALSE)
  }
  if (!is.numeric(x) || anyNA(x) || any(!is.finite(x))) {
    stop("'x' must contain finite numeric values.", call. = FALSE)
  }
  left <- if (interval$closed[1L]) {
    x >= interval$lower
  } else {
    x > interval$lower
  }
  right <- if (interval$closed[2L]) {
    x <= interval$upper
  } else {
    x < interval$upper
  }
  left & right
}
