# Internal S3 proof of concept. Deliberately not exported.
# Existing public functions and their documented return values are untouched.

.vc_new_univariate <- function(fun, domain, study = domain) {
  if (!is.function(fun)) {
    stop("'fun' must be a function.", call. = FALSE)
  }
  if (!inherits(domain, "vc_interval") ||
      !inherits(study, "vc_interval")) {
    stop("'domain' and 'study' must be vc_interval objects.",
         call. = FALSE)
  }

  # Strict interval inclusion, respecting open/closed endpoints.
  left_valid <- study$lower > domain$lower ||
    (study$lower == domain$lower &&
       (!study$closed[1L] || domain$closed[1L]))
  right_valid <- study$upper < domain$upper ||
    (study$upper == domain$upper &&
       (!study$closed[2L] || domain$closed[2L]))
  if (!(left_valid && right_valid)) {
    stop("'study' must be contained in 'domain'.", call. = FALSE)
  }

  structure(
    list(
      fun = fun,
      input_dim = 1L,
      output_dim = 1L,
      variables = "x",
      domain = domain,
      study = study
    ),
    class = c("vc_univariate", "vc_scalar_field", "vc_math")
  )
}

.vc_evaluate <- function(x, at, ...) {
  UseMethod(".vc_evaluate")
}

.vc_evaluate.vc_univariate <- function(x, at, ...) {
  if (!is.numeric(at) || length(at) < 1L ||
      anyNA(at) || any(!is.finite(at))) {
    stop("'at' must be a nonempty finite numeric vector.", call. = FALSE)
  }
  if (any(!.vc_in_interval(at, x$domain))) {
    stop("Evaluation point is outside the mathematical domain.",
         call. = FALSE)
  }

  # Scalar pointwise evaluation accepts non-vectorized user functions.
  ans <- vapply(at, function(z) {
    value <- x$fun(z)
    if (!is.numeric(value) || length(value) != 1L ||
        !is.finite(value)) {
      stop("The function must return a finite numeric scalar.",
           call. = FALSE)
    }
    as.numeric(value)
  }, numeric(1L))

  if (length(at) == 1L) ans[[1L]] else unname(ans)
}
