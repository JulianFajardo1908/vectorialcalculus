# Internal scalar field prototype with a finite box as known domain.

.vc_new_scalar <- function(fun, input_dim, domain, study = domain,
                           variables = NULL, input = c("vector", "coordinates")) {
  input <- match.arg(input)
  if (!is.function(fun) || !is.numeric(input_dim) ||
      length(input_dim) != 1L || !is.finite(input_dim) ||
      input_dim < 1 || input_dim != as.integer(input_dim)) {
    stop("'fun' must be a function and 'input_dim' a positive integer.",
         call. = FALSE)
  }
  n <- as.integer(input_dim)
  if (!inherits(domain, "vc_box") || domain$dimension != n) {
    stop("'domain' must be a vc_box with matching dimension.",
         call. = FALSE)
  }
  if (!inherits(study, "vc_box") || study$dimension != n ||
      !.vc_box_within(study, domain)) {
    stop("'study' must be a vc_box contained in 'domain'.",
         call. = FALSE)
  }
  if (is.null(variables)) variables <- paste0("x", seq_len(n))
  if (!is.character(variables) || length(variables) != n ||
      anyNA(variables) || any(!nzchar(variables)) ||
      anyDuplicated(variables)) {
    stop("'variables' must be n distinct nonempty strings.",
         call. = FALSE)
  }
  structure(
    list(fun = fun, input_dim = n, output_dim = 1L,
         variables = variables, domain = domain, study = study,
         input = input),
    class = c("vc_scalar_field", "vc_math")
  )
}

.vc_evaluate.vc_scalar_field <- function(x, at, ...) {
  if (!is.numeric(at) || length(at) != x$input_dim ||
      anyNA(at) || any(!is.finite(at))) {
    stop("'at' must be a finite numeric point of input_dim.",
         call. = FALSE)
  }
  if (!.vc_box_contains(at, x$domain)) {
    stop("Evaluation point is outside the mathematical domain.",
         call. = FALSE)
  }
  value <- if (identical(x$input, "vector")) {
    x$fun(at)
  } else {
    do.call(x$fun, as.list(unname(at)))
  }
  if (!is.numeric(value) || length(value) != 1L ||
      !is.finite(value)) {
    stop("A scalar field must return one finite numeric value.",
         call. = FALSE)
  }
  as.numeric(value)
}
