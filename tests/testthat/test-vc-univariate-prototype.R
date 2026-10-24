test_that("interval constructor validates bounds and endpoint semantics", {
  new_interval <- vectorialcalculus:::.vc_new_interval
  contains <- vectorialcalculus:::.vc_in_interval
  open_interval <- new_interval(0, 1, closed = c(FALSE, FALSE))

  expect_s3_class(open_interval, "vc_interval")
  expect_identical(contains(c(0, 0.5, 1), open_interval),
                   c(FALSE, TRUE, FALSE))
  expect_error(new_interval(1, 1), "lower < upper")
  expect_error(new_interval(0, 1, closed = c(TRUE, NA)), "without NA")
})

test_that("univariate S3 prototype evaluates only within its domain", {
  interval <- vectorialcalculus:::.vc_new_interval(-2, 2)
  make_function <- vectorialcalculus:::.vc_new_univariate
  evaluate <- vectorialcalculus:::.vc_evaluate

  f <- make_function(function(x) x^2, domain = interval,
                     study = vectorialcalculus:::.vc_new_interval(-1, 1))

  expect_identical(class(f),
                   c("vc_univariate", "vc_scalar_field", "vc_math"))
  expect_equal(evaluate(f, at = 0.5), 0.25)
  expect_equal(evaluate(f, at = c(-2, 0, 2)), c(4, 0, 4))
  expect_error(evaluate(f, at = 2.1), "outside")
})

test_that("study cannot include an excluded domain endpoint", {
  new_interval <- vectorialcalculus:::.vc_new_interval
  make_function <- vectorialcalculus:::.vc_new_univariate

  domain <- new_interval(0, 2, closed = c(FALSE, TRUE))
  expect_error(
    make_function(function(x) log(x), domain,
                  study = new_interval(0, 1, closed = c(TRUE, TRUE))),
    "contained"
  )
  expect_s3_class(
    make_function(function(x) log(x), domain,
                  study = new_interval(0, 1, closed = c(FALSE, TRUE))),
    "vc_univariate"
  )
})

test_that("non-vectorized functions can be evaluated pointwise", {
  new_interval <- vectorialcalculus:::.vc_new_interval
  make_function <- vectorialcalculus:::.vc_new_univariate
  evaluate <- vectorialcalculus:::.vc_evaluate

  f <- make_function(function(x) {
    if (length(x) != 1L) stop("scalar only")
    sqrt(x)
  }, domain = new_interval(0, 4))

  expect_equal(evaluate(f, at = c(0, 1, 4)), c(0, 1, 2))
  expect_error(evaluate(f, at = NA_real_), "finite")
})

test_that("prototype does not attempt plotting or mutate the study interval", {
  new_interval <- vectorialcalculus:::.vc_new_interval
  make_function <- vectorialcalculus:::.vc_new_univariate
  evaluate <- vectorialcalculus:::.vc_evaluate

  domain <- new_interval(0, 10)
  study <- new_interval(2, 3)
  f <- make_function(function(x) x + 1, domain = domain, study = study)

  # A point outside study but within domain is still evaluable.
  expect_equal(evaluate(f, at = 8), 9)
  expect_identical(f$study, study)
})
