test_that("finite boxes preserve dimension and endpoint inclusion", {
  make <- vectorialcalculus:::.vc_new_box
  contains <- vectorialcalculus:::.vc_box_contains
  box <- make(rbind(c(0, 1), c(-2, 2)),
              closed = matrix(c(FALSE, TRUE, TRUE, TRUE), nrow = 2))
  expect_s3_class(box, "vc_box")
  expect_equal(box$dimension, 2L)
  expect_true(contains(c(0.5, 0), box))
  expect_false(contains(c(0, 0), box))
  expect_error(make(matrix(c(0, 0), nrow = 1)), "lower < upper")
})

test_that("scalar fields support vector and coordinate call conventions", {
  make_box <- vectorialcalculus:::.vc_new_box
  make_scalar <- vectorialcalculus:::.vc_new_scalar
  evaluate <- vectorialcalculus:::.vc_evaluate
  box <- make_box(rbind(c(-2, 2), c(-2, 2)))
  f <- make_scalar(function(v) v[1]^2 + v[2]^2,
                   input_dim = 2, domain = box)
  g <- make_scalar(function(x, y) x^2 + y^2,
                   input_dim = 2, domain = box,
                   input = "coordinates")
  expect_s3_class(f, "vc_scalar_field")
  expect_equal(evaluate(f, c(1, 2)), 5)
  expect_equal(evaluate(g, c(1, 2)), 5)
  expect_error(evaluate(f, c(3, 1)), "outside")
  expect_error(evaluate(g, c(1, 2, 3)), "input_dim")
})

test_that("scalar study region respects open domain endpoints", {
  make_box <- vectorialcalculus:::.vc_new_box
  make_scalar <- vectorialcalculus:::.vc_new_scalar
  domain <- make_box(matrix(c(0, 2), nrow = 1),
                     closed = matrix(c(FALSE, TRUE), nrow = 1))
  invalid_study <- make_box(matrix(c(0, 1), nrow = 1))
  expect_error(
    make_scalar(function(v) v[1], 1, domain, invalid_study),
    "contained"
  )
  valid_study <- make_box(matrix(c(0, 1), nrow = 1),
                         closed = matrix(c(FALSE, TRUE), nrow = 1))
  expect_s3_class(
    make_scalar(function(v) v[1], 1, domain, valid_study),
    "vc_scalar_field"
  )
})

test_that("planar Type I and II regions describe the same triangle", {
  type1 <- vectorialcalculus:::.vc_new_region_type1
  type2 <- vectorialcalculus:::.vc_new_region_type2
  type3 <- vectorialcalculus:::.vc_new_region_type3
  contains <- vectorialcalculus:::.vc_region_contains

  a <- type1(c(0, 1), function(x) 0, function(x) 1 - x)
  b <- type2(c(0, 1), function(y) 0, function(y) 1 - y)
  both <- type3(a, b)
  expect_identical(both$equivalence, "unverified")
  expect_s3_class(both, "vc_region_type1")
  expect_true(contains(a, c(0.2, 0.3)))
  expect_true(contains(b, c(0.2, 0.3)))
  expect_true(contains(both, c(0.2, 0.3)))
  expect_false(contains(both, c(0.8, 0.8)))
  expect_error(contains(a, c(-0.1, 0.2)), NA)
})

test_that("Type III uses Type I membership even if supplementary Type II differs", {
  type1 <- vectorialcalculus:::.vc_new_region_type1
  type2 <- vectorialcalculus:::.vc_new_region_type2
  type3 <- vectorialcalculus:::.vc_new_region_type3
  contains <- vectorialcalculus:::.vc_region_contains

  a <- type1(c(0, 1), function(x) 0, function(x) 1 - x)
  b <- type2(c(0, 1), function(y) 0, function(y) 1)
  both <- type3(a, b)
  expect_false(contains(both, c(0.8, 0.8)))
  expect_identical(both$equivalence, "unverified")
})

test_that("inverted planar boundaries are not silently integrated", {
  type1 <- vectorialcalculus:::.vc_new_region_type1
  contains <- vectorialcalculus:::.vc_region_contains
  region <- type1(c(0, 1), function(x) 1, function(x) 0)
  expect_error(contains(region, c(0.5, 0.5)), "lower\\(x\\) > upper\\(x\\)")
})
