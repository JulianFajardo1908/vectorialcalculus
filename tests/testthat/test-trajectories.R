test_that("curve_sample3d preserves the sampled parametrization", {
  curve <- curve_sample3d(
    X = function(t) t,
    Y = function(t) t^2,
    Z = function(t) 0 * t,
    a = 0, b = 1, n_samples = 5
  )

  expect_s3_class(curve, "tbl_df")
  expect_named(curve, c("t", "x", "y", "z"))
  expect_equal(nrow(curve), 5L)
  expect_equal(curve$t, seq(0, 1, length.out = 5L))
  expect_equal(curve$x, curve$t)
  expect_equal(curve$y, curve$t^2)
  expect_equal(curve$z, rep(0, 5L))
})

test_that("arc_length3d calculates unit-circle length without plotting", {
  length_value <- arc_length3d(
    X = function(t) cos(t),
    Y = function(t) sin(t),
    Z = function(t) 0 * t,
    a = 0, b = 2 * pi,
    h = 1e-5,
    method_int = "integrate",
    plot = FALSE
  )

  expect_type(length_value, "double")
  expect_length(length_value, 1L)
  expect_equal(length_value, 2 * pi, tolerance = 1e-7)
})
