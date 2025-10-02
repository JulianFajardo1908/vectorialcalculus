test_that("arc_length3d works with a helix", {
  skip_if_not_installed("pracma")

  f_x <- function(t) cos(t)
  f_y <- function(t) sin(t)
  f_z <- function(t) t

  result <- arc_length3d(f_x, f_y, f_z, 0, 2*pi, n_samples = 200)

  expect_true(is.numeric(result))
  expect_length(result, 1)
  expect_gt(result, 0)
})

test_that("arc_length3d works with a straight line", {
  skip_if_not_installed("pracma")

  f_x <- function(t) t
  f_y <- function(t) 0*t
  f_z <- function(t) 0*t

  result <- arc_length3d(f_x, f_y, f_z, 0, 1)

  expect_equal(round(result, 6), 1)
})
