test_that("gradient_scalar of simple function works (vector input)", {
  # f: R^2 -> R, recibe vector v = (x, y)
  f <- function(v) v[1]^2 + v[2]^2

  # punto
  x0 <- c(1, 1)

  g <- gradient_scalar(f, x0)

  expect_true(is.numeric(g))
  expect_length(g, 2)
  # grad f = (2x, 2y) en (1,1) -> (2,2)
  expect_equal(round(g[1], 3), 2, tolerance = 1e-8)
  expect_equal(round(g[2], 3), 2, tolerance = 1e-8)
})

test_that("gradient_scalar returns finite values", {
  f <- function(v) exp(v[1] + v[2])
  x0 <- c(0.5, -0.5)
  g <- gradient_scalar(f, x0)
  expect_true(all(is.finite(g)))
})
