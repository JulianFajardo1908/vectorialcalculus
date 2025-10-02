test_that("divergence_field of identity field is dimension", {
  # F: R^3 -> R^3, recibe vector v = (x, y, z)
  # F(v) = (x, y, z); div F = 3
  F <- function(v) c(v[1], v[2], v[3])

  x0 <- c(1, 2, 3)
  d <- divergence_field(F, x0)

  expect_true(is.numeric(d))
  expect_length(d, 1)
  expect_equal(d, 3, tolerance = 1e-8)
})

test_that("divergence_field of constant field is zero", {
  # F(v) = (1, 2, 3); div F = 0
  F <- function(v) c(1, 2, 3)
  x0 <- c(10, -4, 0.1)
  d <- divergence_field(F, x0)
  expect_equal(round(d, 6), 0)
})
