test_that("gradient_scalar matches an analytic polynomial gradient", {
  f <- function(v) v[1]^2 + 3 * v[2]^2
  out <- gradient_scalar(f, x0 = c(2, -1), h = 1e-4, plot = FALSE)

  expect_type(out, "double")
  expect_length(out, 2L)
  expect_equal(out, c(4, -6), tolerance = 1e-8)
})

test_that("gradient_scalar rejects a nonpositive finite-difference step", {
  f <- function(v) sum(v^2)
  expect_error(
    gradient_scalar(f, x0 = c(1, 2), h = 0, plot = FALSE),
    "must be positive"
  )
})

test_that("divergence_field matches the trace of a linear Jacobian", {
  field <- function(v) c(v[1], 2 * v[2], 3 * v[3])
  out <- divergence_field(
    field, x0 = c(0.2, -1, 3), h = 1e-4, plot = FALSE
  )

  expect_type(out, "double")
  expect_length(out, 1L)
  expect_equal(out, 6, tolerance = 1e-8)
})

test_that("curl3d matches the analytic curl of a rotating field", {
  field <- function(x, y, z) c(-y, x, 0 * z)
  out <- curl3d(field, x = 1, y = -2, z = 3, h = 1e-4)

  expect_named(out, c("omega_x", "omega_y", "omega_z"))
  expect_equal(unname(out), c(0, 0, 2), tolerance = 1e-8)
})

test_that("curl3d treats time as a fixed parameter", {
  field <- function(x, y, z, t) c(-y, x + t, z)
  out <- curl3d(field, x = 1, y = 2, z = 3, tval = 4, h = 1e-4)
  expect_equal(unname(out), c(0, 0, 2), tolerance = 1e-8)
})
