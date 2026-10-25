# Diagnostic tests only: intentionally expected to reveal known defects in 1.0.5.
# Review test output before making any production changes.

test_that("polar integral of unit disk equals pi", {
  skip_if_not_installed("purrr")
  out <- suppressWarnings(integrate_double_polar(
    function(x, y) 1, 0, 2 * pi,
    function(theta) 0, function(theta) 1,
    n_theta = 2, n_r = 2, plot_domain = FALSE))
  expect_equal(out$integral_value, pi, tolerance = 1e-10)
})

test_that("3D constant-field work includes entire interval", {
  F <- function(x, y, z) c(1, 0, 0)
  a <- line_integral3d_work(F, function(t) c(t, 0, 0),
                            0, 1, plot = FALSE, n_curve = 600, n_field = 2)
  b <- line_integral3d_work(F, function(t) c(1 - t, 0, 0),
                            0, 1, plot = FALSE, n_curve = 600, n_field = 2)
  expect_equal(a$value, 1, tolerance = 1e-8)
  expect_equal(b$value, -1, tolerance = 1e-8)
})

test_that("torsion has the expected orientation", {
  X <- function(t) cos(t); Y <- function(t) sin(t)
  p <- curvature_torsion3d(X, Y, function(t) t, 0, h = 1e-3, plot = FALSE)
  m <- curvature_torsion3d(X, Y, function(t) -t, 0, h = 1e-3, plot = FALSE)
  z <- curvature_torsion3d(X, Y, function(t) 0, 0, h = 1e-3, plot = FALSE)
  expect_equal(p$kappa, 0.5, tolerance = 1e-3)
  expect_equal(p$tau, 0.5, tolerance = 2e-3)
  expect_equal(m$kappa, 0.5, tolerance = 1e-3)
  expect_equal(m$tau, -0.5, tolerance = 2e-3)
  expect_equal(z$kappa, 1, tolerance = 1e-3)
  expect_equal(z$tau, 0, tolerance = 2e-3)
})

test_that("area of unit parametric square ignores unrelated step sizes", {
  skip_if_not_installed("plotly")
  x <- suppressMessages(surface_parametric_area(
    function(u, v) u, function(u, v) v, function(u, v) 0,
    urange = c(0, 1), vrange = c(0, 1),
    nu = 11, nv = 11, h_u = 0.05, h_v = 0.05))
  expect_equal(x$area, 1, tolerance = 1e-8)
})

test_that("bounded critical points never use stale unclipped gradient", {
  r <- critical_points_nd(
    function(x) sum((x - 2)^2), rbind(c(-1, 1), c(-1, 1)),
    start_grid = c(2L, 2L), n_random = 0L,
    max_grid_starts = 4L, maxit = 300, seed = 123)
  expect_equal(nrow(r$critical_points), 0L)
})

test_that("genuine interior critical point remains detectable", {
  r <- critical_points_nd(
    function(x) sum(x^2), rbind(c(-1, 1), c(-1, 1)),
    start_grid = c(3L, 3L), n_random = 0L,
    max_grid_starts = 9L, maxit = 300, seed = 123)
  pts <- r$critical_points
  expect_true(nrow(pts) >= 1L)
  expect_true(any(abs(pts$x1) < 1e-3 & abs(pts$x2) < 1e-3))
})
