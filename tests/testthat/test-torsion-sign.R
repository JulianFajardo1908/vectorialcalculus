test_that("third derivative sign gives positive torsion for a right-handed helix", {
  result <- curvature_torsion3d(
    X = function(t) cos(t),
    Y = function(t) sin(t),
    Z = function(t) t,
    t0 = 1, h = 1e-3, plot = FALSE
  )
  expect_equal(result$kappa, 1 / 2, tolerance = 1e-2)
  expect_equal(result$tau, 1 / 2, tolerance = 1e-2)
})

test_that("reversing helix handedness reverses torsion", {
  result <- curvature_torsion3d(
    X = function(t) cos(t),
    Y = function(t) sin(t),
    Z = function(t) -t,
    t0 = 1, h = 1e-3, plot = FALSE
  )
  expect_equal(result$kappa, 1 / 2, tolerance = 1e-2)
  expect_equal(result$tau, -1 / 2, tolerance = 1e-2)
})

test_that("unit circle in the plane has zero torsion", {
  result <- curvature_torsion3d(
    X = function(t) cos(t),
    Y = function(t) sin(t),
    Z = function(t) 0 * t,
    t0 = 1, h = 1e-3, plot = FALSE
  )
  expect_equal(result$kappa, 1, tolerance = 1e-2)
  expect_equal(result$tau, 0, tolerance = 1e-2)
})
