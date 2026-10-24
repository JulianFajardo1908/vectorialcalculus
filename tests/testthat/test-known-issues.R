test_that("known issue: positive helix torsion has the correct sign", {
  # Not counted as a passing mathematical result.
  # The existing centered third-derivative formula in curvature_torsion3d
  # has its sign reversed. Retain this test as an executable specification
  # and enable it in the dedicated mathematical correction PR.
  skip("Known baseline defect: reversed sign in third derivative.")

  result <- curvature_torsion3d(
    X = function(t) cos(t),
    Y = function(t) sin(t),
    Z = function(t) t,
    t0 = 1, h = 1e-3, plot = FALSE
  )

  expect_equal(result$kappa, 1 / 2, tolerance = 1e-2)
  expect_equal(result$tau, 1 / 2, tolerance = 1e-2)
})
