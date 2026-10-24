test_that("double integration agrees across both triangular descriptions", {
  # Both Type I and Type II describe x >= 0, y >= 0, x + y <= 1.
  # Analytic integral of x + y over this triangle is 1/3.
  skip_if_not_installed("purrr")
  # Current implementation checks ggplot2 even when plot_domain = FALSE.
  skip_if_not_installed("ggplot2")

  integrand <- function(x, y) x + y
  lower <- function(t) 0
  upper <- function(t) 1 - t

  type1 <- integrate_double_xy(
    f = integrand,
    const_min = 0, const_max = 1,
    limit1 = lower, limit2 = upper,
    region_type = "type1",
    n_outer = 40, n_inner = 40,
    plot_domain = FALSE
  )
  type2 <- integrate_double_xy(
    f = integrand,
    const_min = 0, const_max = 1,
    limit1 = lower, limit2 = upper,
    region_type = "type2",
    n_outer = 40, n_inner = 40,
    plot_domain = FALSE
  )

  expect_named(type1, c("integral_value", "domain_plot"))
  expect_null(type1$domain_plot)
  expect_null(type2$domain_plot)
  expect_equal(type1$integral_value, 1 / 3, tolerance = 1e-9)
  expect_equal(type2$integral_value, 1 / 3, tolerance = 1e-9)
  expect_equal(type1$integral_value, type2$integral_value, tolerance = 1e-10)
})

test_that("double integration rejects odd Simpson subdivision counts", {
  skip_if_not_installed("purrr")
  skip_if_not_installed("ggplot2")
  expect_error(
    integrate_double_xy(
      f = function(x, y) x + y,
      const_min = 0, const_max = 1,
      limit1 = function(x) 0,
      limit2 = function(x) 1,
      n_outer = 3, n_inner = 4,
      plot_domain = FALSE
    ),
    "positive even integers"
  )
})
