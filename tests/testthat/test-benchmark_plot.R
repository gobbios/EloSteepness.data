data("removal_global")
# benchmark_plot(removal_global, empirical = TRUE, benchmark = "pij")
x <- benchmark_plot(removal_global, empirical = FALSE, benchmark = "input")
y <- benchmark_plot(removal_global, empirical = TRUE, benchmark = "pij")


test_that("rhos and betas are within range", {
  expect_true(all(x$rho > 0.9))
  expect_true(all(x$beta > 0.5 & x$beta < 1.2))
  expect_true(all(y$rho > 0.8))
  expect_true(all(y$beta > 0.5 & y$beta < 1.2))
})
