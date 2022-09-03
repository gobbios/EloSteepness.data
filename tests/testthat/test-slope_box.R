data("removal_slopes")
x <- slope_box(removal_slopes, empirical = TRUE)
res1 <- all(names(x) %in% levels(removal_slopes$method))

x <- slope_box(removal_slopes, empirical = FALSE)
res2 <- all(names(x) %in% levels(removal_slopes$method))

test_that("summaries for all methods are returned", {
  expect_true(res1)
  expect_true(res2)
})
