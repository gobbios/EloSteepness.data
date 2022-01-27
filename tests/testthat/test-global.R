data("empirical_progress")
data("artificial_progress")
xdata <- rbind(empirical_progress, artificial_progress)
xdata$filename <- paste0(xdata$mname, "@@", xdata$method, ".rds")

x <- list.files(system.file("extdata/removal_experiment", package = "EloSteepness.data"), full.names = TRUE)

r <- sapply(x, unzip, list = TRUE, simplify = FALSE)
r <- do.call("rbind", r)

test_that("removal files are all there", {
  expect_true(all(xdata$filename %in% r$Name))
  expect_true(all(r$Name %in% xdata$filename))
})
