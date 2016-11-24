context("var.wt")

test_that("var.wt matches the covariance matrix from stats::var", {
  expect_equal(stats::var(x = mtcars[, 1:4]), var.wt(x = mtcars[, 1:4]))
})

wts <- (tmp <- runif(nrow(mtcars))) / sum(tmp)

test_that("var.wt matches the weighted covariance matrix from stats::cov.wt", {
  expect_equal(stats::cov.wt(x = mtcars[, 1:6], wt = wts)$cov, var.wt(x =
    mtcars[, 1:6], w = wts))
})

test_that("var.wt returns the covariance matrix of 'x' if 'y' is not supplied",{
  expect_equal(var.wt(x = mtcars[, 1:4]), var.wt(x = mtcars[, 1:4], y = mtcars[,
    1:4]))
})
