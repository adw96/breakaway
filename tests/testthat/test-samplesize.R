context("sample size")
library(breakaway)
data("apples")

est <- breakaway(apples)
test_that("sample size estimate works", {
  # adjusted optional arguments alpha and precision to make test run faster
  expect_null(sample_size_estimate(est$estimate, est$error, prop = 0.1, alpha = 0.15, precision = 20))
})

test_that("sample size plot works", {
  expect_null(sample_size_figure(est$estimate, est$error))
})
