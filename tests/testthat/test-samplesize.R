context("sample size")
library(breakaway)
data("apples")

est <- breakaway(apples)
test_that("sample size estimate works", {
  skip_on_cran()
  expect_null(sample_size_estimate(est$estimate, est$error, prop = 0.1))
})

test_that("sample size plot works", {
  skip_on_cran()
  expect_null(sample_size_figure(est$estimate, est$error))
})
