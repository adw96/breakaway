context("s3 alpha estimate")

library(breakaway)
data(apples)
tmp <- breakaway(apples)


test_that("alpha estimate s3 functions work", { 
  expect_null(print(tmp))
  expect_is(plot(tmp), "ggplot")
})

test_that("alpha estimate function breaks", {
  expect_error(alpha_estimates(c(1,2,3),c(1,2,3)))
})