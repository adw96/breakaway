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

library(phyloseq)
data(GlobalPatterns)
alphas <- breakaway(GlobalPatterns)

alphas_noest <- alphas
fn1 <- function(x) {
  x$estimate <- NA
  return(x)
}
alphas_noest <- lapply(alphas_noest, fn1)
alphas_noest <- alpha_estimates(alphas_noest)

test_that("plot alpha estimates", {
  expect_is(plot(alphas, trim_plot = FALSE), "ggplot")
  expect_is(plot(alphas, physeq = GlobalPatterns, color = "SampleType", shape = "SampleType"), "ggplot")
  expect_error(plot(alphas_noest))
  expect_error(plot(alphas, physeq = GlobalPatterns, color = "SampleType", shape = 1))
  expect_error(plot(alphas, physeq = GlobalPatterns, color = 1))
})
