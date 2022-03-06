context("shannon")
library(breakaway)

test_that("shannon isn't crazy", {
  skip_on_cran()

  expect_equal(true_shannon(c(0.5, 0.5)),
               log(2))

})

test_that("good_turing works", {
  skip_on_cran()

  expect_error(good_turing(c(0.2, 0.8)))
  expect_is(good_turing(apples), "alpha_estimate")

})

test_that("chao_shen works", {
  skip_on_cran()

  expect_error(chao_shen(c(0.2, 0.8)))
  expect_is(chao_shen(apples), "alpha_estimate")

})

