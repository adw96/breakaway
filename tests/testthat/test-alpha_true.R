context("alpha_true")

library(breakaway)

test_that("alpha true functions work", {
  skip_on_cran()
  expect_is(true_shannon(c(.2,.8)), "numeric")
  expect_error(true_shannon(.2))
  expect_error(true_hill(.2), q = 2)
  expect_equal(true_hill(c(.2,.8), q = 1), exp(true_shannon(c(.2,.8))))
  expect_is(true_hill(c(.2,.8), q = 2), "numeric")
  expect_is(true_inverse_simpson(c(.2, .8)), "numeric")
  expect_is(true_gini(c(.2,.8)), "numeric")
  expect_is(true_shannon_e(c(.2,.8)), "numeric")
  expect_error(true_shannon_e(c(.2)))
})
