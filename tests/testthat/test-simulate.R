context("simulate")
library(breakaway)
data("apples")


test_that("simulation code works", {
  expect_is(rnbinomtable(C = 100, size = 20, probability = .1), "data.frame")
  expect_is(rztnbinomtable(C = 100, size = 20, probability = .1), "data.frame")
})