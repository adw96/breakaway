################################################################################
# Use testthat to test that nothing important breaks
################################################################################
library("breakaway")
library("testthat")

################################################################################
test_that("Nothing stupid happens with the apples dataset", {
  data("apples")
  breakaway_apples <- breakaway(apples, print = F, answers = T, plot = F)
  
  expect_that(breakaway_apples$est > 0, is_true(),
              "a species richness estimate is negative!")
  expect_that(breakaway_apples$seest > 0, is_true(),
              "a species richness standard error is negative!")
})

test_that("Nothing stupid happens with the hawaii dataset", {
  data("hawaii")
  breakaway_hawaii <- breakaway(hawaii, print = F, answers = T, plot = F)
  
  expect_that(breakaway_hawaii$est > 0, is_true(),
              "a species richness estimate is negative!")
  expect_that(breakaway_hawaii$seest > 0, is_true(),
              "a species richness standard error is negative!")
})