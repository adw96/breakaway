################################################################################
# Use testthat to test that nothing important breaks
################################################################################
library("breakaway")
library("testthat")

################################################################################
test_that("Can transform_sample_counts of an OTU table that is either orientation", {
  data("apples")
  breakaway_apples <- breakaway(apples, print = F, answers = T, plot = F)
  
  expect_that(breakaway_apples$est > 0, is_true(),
              "a species richness estimate is negative!")
  expect_that(breakaway_apples$seest > 0, is_true(),
              "a species richness standard error is negative!")
})