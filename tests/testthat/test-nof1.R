context("breakaway nof1")
library(breakaway)
library(phyloseq)
data("apples")
apples <- apples[-1,]
# nof1 can't currently take phyloseq!
test_that("breakaway nof1 works", {
  ### Sample 111 was a particular problem after change
  expect_is(breakaway_nof1(apples), "alpha_estimate")
})