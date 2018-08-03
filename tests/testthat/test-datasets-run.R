context("datatsets run")
library(breakaway)

data("apples")
data("hawaii")
data("toy_otu_table")

test_that("datasets load", {
  
  expect_equal(apples[1, 2], 277)
  expect_equal(hawaii[1, 2], 690)
  expect_equal(dim(toy_otu_table)[2], 143) # 143 samples
  
})
