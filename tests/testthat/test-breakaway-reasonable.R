context("breakaway_reasonable")
library(breakaway)

data("toy_otu_table")

set.seed(1)
test_that("breakaway gives reasonable std errors", {
  ### Sample 111 was a particular problem after change
  selected_samples <- c(111, sample(1:143, size = 20, replace = F))
  
  # this failed before  any(is.infinite(ratiovars)) was added
  for (i in selected_samples) {
    y <- breakaway(toy_otu_table[, i])
    expect_true(y$error >= 2 | is.na(y$estimate),
                info = paste("std error is too small for Sample", i))
  }
})

test_that("previously crazy estimates are fine", {
  
  # this is fine: plot is just a disaster
  library(phyloseq)
  data(GlobalPatterns)
  
  bad <- (GlobalPatterns %>% otu_table)[,1] %>% c
  ba <- bad %>% breakaway
  
  expect_true(ba$error >= 5 | is.na(ba$estimate),
              info = paste("std error is too small for first GlobalPatterns"))
  
})
