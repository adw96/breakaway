context("design matrix")
library(breakaway)
library(phyloseq)
data(GlobalPatterns)
test_that("sample size estimate works", {
  skip_on_cran()
  expect_is(make_design_matrix(GlobalPatterns, "SampleType"), "matrix")
})

