context("sampling")
library(breakaway)
library(phyloseq)
data("GlobalPatterns")

OTU1 <- otu_table(matrix(sample(0:5,250,TRUE),25,10), taxa_are_rows=FALSE)
tax1 <- tax_table(matrix("abc", 30, 8))
map1 <- data.frame( matrix(sample(0:3,250,TRUE),25,10),
matrix(sample(c("a","b","c"),150,TRUE), 25, 6) )
map1 <- sample_data(map1)
exam1 <- phyloseq(OTU1, map1, tax1)

test_that("sampling functions work", {
  skip_on_cran()
  expect_is(sample_shannon(GlobalPatterns), "alpha_estimates")
  expect_is(sample_shannon_e(GlobalPatterns), "alpha_estimates")
  expect_is(sample_simpson(GlobalPatterns), "alpha_estimates")
  expect_is(sample_inverse_simpson(GlobalPatterns), "alpha_estimates")
  expect_is(sample_shannon(exam1), "alpha_estimates")
  expect_is(sample_shannon_e(exam1), "alpha_estimates")
  expect_is(sample_simpson(exam1), "alpha_estimates")
  expect_is(sample_inverse_simpson(exam1), "alpha_estimates")
})
