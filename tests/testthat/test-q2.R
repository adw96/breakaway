context("qiime2 examples")
library(breakaway)

test_that("Canonical QIIME2 Example Datasets Work", {

  data("atacama")
  x <- atacama
  a_table <- read.table(text =x,
                        skip = 0,
                        header = F,
                        row.names = NULL)[,-1]
  colnames(a_table) <- colnames(read.csv(text = x, nrows=1, skip=1, sep = "\t"))[-1]
  a_ps <- phyloseq(otu_table(a_table, taxa_are_rows = TRUE))

  expect_warning({a_breakaway <- breakaway(a_ps)})
  expect_is(a_breakaway, "alpha_estimates")
  expect_is(a_breakaway %>% summary, "data.frame")

})
