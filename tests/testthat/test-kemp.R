context("robustify breakaway")
library(breakaway)
library(phyloseq)
library(magrittr)
data("GlobalPatterns")

test_that("alpha_estimates is robust across taxonomy", {


  ps <- y <- GlobalPatterns %>%
    subset_samples(SampleType %in% c("Mock")) %>%
    tax_glom("Phylum")

  sum <- ps %>% breakaway %>% summary
  
  expect_true(all(!is.na(sum$estimate)))

  plot(GlobalPatterns %>% breakaway)
  
})