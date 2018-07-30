library(breakaway)
library(phyloseq)

data("GlobalPatterns")

test_that("alpha_estimates works", {
  
  expect_silent({y <- alpha_estimates(breakaway(apples), 
                                      chao_bunge(apples))})
  expect_is(y, "alpha_estimates")
  expect_is(summary(y), "tbl")
  expect_is(plot(y), "ggplot")
  expect_is(plot(y, symmetric = FALSE), "ggplot")
  expect_is(plot(y, xaxis = 3:4, symmetric = FALSE), "ggplot")
  expect_equal(summary(y)[1,1] %>% unlist %>% unname, 
               breakaway(apples)$estimate)
  
})

test_that("breakaway runs on phyloseq objects", {
  
  ## this one
  expect_is(breakaway(GlobalPatterns %>% 
                        subset_samples(SampleType == "Mock")), 
            "alpha_estimates")
})
  
test_that("breakaway runs on matrices objects", {
  # because it is a data frame
  expect_warning({y <- breakaway(toy_otu_table[, 1:2])})
  expect_is(y, "alpha_estimates")
  # bc it's a matrix
  expect_warning({y <- breakaway(as.matrix(toy_otu_table[, 1:2]))})
  expect_is(y, "alpha_estimates")
  
})


test_that("alpha_estimates is robust across taxonomy", {
  
  ## this one
  expect_silent({y <- GlobalPatterns %>% 
    subset_samples(SampleType %in% c("Mock")) %>%
    tax_glom("Phylum") %>%
    sample_richness})
  
  expect_is(y, "alpha_estimates")
  expect_is(summary(y), "tbl")
  expect_is(plot(y), "ggplot")
  
  ps <- GlobalPatterns %>% 
    subset_samples(SampleType %in% c("Mock")) %>%
    tax_glom("Order")
  expect_silent({z <- ps %>%
    breakaway})
  expect_is(z, "alpha_estimates")
  expect_is(summary(z), "tbl")
  expect_is(plot(z, data = ps, color = "SampleType"), "ggplot")
  
  
  expect_is(GlobalPatterns %>% 
              subset_samples(X.SampleID %in% c("Even2")) %>%
              tax_glom("Order") %>% breakaway, 
            "alpha_estimates")
  
})

