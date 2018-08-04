context("s3 alpha estimate")

library(breakaway)
library(phyloseq)
data(apples)
tmp <- breakaway(apples)
data(GlobalPatterns)


test_that("alpha estimate s3 functions work", { 
  expect_null(print(tmp))
  expect_is(plot(tmp), "ggplot")
})

test_that("alpha estimate function breaks", {
  expect_error(alpha_estimates(c(1,2,3),c(1,2,3)))
})


fn1 <- function(x) {
  x$estimate <- NA
  return(x)
}

test_that("plot alpha estimates", {
  alphas <- breakaway(GlobalPatterns %>% 
                        subset_samples(SampleType == "Soil"))
  
  expect_is(plot(alphas, trim_plot = FALSE), "ggplot")
  
  ## help
  expect_is(plot(alphas, physeq = GlobalPatterns %>% 
                   subset_samples(SampleType == "Soil"), 
                 color = "SampleType", shape = "SampleType"), "ggplot")
  
  alphas_noest <- alphas
  alphas_noest <- lapply(alphas_noest, fn1)
  alphas_noest <- alpha_estimates(alphas_noest)
  expect_error(plot(alphas_noest))
  
  expect_error(plot(alphas, physeq = GlobalPatterns %>% 
                      subset_samples(SampleType == "Soil"), 
                    color = "SampleType", shape = 1))
  
  expect_error(plot(alphas, physeq = GlobalPatterns, color = 1))
})
