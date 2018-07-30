context("phyloseq")
library(breakaway)
library(phyloseq)
data("GlobalPatterns")

test_that("estimates accept phyloseq objects", {
  
  
  richness_estimates <- list(sample_richness,
                             chao1,
                             chao1_bc,
                             breakaway,
                             chao_bunge,
                             wlrm_transformed,
                             wlrm_untransformed,
                             chao1_bc,
                             good_turing)
  # breakaway_nof1
  
  for (richness_estimate in richness_estimates) {
    
    y <- richness_estimate(GlobalPatterns) 
    
    expect_is(y, "alpha_estimates")
    expect_is(plot(y), "ggplot")
    
  }
  
})
