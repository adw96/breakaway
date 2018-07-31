context("datatsets run")
library(breakaway)

data("apples")
data("hawaii")
data("toy_otu_table")

random_samples <- sample(x = 1:143, size = 5, replace = F)
tables <- apply(toy_otu_table[,random_samples], 2, make_frequency_count_table)
datasets <- list(apples, hawaii)
datasets <- append(datasets, tables)
dataset_names <- c("apples", 
                   "hawaii", 
                   paste("toy_otu_table_sample", 
                         random_samples, sep = ""))


test_that("datasets load", {
  
  expect_equal(apples[1, 2], 277)
  expect_equal(hawaii[1, 2], 690)
  
  expect_equal(dim(toy_otu_table)[2], 143) # 143 samples
  
})

richness_estimates <- list(chao1,
                           breakaway,
                           chao_bunge,
                           wlrm_transformed,
                           wlrm_untransformed,
                           chao1_bc)

test_that("richness for the inbuilt datasets", {
  
  expect_equal(length(datasets), length(dataset_names))
  
  for (i in 1:length(datasets)) {
    
    dataset <- datasets[[i]]
    lower_bound <- sum(dataset[,2])
    
    for (richness_estimate in richness_estimates) {
      
      y <- richness_estimate(dataset) 
      
      expect_is(y, "alpha_estimate")
      
      expect_true(y$estimate >= lower_bound | is.na(y$estimate), 
                  info = paste("estimate is negative for dataset", 
                               dataset_names[i]))
      
      expect_true(y$error >= 0 | is.na(y$estimate), 
                  info = paste("std error is negative for dataset", 
                               dataset_names[i]))
      
      expect_true(sample_shannon(dataset)$estimate >= 0, 
                  info = paste("Plug-in shannon is negative for dataset", dataset_names[i]))
      
    }
  }
})




