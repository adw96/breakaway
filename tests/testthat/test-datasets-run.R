data("apples")
data("hawaii")
data("toy_otu_table")

random_samples <- sample(x=1:143, size=5, replace=F)
tables <- apply(toy_otu_table[,random_samples], 2, make_frequency_count_table)
datasets <- list(apples, hawaii)
datasets <- append(datasets, tables)
dataset_names <- c("apples", "hawaii", paste("toy_otu_table_sample", random_samples, sep=""))


test_that("datasets load", {
  
  expect_equal(apples[1, 2], 277)
  expect_equal(hawaii[1, 2], 690)
  
  expect_equal(dim(toy_otu_table)[2], 143) # 143 samples
  
})



test_that("richness for the inbuilt datasets", {
  
  expect_equal(length(datasets), length(dataset_names))
  
  for (i in 1:length(datasets)) {
    
    dataset <- datasets[[i]]
    
    breakaway_dataset <- breakaway(dataset, output = F, answers = T, plot = F)
    cb_dataset <- chao_bunge(dataset, answers = T)
    wlrmt_dataset <- wlrm_transformed(dataset)
    wlrmut_dataset <- wlrm_untransformed(dataset, answers = T)
    
    lower_bound <- sum(dataset[,2])
    
    #### Richness
    expect_true(breakaway_dataset$est >= lower_bound, 
                info = paste("breakaway estimate is negative for dataset", dataset_names[i]))
    
    expect_true(cb_dataset$est >= lower_bound, 
                info = paste("chao_bunge estimate is negative for dataset", dataset_names[i]))
    
    expect_true(wlrmt_dataset$seest >= 0, 
                info = paste("wlrm_transformed std error is negative for dataset", dataset_names[i]))
    
    expect_true(wlrmut_dataset$seest >= 0, 
                info = paste("wlrm_untransformed std error is negative for dataset", dataset_names[i]))
    
    #### Shannon
    shannon_dataset <- shannon(dataset)
    expect_true(shannon_dataset >= 0, 
                info = paste("shannon is negative for dataset", dataset_names[i]))
    
    # expect_true(breakaway_dataset$seest > breakaway_dataset$est*0.05, 
    # info = "a species richness standard error is too low!")
    
    # expect_true(shannon_dataset$standard_error > 0, 
    #             info = paste("shannon std error is negative for dataset", dataset_names[i]))
  }
})




