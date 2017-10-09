################################################################################
# Use testthat to test that nothing important breaks
################################################################################
library("breakaway")
library("testthat")

################################################################################
test_that("alpha diversity for the inbuilt datasets", {
  data("apples")
  data("hawaii")
  data("toy_otu_table")
  tables <- apply(toy_otu_table[,1:5], 2, make_frequency_count_table)
  
  datasets <- list(apples, hawaii)
  datasets <- append(datasets, tables)
  for (i in 1:length(datasets)) {
    dataset <- datasets[[i]]
    breakaway_dataset <- breakaway(dataset, output = F, answers = T, plot = F)
    lower_bound <- sum(dataset[,2])
    shannon_dataset <- shannon_better(dataset)
    
    expect_that(breakaway_dataset$est > lower_bound, is_true(),
                "a species richness estimate is negative!")
    expect_that(breakaway_dataset$seest > 0, is_true(),
                "a species richness standard error is negative!")
    expect_that(breakaway_dataset$seest > breakaway_dataset$est*0.05, is_true(),
                "a species richness standard error is too low!")
    expect_that(shannon_dataset$estimate > 0, is_true(),
                "a Shannon diversity estimate is negative!")
    expect_that(shannon_dataset$standard_error > 0, is_true(),
                "a Shannon diversity standard error is negative!")
  }
})
