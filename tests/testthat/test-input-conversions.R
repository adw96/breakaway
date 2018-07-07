#### Check we can move between different types of input data and return the same values

data("apples")
data("hawaii")

test_that("Initial conversions work", {
  
  expect_equal(convert(apples), apples)
  expect_equal(convert(hawaii), hawaii)
  # convert("data/butterfly.csv")
  expect_error(convert(toy_otu_table))
  expect_error(convert(toy_metadata))
  expect_error(convert(toy_taxonomy))
  
  
  expect_error(suppressWarnings({convert(cbind(c(1:5,7,6), 20:14))}))
  
})

datasets_to_test <- list(apples,
                         hawaii)

richness_estimates <- list(chao1,
                           breakaway,
                           chao_bunge,
                           wlrm_transformed,
                           wlrm_untransformed,
                           # objective_bayes_negbin(apples, iterations=100, burn.in=50),
                           # objective_bayes_poisson(apples, iterations=100, burn.in=50),
                           chao1_bc)

test_that("Richness estimates output the correct type", {
  
  ### apply all functions in richness_estimates to all datasets in datasets_to_test
  
  for (dataset in datasets_to_test) {
    for (richness_estimate in richness_estimates) {
      expect_is(richness_estimate(dataset), "alpha_estimate")
    }
  }
  
})
