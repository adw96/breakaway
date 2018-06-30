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
  
})

# datasets_to_test <- list(apples, 
#                          hawaii)
# 
# richness_estimates <- list(breakaway, 
#                            chao_bunge,
#                            chao1, 
#                            chao1_bc,
#                            wlrm_transformed,
#                            wlrm_untransformed,
#                            objective_bayes_negbin(apples, iterations=100, burn.in=50),
#                            objective_bayes_poisson(apples, iterations=100, burn.in=50))
# 
# 
# test_that("Richness estimates accept frequency count tables", {
#   
#   ### apply all functions in richness_estimates to all datasets in datasets_to_test
#   
#   for (dataset in datasets_to_test) {
#     for (richness_estimate in richness_estimates) {
#       expect_is(richness_estimate(dataset), "diversity_estimate")
#     }
#   }
#   
# })
# 
# test_that("Richness estimates accept count vectors", {
#   
#   # convert all datasets to count vectors and test
#   
#   for (dataset in datasets_to_test) {
#     count_vector_to_test <- to_count_vector(dataset)
#     for (richness_estimate in richness_estimates) {
#       expect_is(richness_estimate(count_vector_to_test), "diversity_estimate")
#     }
#   }
# })
