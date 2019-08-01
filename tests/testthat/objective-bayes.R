test_that("objective Bayes functions work", {
  
  lower_bound <- sum(apples[11:20,2])
  subset_apples <- apples[c(1, 11:20), ]
  expect_is({y <- objective_bayes_geometric(apples[11:20, ], iterations=5, plot=F, answers=T, output=F); y}, 
            "list")
  expect_gte(y$results["mean.C"], 
            "list")
  
  
})




