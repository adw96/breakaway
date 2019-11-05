context("test-pathological")
library(breakaway)

test_that("issues that are fixed on github stay fixed", {
  
  issue70 <- cbind(c(1, 2, 3, 4, 5, 6, 9, 10, 14, 16, 18, 21, 22, 24, 25, 36, 51, 52, 80, 91, 774, 1138, 1205, 24789), 
                   c(2, 1, 2, 1, 1, 2, rep(1, 6), 2, rep(1, 11))) 
  expect_is(breakaway(issue70), "alpha_estimate")
  expect_warning(chao_bunge(issue70, cutoff=10))

  issue68 <- structure(list(Var1 = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 
                                16, 19, 20, 22, 23, 24, 25, 28), Freq = c(69, 15, 8, 23, 4, 1, 
                                                                          5, 3, 3, 2, 3, 1, 2, 1, 1, 1, 1, 2, 1, 2)), row.names = c(NA, 
                                                                                                                                    -20L), class = "data.frame")
  expect_true(breakaway(issue68)$estimate > 200)

})
