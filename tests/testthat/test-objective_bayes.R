test_that("objective bayes", {
  z <- rnbinomtable(20, 5, 0.5)
  expect_is(z, "data.frame")
  expect_is(objective_bayes_geometric(z, 
                                      iterations = 10, 
                                      burn.in = 2,
                                      plot = F,
                                      output = F,
                                      answers = T),
            "list")
  expect_is(objective_bayes_negbin(z, 
                                   iterations = 10, 
                                   burn.in = 2,
                                   plot = F,
                                   output = F,
                                   answers = T),
            "list")
  expect_is(objective_bayes_poisson(z, 
                                    iterations = 10, 
                                    burn.in = 2,
                                    plot = F,
                                    output = F,
                                    answers = T),
            "list")
  # expect_is(objective_bayes_geometric(z, 
  #                                     iterations = 10, 
  #                                     burn.in = 2,
  #                                     plot = F,
  #                                     output = F,
  #                                     answers = T),
  #           "list")
})
