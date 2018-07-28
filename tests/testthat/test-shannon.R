library(breakaway)

test_that("shannon isn't crazy", {
  
  expect_warning(a <- shannon(c(0.02, 0.08)))
  
  expect_equal(shannon(c(0.2, 0.8)), a)
  
})

test_that("good_turing works", {
  
  expect_error(good_turing(c(0.2, 0.8)))
  expect_is(good_turing(apples), "alpha_estimate")
  
})

test_that("chao_shen works", {
  
  expect_error(chao_shen(c(0.2, 0.8)))
  expect_is(chao_shen(apples), "alpha_estimate")
  
})
