test_that("shannon isn't crazy", {
  
  expect_warning(a <- shannon(c(0.02, 0.08)))
  
  expect_equal(shannon(c(0.2, 0.8)), a)
})
