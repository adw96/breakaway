# A list of things that Amy wants to fix

test_that("toy_otu_table inherits from phyloseq", {
  
  data("toy_otu_table")
  expect_true(toy_otu_table %>% class %in% c("phyloseq", "otu_table"))
  
})

test_that("make aggregated phyloseq object for testing", {
  expect_true(FALSE)
})

test_that("breakaway on phyloseq objects runs in parallel", {
  expect_true(FALSE)
})

test_that("GlobalPatterns runs without warnings", {
  x <- GlobalPatterns %>% breakaway
  expect_is(x,
              "alpha_estimates")
  
  ### All estimates are sensible
})
