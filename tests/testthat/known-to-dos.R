# A list of things that Amy wants to fix

test_that("toy_otu_table inherits from phyloseq", {
  
  data("toy_otu_table")
  expect_true(toy_otu_table %>% class %in% c("phyloseq", "otu_table"))
  
})

test_that("toy_otu_table inherits from phyloseq", {
  
  data("toy_otu_table")
  expect_true(toy_otu_table %>% class %in% c("phyloseq", "otu_table"))
  
})

test_that("breakaway runs on phyloseq objects and matrices", {
  
  expect_is(breakaway(toy_otu_table), "alpha_estimate")
  
  
})

test_that("crazy estimates give warnings", {
  
  expect_warning((GlobalPatterns %>% otu_table)[,1] %>% c %>% breakaway)

})
