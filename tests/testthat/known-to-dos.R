# A list of things that Amy wants to fix

test_that("toy_otu_table inherits from phyloseq", {
  
  data("toy_otu_table")
  expect_true(toy_otu_table %>% class %in% c("phyloseq", "otu_table"))
  
})

test_that("toy_otu_table inherits from phyloseq", {
  
  data("toy_otu_table")
  expect_true(toy_otu_table %>% class %in% c("phyloseq", "otu_table"))
  
})

test_that("wlrm_untransformed is occasionally NULL", {
  
  dataset <- make_frequency_count_table(toy_otu_table[, 126])
  wlrmut_dataset <- wlrm_untransformed(fc, answers = T, print=F)
  expect_false(wlrmut_dataset, is.null)
  
  dataset <- make_frequency_count_table(toy_otu_table[, 31])
  wlrmut_dataset <- wlrm_untransformed(fc, answers = T, print=F)
  expect_false(wlrmut_dataset, is.null)
  
  
})