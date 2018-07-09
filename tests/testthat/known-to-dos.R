# A list of things that Amy wants to fix

test_that("toy_otu_table inherits from phyloseq", {
  
  data("toy_otu_table")
  expect_true(toy_otu_table %>% class %in% c("phyloseq", "otu_table"))
  
})


test_that("use aes_string()", {
  # check out
  # https://stackoverflow.com/questions/9439256/how-can-i-handle-r-cmd-check-no-visible-binding-for-global-variable-notes-when
  expect_true(FALSE)
})

test_that("make plot alpha_estimates have x axis", {
  expect_true(FALSE)
})

test_that("make aggregated phyloseq object for testing", {
  expect_true(FALSE)
})

test_that("breakaway on phyloseq objects runs in parallel", {
  expect_true(FALSE)
})
