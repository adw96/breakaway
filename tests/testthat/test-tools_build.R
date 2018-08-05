context("sampling")
library(breakaway)
library(phyloseq)

values <- sample(0:5,250,TRUE)

OTU1 <- otu_table(values %>% matrix(nrow = 10, ncol = 25), 
                  taxa_are_rows = FALSE)
OTU2 <- otu_table(values %>% matrix(nrow = 10, ncol = 25), 
                  taxa_are_rows = TRUE)

test_that("tools_build works as intended", {
  
  expect_is(build_frequency_count_tables(OTU1), "list")
  expect_is(build_frequency_count_tables(OTU2), "list")
  
  expect_true(all(build_frequency_count_tables(OTU1) %>% names ==
                    paste("sa", 1:10, sep = "")))
  expect_true(all(build_frequency_count_tables(OTU2) %>% names ==
                    paste("sa", 1:10, sep = "")))

  expect_is(proportions_instead(values), "numeric")
})

