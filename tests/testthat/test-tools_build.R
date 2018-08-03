context("sampling")
library(breakaway)
library(phyloseq)
data("GlobalPatterns")

OTU1 <- otu_table(matrix(sample(0:5,250,TRUE),25,10), taxa_are_rows = FALSE)
OTU2 <- otu_table(matrix(sample(0:5,250,TRUE),10,25), taxa_are_rows = TRUE)

my_freq <- build_frequency_count_tables(OTU1)

test_that("build frequency count table", {
  expect_is(my_freq, "list")
  expect_is(build_frequency_count_tables(OTU2), "list")
})

tmp <- my_freq$sp1
tmp2 <- tmp
tmp2[,1] <- as.character(tmp[,1])
test_that("count_freq_indices works", {
  expect_is(breakaway:::convert_freq_indices(tmp), "matrix")
  expect_is(breakaway:::convert_freq_indices(tmp2), "matrix")
  expect_equal(breakaway:::convert_freq_indices(tmp), breakaway:::convert_freq_indices(tmp2))
})

test_that("proportions_instead works", {
  expect_is(proportions_instead(tmp), "numeric")
  expect_is(proportions_instead(matrix(sample(0:5,250,TRUE),25,10)), "numeric")
})



tax1 <- tax_table(matrix("abc", 30, 8))
map1 <- data.frame( matrix(sample(0:3,250,TRUE),25,10),
                    matrix(sample(c("a","b","c"),150,TRUE), 25, 6) )
map1 <- sample_data(map1)
exam1 <- phyloseq(OTU1, map1, tax1)