################################################################################
# Use testthat to test that nothing important breaks
################################################################################

library(testthat)
library(breakaway)
library(phyloseq)
library(magrittr)

test_check("breakaway")
