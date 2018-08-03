context("breakaway")
library(breakaway)
library(phyloseq)
library(magrittr)
# library(testthat)
data("GlobalPatterns")

#### #### #### #### #### #### #### #### 
#### All the datasets
#### #### #### #### #### #### #### #### 
# see helper-data.R

gp_sub <- GlobalPatterns %>% 
  subset_samples(SampleType %in% c("Mock", "Even2"))

gp_order <- gp_sub %>%
  tax_glom("Order")

gp_phylum <- gp_order %>%
  tax_glom("Phylum")

datasets_ps <- list(GlobalPatterns,
                    gp_order %>% 
                      subset_samples(X.SampleID %in% c("Even2")),
                    gp_order %>% 
                      subset_samples(SampleType %in% c("Mock")),
                    gp_phylum %>% 
                      subset_samples(X.SampleID %in% c("Even2")),
                    gp_phylum %>% 
                      subset_samples(SampleType %in% c("Mock")),
                    gp_order)


#### #### #### #### #### #### #### #### 
#### All the checks
#### #### #### #### #### #### #### #### 
warnings_thrown <- function(est) {
  working <- TRUE
  if (is.nan(est$error) & is.null(est$warnings)) {
    working <- FALSE
  } 
  working
}

finite_ci <- function(est) {
  working <- TRUE
  if (is.nan(est$ci[2])) {
    working <- FALSE
  } 
  working
}

satisfies_bound <- function(est) {
  
  working <- FALSE
  if (est$estimate >= 0 | is.na(est$estimate)) {
    working <- TRUE
  } 
  working
  
}


test_that("All estimates", {
  
  mm <- lapply(X = datasets, breakaway) 
  
  mm_ps <- lapply(X = datasets_ps, breakaway)
  
  # correct_class
  lapply(X = mm, 
         FUN = expect_is, class = "alpha_estimate")
  
  lapply(X = mm_ps, 
         FUN = expect_is, class = "alpha_estimates")
  
  # valid estimates
  summaries <- lapply(X = c(mm_ps, mm), 
                      FUN = summary)

  # shouldn't be true for breakaway_nof1, the 9th estimate
  expect_true(lapply(X = lapply(mm, 
                                function(x) x$estimate), FUN = is.na) %>%
                lapply(function(x) !x) %>%
                unlist %>% all)
  
  expect_true(lapply(FUN = warnings_thrown, mm) %>%
                unlist %>% all)
  
  expect_true(lapply(FUN = satisfies_bound, mm) %>%
                unlist %>% all)
  
  ## dataset 4
  expect_true(lapply(FUN = finite_ci, mm) %>%
                unlist %>% all)
  
})
