context("test all estimates")
library(breakaway)
library(magrittr)
library(phyloseq)
data("GlobalPatterns")
set.seed(1)
# library(testthat)

#### #### #### #### #### #### #### #### 
#### All the estimates
#### #### #### #### #### #### #### #### 
richness_estimates <- list(sample_richness,
                           chao1,
                           chao1_bc,
                           chao_bunge,
                           wlrm_transformed,
                           wlrm_untransformed,
                           chao1_bc,
                           breakaway_nof1)
# Missing: breakaway_nof1, good_turing

#### #### #### #### #### #### #### #### 
#### All the datasets
#### #### #### #### #### #### #### #### 
# 
# random_samples <- sample(x = 1:143, size = 5, replace = F)
# tables <- apply(toy_otu_table[,random_samples], 2, make_frequency_count_table)
# 
# datasets <- list(hawaii,
#                  apples)
# 
# datasets %<>% c(tables)
# 
datasets_ps <- list(GlobalPatterns,
                    GlobalPatterns %>%
                      subset_samples(SampleType %in% c("Mock")) %>%
                      tax_glom("Phylum"))

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

execute <- function(est, df) {
  est(df)
}

test_that("All estimates", {
  
  mm <- lapply(richness_estimates, 
               function(a){ lapply(datasets, function(b){execute(a,b)}) }) %>%
    unlist(recursive = F)
  
  mm_ps <- lapply(richness_estimates, 
                  function(a){ lapply(datasets_ps, function(b){
                    execute(a,b)}) }) %>%
    unlist(recursive = F)
  
  # correct_class
  lapply(X = mm, 
         FUN = expect_is, class = "alpha_estimate")
  
  lapply(X = mm_ps, 
         FUN = expect_is, class = "alpha_estimates")
  
  # valid estimates
  summaries <- lapply(X = c(mm_ps, mm), 
                      FUN = summary)
  
  expect_true(lapply(FUN = satisfies_bound, mm) %>%
                unlist %>% all)
  
})
