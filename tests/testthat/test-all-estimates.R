context("Test everything")
library(breakaway)
library(phyloseq)
library(magrittr)
# library(testthat)
data("GlobalPatterns")

#### #### #### #### #### #### #### #### 
#### All the estimates
#### #### #### #### #### #### #### #### 
richness_estimates <- list(sample_richness,
                           chao1,
                           chao1_bc,
                           breakaway,
                           chao_bunge,
                           wlrm_transformed,
                           wlrm_untransformed,
                           chao1_bc,
                           good_turing)
# Missing: breakaway_nof1

#### #### #### #### #### #### #### #### 
#### All the datasets
#### #### #### #### #### #### #### #### 

df1 <- cbind(c(1,2,3,6), c(8,4,2,1))

random_samples <- sample(x = 1:143, size = 5, replace = F)
tables <- apply(toy_otu_table[,random_samples], 2, make_frequency_count_table)

datasets <- list(hawaii,
                 apples,
                 df1)
datasets %<>% c(tables)

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
  
  expect_true(lapply(X = lapply(mm[-43], function(x) x$estimate), FUN = is.na) %>%
                lapply(function(x) !x) %>%
                unlist %>% all)
  
  expect_true(lapply(FUN = warnings_thrown, mm) %>%
                unlist %>% all)
  
  expect_true(lapply(FUN = satisfies_bound, mm) %>%
                unlist %>% all)
  
  expect_true(lapply(FUN = finite_ci, mm) %>%
                unlist %>% all)
  
})


#### #### #### #### #### #### #### #### 
#### Singleton-missing checks
#### #### #### #### #### #### #### #### 

df2 <- cbind(c(3,6,9,12,15,16), c(3,2,1,15,1,5))
datasets_nof1 <- list(df2)
breakaway(df2)

test_that("Harder estimates", {
  
  execute <- function(est, df) {
    est(df)
  }
  
  mm <- lapply(richness_estimates, 
               function(a){ lapply(datasets_nof1, function(b){execute(a,b)}) }) %>%
    unlist(recursive = F)
  
  # correct_class
  lapply(X = mm, 
         FUN = expect_is, class = "alpha_estimate")
  
  # valid estimates
  expect_true(lapply(FUN = warnings_thrown, mm) %>%
                unlist %>% all)
  
})
