#' Wrapper for phyloseq objects
#' 
#' @param fn alpha diversity estimator function with \code{breakaway} to be applied to \code{physeq}
#' @param physeq \code{phyloseq} object
#' @param ... Additional arguments for fn
#' 
#' @return Object of class \code{alpha_estimates}
physeq_wrap <- function(fn, physeq, ...) {
  
  if (physeq %>% otu_table %>% taxa_are_rows) {
    otus <- physeq %>% 
      get_taxa %>% 
      t
  } else {
    otus <- physeq %>% 
      otu_table
  }
  
  ests <- otus %>%
    apply(1, function(x) fn(make_frequency_count_table(x), ...)) %>%
    alpha_estimates
  
  names(ests) <- physeq %>% sample_names
  
  ests
}


