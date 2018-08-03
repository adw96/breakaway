#' Wrapper for phyloseq objects
#' 
#' @param fn Richness estimator function with \code{breakaway} to be applied to \code{physeq}
#' @param physeq \code{phyloseq} object
#' 
#' @return Object of class \code{alpha_estimates}
#' @export
physeq_wrap <- function(fn, physeq) {
  if (physeq %>% otu_table %>% taxa_are_rows) {
    return(physeq %>% 
             get_taxa %>%
             apply(2, function(x) fn(make_frequency_count_table(x))) %>%
             alpha_estimates)
  } else {
    return(physeq %>% 
             otu_table %>%
             apply(1, function(x) fn(make_frequency_count_table(x))) %>%
             alpha_estimates)
  }
}


