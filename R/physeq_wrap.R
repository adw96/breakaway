#' Wrapper for phyloseq objects
#' 
#' @param fn alpha diversity estimator function with \code{breakaway} to be applied to \code{physeq}
#' @param physeq A \code{phyloseq} object, or an object of class \code{otu_table} 
#' @param ... Additional arguments for fn
#' 
#' @return Object of class \code{alpha_estimates}
physeq_wrap <- function(fn, physeq, ...) {
  
  if (class(physeq) == "otu_table") {
    ot <- physeq
  } else if (class(physeq) == "phyloseq") {
    ot <- physeq %>% otu_table
  } else {
    stop(paste("Unknown type passed to physeq_wrap; object is of class",
               class(physeq)))
  }
  
  if (ot %>% taxa_are_rows) {
    otus <- physeq %>% get_taxa %>% t
  } else {
    otus <- ot
  }
  
  ests <- otus %>%
    apply(1, function(x) fn(make_frequency_count_table(x), ...)) %>%
    alpha_estimates
  
  names(ests) <- physeq %>% sample_names
  
  ests
}


