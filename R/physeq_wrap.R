#' Wrapper for phyloseq objects
#' 
#' @param input_data \code{phyloseq} object
#' 
#' @return Object of class \code{alpha_estimates}
#' @export
physeq_wrap <- function(input_data) {
  if (input_data %>% otu_table %>% taxa_are_rows) {
    return(input_data %>% 
             get_taxa %>%
             apply(2, function(x) chao1(make_frequency_count_table(x))) %>%
             alpha_estimates)
  } else {
    return(input_data %>% 
             otu_table %>%
             apply(1, function(x) chao1(make_frequency_count_table(x))) %>%
             alpha_estimates)
  }
}


