#' Sample richness estimator
#' 
#' This function implements the sample richness estimate, which is the number of non-zero taxa per sample.
#' 
#' 
#' @param input_data An input type that can be processed by \code{convert()} or a \code{phyloseq} object
#'
#' @return An object of class \code{alpha_estimate}, or \code{alpha_estimates} for \code{phyloseq} objects
#' @examples
#' 
#' 
#' sample_richness(apples)
#' 
#' 
#' @export 
sample_richness <- function(input_data) {
  if (class(input_data) == "phyloseq") {
    if (input_data %>% otu_table %>% taxa_are_rows) {
      return(input_data %>% 
               get_taxa %>%
               apply(2, function(x) sample_richness(make_frequency_count_table(x))) %>%
               alpha_estimates)
    } else {
      return(input_data %>% 
               otu_table %>%
               apply(1, function(x) sample_richness(make_frequency_count_table(x))) %>%
               alpha_estimates)
    }
  }
  
  my_est <- as.numeric(colSums(convert(input_data))["frequency"])
  
  ## construct diversity_estimate
  alpha_estimate(estimate = my_est,
                 error = 0,
                 estimand = "richness",
                 name = "Plug-in",
                 interval = c(my_est, my_est),
                 type = NA,
                 model = "none",
                 frequentist = TRUE,
                 parametric = FALSE,
                 reasonable = FALSE,
                 interval_type = NA)
}