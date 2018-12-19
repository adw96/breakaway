#' The Good-Turing estimate of species richness
#' 
#' 
#' @param input_data An input type that can be processed by \code{convert()} or a \code{phyloseq} object
#'
#' @return An object of class \code{alpha_estimate}, or \code{alpha_estimates} for \code{phyloseq} objects
#' @export 
good_turing  <- function(input_data) {
  
  if (class(input_data) == "phyloseq") {
    if (input_data %>% otu_table %>% taxa_are_rows) {
      return(input_data %>% 
               get_taxa %>%
               apply(2, function(x) good_turing(make_frequency_count_table(x))) %>%
               alpha_estimates)
    } else {
      return(input_data %>% 
               otu_table %>%
               apply(1, function(x) good_turing(make_frequency_count_table(x))) %>%
               alpha_estimates)
    }
  }
  
  cleaned_data <- convert(input_data)
  
  the_warning <- NULL
  if (cleaned_data[1,1]!=1 || cleaned_data[1,2]==0) {
    the_warning <- "You don't have an observed singleton count.\n Chao-Shen isn't built for that data structure.\n"
  } 
  
  cc <- sum(cleaned_data[, 2])
  n <- sum(cleaned_data[, 2] * cleaned_data[, 1])
  f1 <- ifelse(cleaned_data[1,1] == 1, cleaned_data[1,2], 0)
  
  chat <- cc / (1 - f1/n)
  alpha_estimate(estimate = chat, 
                 error = NA,
                 estimand = "Shannon",
                 name = "good_turing",
                 parametric = FALSE,
                 frequentist = TRUE,
                 warnings = the_warning)
}
