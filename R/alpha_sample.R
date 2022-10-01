fc_to_proportions <- function(fc) {
  
  counts_long <- rep(x = fc$index, times = fc$frequency)
  nn <- sum(counts_long)
  
  stopifnot(sum(fc$frequency * fc$index) == nn)
  
  counts_long/nn
  
}

#' Plug-in Shannon diversity
#' 
#' This function implements the plug-in Shannon diversity
#' 
#' @param input_data An input type that can be processed by \code{convert()} or a \code{phyloseq} object
#'
#' @return An object of class \code{alpha_estimate}, or \code{alpha_estimates} for \code{phyloseq} objects
#' @examples
#' sample_shannon(apples)
#' 
#' @export 
sample_shannon <- function(input_data) {
  
  if ("phyloseq" %in% class(input_data)) {
    if (input_data %>% otu_table %>% taxa_are_rows) {
      return(input_data %>% 
               get_taxa %>%
               apply(2, function(x) sample_shannon(make_frequency_count_table(x))) %>%
               alpha_estimates)
    } else {
      return(input_data %>% 
               otu_table %>%
               apply(1, function(x) sample_shannon(make_frequency_count_table(x))) %>%
               alpha_estimates)
    }
  }
  
  obsd_props <- input_data %>% convert %>% fc_to_proportions
  
  my_est <- true_shannon(obsd_props)
  
  ## construct diversity_estimate
  alpha_estimate(estimate = my_est,
                 error = 0,
                 estimand = "shannon",
                 name = "Plug-in",
                 interval = NA,
                 type = NA,
                 model = "none",
                 frequentist = TRUE,
                 parametric = FALSE,
                 reasonable = FALSE,
                 interval_type = NA)
}

#' Plug-in Shannon's E ("Equitability")
#' 
#' This function implements the plug-in Shannon's E
#' 
#' @param input_data An input type that can be processed by \code{convert()} or a \code{phyloseq} object
#'
#' @return An object of class \code{alpha_estimate}, or \code{alpha_estimates} for \code{phyloseq} objects
#' @examples
#' sample_shannon_e(apples)
#' 
#' @export 
sample_shannon_e <- function(input_data) {
  
  if ("phyloseq" %in% class(input_data)) {
    if (input_data %>% otu_table %>% taxa_are_rows) {
      return(input_data %>% 
               get_taxa %>%
               apply(2, function(x) sample_shannon_e(make_frequency_count_table(x))) %>%
               alpha_estimates)
    } else {
      return(input_data %>% 
               otu_table %>%
               apply(1, function(x) sample_shannon_e(make_frequency_count_table(x))) %>%
               alpha_estimates)
    }
  }
  
  obsd_props <- input_data %>% convert %>% fc_to_proportions
  
  my_est <- true_shannon_e(obsd_props)
  
  ## construct diversity_estimate
  alpha_estimate(estimate = my_est,
                 error = 0,
                 estimand = "equitability",
                 name = "Plug-in",
                 interval = NA,
                 type = NA,
                 model = "none",
                 frequentist = TRUE,
                 parametric = FALSE,
                 reasonable = FALSE,
                 interval_type = NA)
}

#' Plug-in Simpson diversity
#' 
#' This function implements the plug-in Simpson diversity
#' 
#' @param input_data An input type that can be processed by \code{convert()} or a \code{phyloseq} object
#'
#' @return An object of class \code{alpha_estimate}, or \code{alpha_estimates} for \code{phyloseq} objects
#' @examples
#' sample_simpson(apples)
#' 
#' @export 
sample_simpson <- function(input_data) {
  
  if ("phyloseq" %in% class(input_data)) {
    if (input_data %>% otu_table %>% taxa_are_rows) {
      return(input_data %>% 
               get_taxa %>%
               apply(2, function(x) sample_simpson(make_frequency_count_table(x))) %>%
               alpha_estimates)
    } else {
      return(input_data %>% 
               otu_table %>%
               apply(1, function(x) sample_simpson(make_frequency_count_table(x))) %>%
               alpha_estimates)
    }
  }
  
  obsd_props <- input_data %>% convert %>% fc_to_proportions
  
  my_est <- true_simpson(obsd_props)
  
  ## construct diversity_estimate
  alpha_estimate(estimate = my_est,
                 error = 0,
                 estimand = "simpson",
                 name = "Plug-in",
                 interval = NA,
                 type = NA,
                 model = "none",
                 frequentist = TRUE,
                 parametric = FALSE,
                 reasonable = FALSE,
                 interval_type = NA)
}

#' Plug-in Inverse Simpson diversity
#' 
#' This function implements the plug-in Inverse Simpson diversity
#' 
#' @param input_data An input type that can be processed by \code{convert()} or a \code{phyloseq} object
#'
#' @return An object of class \code{alpha_estimate}, or \code{alpha_estimates} for \code{phyloseq} objects
#' @examples
#' sample_inverse_simpson(apples)
#' 
#' @export 
sample_inverse_simpson <- function(input_data) {
  
  if ("phyloseq" %in% class(input_data)) {
    if (input_data %>% otu_table %>% taxa_are_rows) {
      return(input_data %>% 
               get_taxa %>%
               apply(2, function(x) sample_inverse_simpson(make_frequency_count_table(x))) %>%
               alpha_estimates)
    } else {
      return(input_data %>% 
               otu_table %>%
               apply(1, function(x) sample_inverse_simpson(make_frequency_count_table(x))) %>%
               alpha_estimates)
    }
  }
  
  obsd_props <- input_data %>% convert %>% fc_to_proportions
  
  my_est <- true_inverse_simpson(obsd_props)
  
  ## construct diversity_estimate
  alpha_estimate(estimate = my_est,
                 error = 0,
                 estimand = "inverse simpson",
                 name = "Plug-in",
                 interval = NA,
                 type = NA,
                 model = "none",
                 frequentist = TRUE,
                 parametric = FALSE,
                 reasonable = FALSE,
                 interval_type = NA)
}