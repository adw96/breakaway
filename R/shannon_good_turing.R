#' The Good-Turing estimate of Shannon diversity
#' 
#' 
#' @param input_data An input type that can be processed by \code{convert()}
#' @return An object of class \code{alpha_estimate} 
#' @export good_turing
good_turing  <- function(input_data) {
  
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
