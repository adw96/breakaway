#' @export good_turing
good_turing  <- function(my_data) {
  
  cleaned_data <- check_format(my_data)
  
  if (cleaned_data[1,1]!=1 || cleaned_data[1,2]==0) {
    warning("You don't have an observed singleton count.\n Chao-Shen isn't built for that data structure.\n")
  } 
  
  cc <- sum(cleaned_data[, 2])
  n <- sum(cleaned_data[, 2] * cleaned_data[, 1])
  f1 <- ifelse(cleaned_data[1,1] == 1, cleaned_data[1,2], 0)
  
  chat <- cc / (1 - f1/n)
  list("est" = chat, 
       "se" = NA)
}
