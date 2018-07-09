#' The Chao-Shen estimate of Shannon diversity
#' 
#' 
#' @param input_data An input type that can be processed by \code{convert()}
#' @return An object of class \code{alpha_estimate} 
#' @export chao_shen
chao_shen  <- function(input_data) {
  
  cleaned_data <- convert(input_data)
  
  the_warning <- NULL
  if (cleaned_data[1,1]!=1 || cleaned_data[1,2]==0) {
    the_warning <- "You don't have an observed singleton count.\n Chao-Shen isn't built for that data structure.\n"
  } 
  
  estimate <- chao_shen_estimate(cleaned_data)  
  cc <- sum(cleaned_data[, 2])
  
  j_max <- length(cleaned_data[, 2])  
  
  derivative <- vector("numeric", j_max)
  for (i in 1:j_max) {
    perturbed_table <- cleaned_data
    perturbed_table[i, 2] <- perturbed_table[i, 2] + 1
    upper <- chao_shen_estimate(perturbed_table)
    perturbed_table[i, 2] <- perturbed_table[i, 2] - 2
    lower <- chao_shen_estimate(perturbed_table)
    derivative[i] <- (upper - lower)/2
  }
  
  variance_estimate <- t(derivative) %*% multinomial_covariance(cleaned_data, cc/estimate) %*% derivative
  
  alpha_estimate(estimate = estimate, 
                 error = c(ifelse(variance_estimate < 0, 0, sqrt(variance_estimate))),
                 estimand = "Shannon",
                 name = "chao_shen",
                 parametric = FALSE,
                 frequentist = TRUE,
                 warnings = the_warning)
  
}

chao_shen_estimate <- function(cleaned_data) {
  n <- sum(cleaned_data[, 2] * cleaned_data[, 1])
  f1 <- ifelse(cleaned_data[1,1] == 1, cleaned_data[1,2], 0)
  
  p_hat <- to_proportions(cleaned_data, type="frequency count")
  chat <- 1 - f1/n
  p_tilde <- chat * p_hat
  
  -sum(p_tilde * log(p_tilde) / (1 - (1 - p_tilde)^n))
  
}

multinomial_covariance <- function(my_data, chat) {
  frequencies <- my_data[, 2]
  j_max <- length(frequencies)  
  
  # divide diagonal by 2 so that when we make it symmetric we don't double count
  estimated_covariance <- diag(frequencies*(1 - frequencies / chat) / 2) 
  
  for (i in 1:(j_max-1)) {
    estimated_covariance[i, (i + 1):j_max] <- -frequencies[i]*frequencies[(i + 1):j_max]/chat
  }
  
  estimated_covariance + t(estimated_covariance)
  
}



