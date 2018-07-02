#' Bias-corrected Chao1 species richness estimator
#' 
#' This function implements the bias-corrected Chao1 richness estimate.
#' 
#' 
#' @param input_data An input type that can be processed by \code{convert()}
#' @param output Deprecated; only for backwards compatibility
#' @param answers Deprecated; only for backwards compatibility
#'
#' @return An object of class \code{alpha_estimate}
#' @note The authors of this package strongly discourage the use of this
#' estimator.  It is only valid when you wish to assume that every taxa has
#' equal probability of being observed. You don't really think that's possible,
#' do you?
#' @author Amy Willis
#' @examples
#' 
#' 
#' chao1_bc(apples)
#' 
#' 
#' @export chao1_bc
chao1_bc <- function(input_data, output=NULL, answers=NULL) {
  
  my_data <- convert(input_data)
  
  # TODO: this is a stupid way of doing it, find a better one 
  index  <- 1:max(my_data[,1])
  frequency_index <- rep(0, length(index))
  frequency_index[my_data[,1]] <- my_data[,2]
  f1  <- frequency_index[1]
  f2 <- frequency_index[2]
  n <- sum(frequency_index)
  
  f0 <- f1*(f1-1)/(2*(f2+1))
  diversity <- n + f0
  
  diversity_se <- sqrt(f1*(f1-1)/(2*(f2+1)) + f1*(2*f1-1)^2/(4*(f2+1)^2) + f1^2*f2*(f1-1)^2/(4*(f2+1)^4))
  
  # TODO: write a function to do this
  d <- exp(1.96*sqrt(log(1 + diversity_se^2 / f0)))
  
  ## construct diversity_estimate
  alpha_estimate(estimate = diversity,
                 error = diversity_se,
                 estimand = "richness",
                 name = "chao1_bc",
                 interval = c(n + f0/d, n + f0*d),
                 type = "parametric",
                 model = "Poisson (homogeneous)",
                 frequentist = TRUE,
                 parametric = TRUE,
                 reasonable = FALSE,
                 interval_type = "Approximate: log-normal")
  
}
