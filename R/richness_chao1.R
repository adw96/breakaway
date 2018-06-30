#' Chao1 species richness estimator
#' 
#' This function implements the Chao1 richness estimate, which is often
#' mistakenly referred to as an index.
#' 
#' 
#' @param input_data The sample frequency count table for the population of interest.
#' See dataset apples for sample formatting.
#' @param output Logical: whether the results should be printed to screen.
#' @param answers Should the answers be returned as a list?
#' @return The results of the estimator, including standard error.
#' @note The authors of this package strongly discourage the use of this
#' estimator.  It is only valid when you wish to assume that every taxa has
#' equal probability of being observed. You don't really think that's possible,
#' do you?
#' @author Amy Willis
#' @examples
#' 
#' 
#' chao1(apples)
#' 
#' 
#' @export chao1
chao1 <- function(input_data, output=TRUE, answers=FALSE) {
  
  my_data <- convert(input_data)
  
  index  <- 1:max(my_data[,1])
  frequency_index <- rep(0, length(index))
  frequency_index[my_data[,1]] <- my_data[,2]
  f1  <- frequency_index[1]
  f2 <- frequency_index[2]
  n <- sum(frequency_index)
  
  f0 <- f1^2/(2*f2)
  diversity <- n + f0
  
  diversity_se <- sqrt(f2*(0.5*(f1/f2)^2 + (f1/f2)^3 + 0.25*(f1/f2)^4))
  
  if(output) {
    cat("################## Chao1 ##################\n")
    cat("\tThe estimate of total diversity is", round(diversity),
        "\n \t with std error",round(diversity_se),"\n")
    cat("You know that this estimate is only valid if all taxa are equally abundant, right?\n")
  }
  if(answers) {
    result <- list()
    result$name <- "Chao1"
    result$est <- diversity
    result$seest <- diversity_se
    d <- exp(1.96*sqrt(log(1+result$seest^2/f0)))
    result$ci <- c(n+f0/d,n+f0*d)
    return(result)
  }
}

