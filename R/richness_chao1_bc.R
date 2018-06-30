#' Bias-corrected Chao1 species richness estimator
#' 
#' This function implements the bias-corrected Chao1 richness estimate.
#' 
#' 
#' @param data The sample frequency count table for the population of interest.
#' See dataset apples for sample formatting.
#' @param output Logical: whether the results should be printed to screen.
#' @param answers Should the answers be returned as a list?
#' @return The results of the estimator, including standard error.
#' @note The authors of this package strongly discourage the use of this
#' estimator. It is underpinned by totally implausible assumptions that are not
#' made by other richness estimators.  Bias correcting Chao1 is the least of
#' your problems.
#' @author Amy Willis
#' @examples
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' chao1_bc(apples)
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' @export chao1_bc
chao1_bc <- function(data, output=TRUE, answers=FALSE) {
  
  if( !(is.matrix(data) || is.data.frame(data))) {
    filename <- data
    ext <- substr(filename, nchar(filename)-2, nchar(filename))
    if (ext == "csv") {
      data <- read.table(file=filename, header=0,sep=",")
      if( data[1,1] !=1) data <- read.table(filename, header=1,sep=",")
    } else if (ext == "txt") {
      data <- read.table(file=filename, header=0)
    } else cat("Please input your data as a txt or csv file,
               or as an R dataframe or matrix.")
  }
  
  if ( is.factor(data[,1]) ) {
    fs <- as.numeric(as.character(data[,1]))
    data <- cbind(fs,data[,2])
    data <- data[data[,1]!=0,]
  }
  
  
  
  index  <- 1:max(data[,1])
  frequency_index <- rep(0, length(index))
  frequency_index[data[,1]] <- data[,2]
  f1  <- frequency_index[1]
  f2 <- frequency_index[2]
  n <- sum(frequency_index)
  
  f0 <- f1*(f1-1)/(2*(f2+1))
  diversity <- n + f0
  
  diversity_se <- sqrt(f1*(f1-1)/(2*(f2+1)) + f1*(2*f1-1)^2/(4*(f2+1)^2) + f1^2*f2*(f1-1)^2/(4*(f2+1)^4))
  
  if(output) {
    cat("################## Bias-corrected Chao1 ##################\n")
    cat("\tThe estimate of total diversity is", round(diversity),
        "\n \t with std error",round(diversity_se),"\n")
  }
  if(answers) {
    result <- list()
    result$name <- "Chao1_bc"
    result$est <- diversity
    result$seest <- diversity_se
    d <- exp(1.96*sqrt(log(1+result$seest^2/f0)))
    result$ci <- c(n+f0/d,n+f0*d)
    return(result)
  }
}
