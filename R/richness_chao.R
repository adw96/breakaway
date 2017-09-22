#' Chao-Bunge species richness estimator
#' 
#' This function implements the species richness estimation procedure outlined
#' in Chao & Bunge (2002).
#' 
#' 
#' @param data The sample frequency count table for the population of interest.
#' See dataset apples for sample formatting.
#' @param cutoff The maximum frequency to use in fitting.
#' @param print Logical: whether the results should be printed to screen.
#' @param answers Should the answers be returned as a list?
#' @return The results of the estimator, including standard error.
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
#' chao_bunge(apples)
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
#' @export chao_bunge
chao_bunge <- function(data, cutoff=10, print=TRUE, answers=FALSE) {
  
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
  n <- sum(frequency_index)
  
  d_a <- sum(data[data[,1]>cutoff,2])
  k <- 2:cutoff
  m <- 1:cutoff
  numerator <- frequency_index[k]
  denominator <- 1 - f1*sum(m^ 2*frequency_index[m])/(sum(m*frequency_index[m]))^ 2
  diversity  <- d_a + sum(numerator /denominator)
  f0 <- diversity - sum(frequency_index)
  
  ## standard error is to be made available in version 3.1; please contact the author with requests
  diversity_se <- NA
  
  if(print) {
    cat("################## Chao-Bunge ##################\n")
    cat("\tThe estimate of total diversity is", round(diversity),
        "\n \t with std error",round(diversity_se),"\n")
  }
  if(answers) {
    result <- list()
    result$name <- "Chao-Bunge"
    result$est <- diversity
    result$seest <- as.vector(diversity_se)
    d <- exp(1.96*sqrt(log(1+result$seest^2/f0)))
    result$ci <- c(n+f0/d,n+f0*d)
    return(result)
  }
}





























#' Chao1 species richness estimator
#' 
#' This function implements the Chao1 richness estimate, which is often
#' mistakenly referred to as an index.
#' 
#' 
#' @param data The sample frequency count table for the population of interest.
#' See dataset apples for sample formatting.
#' @param print Logical: whether the results should be printed to screen.
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
#' chao1(apples)
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
#' @export chao1
chao1 <- function(data, print=TRUE, answers=FALSE) {

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

  f0 <- f1^2/(2*f2)
  diversity <- n + f0

  diversity_se <- sqrt(f2*(0.5*(f1/f2)^2 + (f1/f2)^3 + 0.25*(f1/f2)^4))

  if(print) {
    cat("################## Chao1 ##################\n")
    cat("\tThe estimate of total diversity is", round(diversity),
        "\n \t with std error",round(diversity_se),"\n")
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





























#' Bias-corrected Chao1 species richness estimator
#' 
#' This function implements the bias-corrected Chao1 richness estimate.
#' 
#' 
#' @param data The sample frequency count table for the population of interest.
#' See dataset apples for sample formatting.
#' @param print Logical: whether the results should be printed to screen.
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
#' @export chao1_bc
chao1_bc <- function(data, print=TRUE, answers=FALSE) {
  
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
  
  if(print) {
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
