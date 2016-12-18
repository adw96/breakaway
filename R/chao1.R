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
