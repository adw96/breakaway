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
