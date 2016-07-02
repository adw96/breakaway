wlrm_untransformed  <- function(data, print=TRUE, plot=FALSE, answers=FALSE, cutoff=NA) {

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

  if(length(data)>1) {
    if (data[1,1]!=1 || data[1,2]==0) {
      if(print) cat("You don't have an observed singleton count.\n The WLRM isn't built for that data structure.\n Consider using breakaway_nof1 instead.\n")
    } else {
      data <- data[!(data[,2]==0 | is.na(data[,2])),]
      orig_data <- data
      n <- sum(orig_data[,2])
      f1 <- data[1,2]

      if (is.na(cutoff)) {
        cutoff <- ifelse(is.na(which(data[-1,1]-data[-length(data[,1]),1]>1)[1]),length(data[,1]),which(data[-1,1]-data[-length(data[,1]),1]>1)[1])
      }
      data <- data[1:cutoff,]
      ys <- (data[1:(cutoff-1),1]+1)*data[2:cutoff,2]/data[1:(cutoff-1),2]
      xs <- 1:(cutoff-1)
      xbar <- mean(xs)
      lhs <- list("x"=xs-xbar,"y"=data[2:cutoff,2]/data[1:(cutoff-1),2])

      result <- list()

      weights_untrans <- (data[-1,1]^2*data[-1,2]*data[-cutoff,2]^-2*(data[-1,2]*data[-cutoff,2]^-1+1))^-1

      lm2 <- lm(ys~xs,weights=weights_untrans)
      b0_hat <- summary(lm2)$coef[1,1]

      if(b0_hat > 0) {
        b0_se <- summary(lm2)$coef[1,2]
        f0 <- f1/b0_hat
        diversity <- f0 + n
        f0_se <- sqrt( f1*(1-f1/n)*b0_hat^-2 + f1^2*b0_hat^-4*b0_se^2   ) #1st order d.m.
        diversity_se <- sqrt(f0_se^2+n*f0/(n+f0))

        if(print) {
          cat("################## utWLRM ##################\n")
          cat("\tThe best estimate of total diversity is", round(diversity),
              "\n \t with std error",round(diversity_se),"\n")
        }
        if(answers) {
          result$name <- "utWLRM"
          result$para <- summary(lm2)$coef[,1:2]
          result$est <- diversity
          result$seest <- as.vector(diversity_se)
          result$full <- lm2
          d <- exp(1.96*sqrt(log(1+result$seest^2/f0)))
          result$ci <- c(n+f0/d,n+f0*d)
          return(result)
        }

        if(plot) {
          yhats <- fitted(lm2)/(data[1:(cutoff-1),1]+1)
          par(new=FALSE)
          plot(xs,lhs$y,xlim=c(0,max(xs)+1),ylim=c(min(lhs$y,yhats),max(lhs$y)*1.05),
               ylab="f(x+1)/f(x)",xlab="x",main="Plot of ratios and fitted values");
          points(xs,yhats,pch=18)
          points(0,b0_hat,col="red",pch=18)
          legend(0,max(lhs$y),c("Fitted values", "Prediction"),pch=c(18,18),col=c("black", "red"),cex=0.8,bty = "n")
        }
      } else {
        if(print) cat("The utWLRM failed to produce an estimate.\n")
      }
    }
  }
}
