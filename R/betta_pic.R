#' function for plotting total diversity
#' 
#' A simple plotting interface for comparing total diversity across samples or
#' a covariate gradient.
#' 
#' 
#' @param y A vector of estimates of total diversity. Other parameter estimates
#' are accessible; this method may be used for plotting any parameter
#' estimates..
#' @param se The standard errors in \samp{y}, the diversity (or other
#' parameter's) estimates.
#' @param x A vector of covariates to form the x-coordinates of the intervals.
#' If no argument is given, defaults to the order.
#' @param ylimu The upper endpoint of the y-axis.
#' @param myy A label for the y-axis. Default is none.
#' @param mymain A main title for the plot. Default is none.
#' @param mycol Colors for the plotting points. Default is black.
#' @param labs x-axis labels. Default is none (non x-axis plotted).
#' @param mypch Plotting characters for the estimates. Defaults to circles.
#' @param myxlim A vector of x-axis limits. Default is 0 to 10\% greater than
#' the maximum x value.
#' @author Amy Willis
#' @seealso \code{\link{betta}}
#' @references Willis, A., Bunge, J., and Whitman, T. (2015). Inference for
#' changes in biodiversity. arXiv preprint.
#' @keywords diversity
#' @examples
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' betta_pic(c(1552, 1500, 884), c(305, 675, 205), mymain = "Example title")
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' @export betta_pic
betta_pic <- function(y, se, x=1:length(y), ylimu=NA, myy=NA, mymain=NA, mycol=rep("black", length(y)), labs=NA, mypch=rep(16, length(y)), myxlim=c(0.8*min(x, na.rm=T), 1.2*max(x, na.rm=T))) {
  n <- length(y)
  ylimu <- ifelse(is.na(ylimu), max(y+2*se, na.rm = T), ylimu)
  par(xpd=NA)
  plot(0, 0, type="n", xlim=myxlim, ylim=c(0, ylimu), xlab="", bty="n", ylab=myy, main=mymain)
  for (i in 1:n) {
    if(!is.na(y[i]) & !is.na(x[i])) {
      points(x[i], y[i], pch=mypch[i], col=mycol[i])
      lines(c(x[i], x[i]), c(max(0, y[i]-1.96*se[i], na.rm = T), min(y[i]+1.96*se[i], ylimu)), col=mycol[i])
    }
  }
  if(!is.na(labs)) axis(1, at=1:length(y), labels=labs, las=2, cex=0.8)
}


#' @export hill_pic
hill_pic <- function(input, res = 20, method = "resample") {
  
  xx <- seq(from = 0, to = 4, length.out = res)
  my_results <- alpha_better(input = input, xx)
  raw <- hill(input = input, q = xx)
  ests <- my_results$Estimate
  
  if (method  == "delta") {
    ses <- my_results$StdError
    upper <- ests + 2 * ses
    lower <- ests - 2 * ses
  } else {
    if (method != "resample") warning("Using resampled standard errors")
    my_sample_size <- sum(input[, 1] * input[,2])
    y1 <- replicate(n = 20, 
                    subsample_otu(input), 
                    simplify = "list")
    y2 <- apply(y1, 2, function(x) data.frame("Index"=as.numeric(as.character(unlist(x[1]))), "Freq"=as.numeric(as.character(unlist(x[2])))))
    y3 <- lapply(y2, alpha_better, q = xx)
    y4 <- lapply(y3, function(x) x[,2])
    y5 <- do.call(rbind, y4)
    ests <- apply(y5, 2, function(x) quantile(x, 0.5))
    lower <- apply(y5, 2, function(x) quantile(x, 0.025))
    upper <- apply(y5, 2, function(x) quantile(x, 0.975))
  }
  
  plot(c(xx, xx), c(lower, upper), type = "n", xlab = "q", ylab = "Hill Number estimate")
  points(xx, ests, pch = 16)
  lines(xx, lower, lty = 2)
  lines(xx, upper, lty = 2)
  points(xx, raw, pch = 16, col = "red")
}


