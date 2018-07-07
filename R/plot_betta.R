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
#' @param myy Deprecated, for backwards compatibility
#' @param mymain Deprecated, for backwards compatibility
#' @param mycol Deprecated, for backwards compatibility
#' @param labs Deprecated, for backwards compatibility
#' @param mypch Deprecated, for backwards compatibility
#' @param myxlim Deprecated, for backwards compatibility
#' @author Amy Willis
#' @seealso \code{\link{betta}}
#' @references Willis, A., Bunge, J., and Whitman, T. (2015). Inference for
#' changes in biodiversity. arXiv preprint.
#' @keywords diversity
#' @import ggplot2
#' @examples
#' betta_pic(c(1552, 1500, 884), c(305, 675, 205), mymain = "Example title")
#' 
#' @export betta_pic
betta_pic <- function(y, se, 
                      x = 1:length(y), 
                      ylimu=NULL, myy=NULL, 
                      mymain=NULL, 
                      mycol=NULL,
                      labs=NULL, mypch=NULL,
                      myxlim=NULL) {
  
  df <- data.frame("x" = x, "y" = y, "lower" = y - 1.96*se, "upper" = y + 1.96*se)
  
  ggplot(df, aes(x = x)) +
    ylab("Alpha-diversity estimate") +
    geom_point(aes(y = y)) +
    geom_linerange(aes(ymin = lower, ymax = upper)) +
    theme_bw()
    
  
  # n <- length(y)
  # ylimu <- ifelse(is.na(ylimu), max(y+2*se, na.rm = T), ylimu)
  # par(xpd=NA)
  # plot(0, 0, type="n", xlim=myxlim, ylim=c(0, ylimu), xlab="", bty="n", ylab=myy, main=mymain)
  # for (i in 1:n) {
  #   if(!is.na(y[i]) & !is.na(x[i])) {
  #     points(x[i], y[i], pch=mypch[i], col=mycol[i])
  #     lines(c(x[i], x[i]), c(max(0, y[i]-1.96*se[i], na.rm = T), min(y[i]+1.96*se[i], ylimu)), col=mycol[i])
  #   }
  # }
  # if(!is.na(labs)) axis(1, at=1:length(y), labels=labs, las=2, cex=0.8)
}
