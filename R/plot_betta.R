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
#' 
#' @return A ggplot object.
#' 
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
  
  ggplot(df, aes_string(x = "x")) +
    ylab("Alpha-diversity estimate") +
    geom_point(aes_string(y = "y")) +
    geom_linerange(aes_string(ymin = "lower", 
                              ymax = "upper")) +
    theme_bw()
    
}
