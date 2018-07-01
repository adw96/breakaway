
#' The untransformed weighted linear regression estimator for species richness estimation
#' 
#' This function implements the untransformed version of the species richness
#' estimation procedure outlined in Rocchetti, Bunge and Bohning (2011).
#' 
#' @param input_data An input type that can be processed by \code{convert()}
#' @param cutoff Maximum frequency count to use
#' @param print Deprecated; only for backwards compatibility
#' @param plot Deprecated; only for backwards compatibility
#' @param answers Deprecated; only for backwards compatibility
#' @return An object of class \code{alpha_estimate}
#' 
#' @note This estimator is based on the negative binomial model and for that
#' reason generally produces poor fits to microbial data. The result is usually
#' artificially low standard errors. Caution is advised.
#' 
#' @author Amy Willis
#' @seealso \code{\link{breakaway}}; \code{\link{apples}};
#' \code{\link{wlrm_transformed}}
#' @references Rocchetti, I., Bunge, J. and Bohning, D. (2011). Population size
#' estimation based upon ratios of recapture probabilities. \emph{Annals of
#' Applied Statistics}, \bold{5}.
#' @keywords diversity models
#' @importFrom graphics legend points
#' @importFrom utils head read.table
#' 
#' @import ggplot2 
#' 
#' @examples
#' 
#' wlrm_untransformed(apples)
#' 
#' @export wlrm_untransformed
wlrm_untransformed  <- function(input_data, 
                                cutoff=NA,
                                print=NULL, 
                                plot=NULL, 
                                answers=NULL) {
  my_data <- convert(input_data)
  
  if (my_data[1,1] != 1 || my_data[1,2] == 0) {
    diversity <- NA
    diversity_se <- NA
    my_warning <- "no singletons"
  } else {
    
    n <- sum(my_data[,2])
    # TODO fix 
    f1 <- my_data[1,2]
    
    if (is.na(cutoff)) {
      cutoff <- ifelse(is.na(which(my_data[-1,1]-my_data[-length(my_data[,1]),1]>1)[1]),length(my_data[,1]),which(my_data[-1,1]-my_data[-length(my_data[,1]),1]>1)[1])
    }
    
    my_data <- my_data[1:cutoff,]
    ys <- (my_data[1:(cutoff-1),1]+1)*my_data[2:cutoff,2]/my_data[1:(cutoff-1),2]
    xs <- 1:(cutoff-1)
    xbar <- mean(xs)
    lhs <- list("x"=xs-xbar,
                "y"=my_data[2:cutoff, 2] / my_data[1:(cutoff-1), 2])
    
    weights_untrans <- 1/(my_data[-1, 1]^2 * my_data[-1, 2] * my_data[-cutoff, 2]^-2 * 
                            (my_data[-1, 2] * my_data[-cutoff,2]^-1 + 1))
    
    lm2 <- lm(ys ~ xs, 
              weights = weights_untrans)
    b0_hat <- summary(lm2)$coef[1,1]
    
    if(b0_hat > 0) {
      
      b0_se <- summary(lm2)$coef[1,2]
      f0 <- f1/b0_hat
      diversity <- f0 + n
      f0_se <- sqrt( f1*(1-f1/n)*b0_hat^-2 + f1^2*b0_hat^-4*b0_se^2   ) #1st order d.m.
      diversity_se <- sqrt(f0_se^2+n*f0/(n+f0))
      d <- exp(1.96*sqrt(log(1+diversity_se^2/f0)))
      
      
      
      # if(print) {
      #   cat("################## utWLRM ##################\n")
      #   cat("\tThe best estimate of total diversity is", round(diversity),
      #       "\n \t with std error",round(diversity_se),"\n")
      # }
      # if(answers) {
      #   result <- list()
      #   result$name <- "utWLRM"
      #   result$para <- summary(lm2)$coef[,1:2]
      #   result$est <- diversity
      #   result$seest <- as.vector(diversity_se)
      #   result$full <- lm2
      #   result$ci <- c(n+f0/d,n+f0*d)
      #   return(result)
      # }
      
      yhats <- fitted(lm2)/(my_data[1:(cutoff-1),1]+1)
      
      plot_data <- rbind(data.frame("x" = xs, "y" = lhs$y, 
                                    "type" = "Observed"),
                         data.frame("x" = xs, "y" = yhats, 
                                    "type" = "Fitted"),
                         data.frame("x" = 0, "y" = b0_hat, 
                                    "type" = "Prediction"))
      
      my_plot <- ggplot(plot_data, 
                        aes(x = x, 
                            y = y,
                            col = type, pch = type)) +
        geom_point() +
        labs(x = "x", y = "f(x+1)/f(x)", title = "Plot of ratios and fitted values: tWLRLM") +
        theme_bw()
      
      # if(plot) {
      #   yhats <- fitted(lm2)/(my_data[1:(cutoff-1),1]+1)
      #   par(new=FALSE)
      #   plot(xs,lhs$y,xlim=c(0,max(xs)+1),ylim=c(min(lhs$y,yhats),max(lhs$y)*1.05),
      #        ylab="f(x+1)/f(x)",xlab="x",main="Plot of ratios and fitted values");
      #   points(xs,yhats,pch=18)
      #   points(0,b0_hat,col="red",pch=18)
      #   legend(0,max(lhs$y),c("Fitted values", "Prediction"),pch=c(18,18),col=c("black", "red"),cex=0.8,bty = "n")
      # }
      my_warning <- NULL
      
    } else {
      diversity <- NA
      diversity_se <- NA
      my_warning <- "negative richness estimate"
      my_plot <- NULL
      f0 <- NA; d <- NA
      # if(print) cat("The utWLRM failed to produce an estimate.\n")
    }
  }
  alpha_estimate(estimate = diversity,
                 error = diversity_se,
                 estimand = "richness",
                 name = "wlrm_untransformed",
                 interval = c(n + f0/d, n + f0*d),
                 type = "parametric",
                 model = "Negative Binomial",
                 frequentist = TRUE,
                 parametric = TRUE,
                 reasonable = FALSE,
                 interval_type = "Approximate: log-normal",
                 plot = my_plot, 
                 warnings = my_warning,
                 other = list(para = summary(lm2)$coef[,1:2],
                              full = lm2,
                              cutoff = cutoff),
                 est = diversity,
                 seest = diversity_se)
}