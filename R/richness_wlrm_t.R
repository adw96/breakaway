#' The transformed weighted linear regression estimator for species richness estimation
#' 
#' This function implements the transformed version of the species richness
#' estimation procedure outlined in Rocchetti, Bunge and Bohning (2011).
#' 
#' @param input_data An input type that can be processed by \code{convert()} or a \code{phyloseq} object
#' @param cutoff Maximum frequency count to use
#' @param print Deprecated; only for backwards compatibility
#' @param plot Deprecated; only for backwards compatibility
#' @param answers Deprecated; only for backwards compatibility
#'
#' @return An object of class \code{alpha_estimate}, or \code{alpha_estimates} for \code{phyloseq} objects
#' 
#' @note While robust to many different structures, model is almost always
#' misspecified. The result is usually implausible diversity estimates with
#' artificially low standard errors. Extreme caution is advised.
#' @author Amy Willis
#' @seealso \code{\link{breakaway}}; \code{\link{apples}};
#' \code{\link{wlrm_untransformed}}
#' @references Rocchetti, I., Bunge, J. and Bohning, D. (2011). Population size
#' estimation based upon ratios of recapture probabilities. \emph{Annals of
#' Applied Statistics}, \bold{5}.
#' @keywords diversity models
#' 
#' @import ggplot2 
#' @import stats 
#' 
#' @examples
#' 
#' wlrm_transformed(apples)
#' wlrm_transformed(apples, plot = FALSE, print = FALSE, answers = TRUE)
#' 
#' @export
wlrm_transformed <- function(input_data, 
                             cutoff = NA,
                             print = NULL, 
                             plot = NULL, 
                             answers = NULL) {
  
  
  if ((intersect(class(input_data), 
                 c("phyloseq", "otu_table")) %>% length) > 0) {
    return(physeq_wrap(fn = wlrm_transformed, physeq = input_data,
                       cutoff, print, plot, answers))
  }
  
  my_data <- convert(input_data)
  
  if (my_data[1,1] != 1 || my_data[1,2] == 0) {
    wlrm_alpha_estimate <- alpha_estimate(estimand = "richness",
                                          estimate = NA,
                                          error = NA,
                                          name = "wlrm_transformed",
                                          interval = c(NA, NA),
                                          type = "parametric",
                                          model = "Negative Binomial",
                                          frequentist = TRUE, 
                                          parametric = TRUE,
                                          reasonable = FALSE,
                                          warnings = "no singleton count")
    
  } else {
    
    n <- sum(my_data[,2])
    # TODO fix 
    f1 <- my_data[1,2]
    
    cutoff <- cutoff_wrap(my_data, requested = cutoff) 
    
    if (cutoff < 4) {
      wlrm_alpha_estimate <- alpha_estimate(estimand = "richness",
                                            estimate = NA,
                                            error = NA,
                                            model = "Negative Binomial",
                                            name = "wlrm_transformed",
                                            frequentist = TRUE,
                                            interval = c(NA, NA),
                                            parametric = TRUE,
                                            reasonable = FALSE,
                                            warnings = "insufficient contiguous frequencies")
    } else {
      
      my_data <- my_data[1:cutoff,]
      ys <- (my_data[1:(cutoff-1),1]+1)*my_data[2:cutoff,2]/my_data[1:(cutoff-1),2]
      xs <- 1:(cutoff-1)
      xbar <- mean(xs)
      lhs <- list("x"=xs-xbar,
                  "y"=my_data[2:cutoff,2]/my_data[1:(cutoff-1),2])
      
      weights_trans <- 1 / (1 / my_data[-1,2] + 1 / my_data[-cutoff,2])
      lm1 <- lm(log(ys) ~ xs,
                weights = weights_trans)
      b0_hat <- summary(lm1)$coef[1,1]
      b0_se <- summary(lm1)$coef[1,2]
      f0 <- f1*exp(-b0_hat)
      diversity <- f0 + n
      f0_se <- sqrt( (exp(-b0_hat))^2*f1*(b0_se^2*f1+1) )
      
      if (is.nan(f0_se)) {
        my_warning <- "infinite std error"
        diversity_se <- NaN
        d <- NaN
      } else {
        diversity_se <- sqrt(f0_se^2+n*f0/(n+f0))
        d <- exp(1.96*sqrt(log(1+diversity_se^2/f0)))
        my_warning <- NULL
      }
      
      yhats <- exp(fitted(lm1))/(my_data[1:(cutoff-1),1]+1)
      
      plot_data <- rbind(data.frame("x" = xs, "y" = lhs$y, 
                                    "type" = "Observed"),
                         data.frame("x" = xs, "y" = yhats, 
                                    "type" = "Fitted"),
                         data.frame("x" = 0, "y" = exp(b0_hat), 
                                    "type" = "Prediction"))
      
      my_plot <- ggplot(plot_data, 
                        aes(x = x, y = y, col = type, pch = type)) +
        geom_point() +
        labs(x = "x", y = "f(x+1)/f(x)", title = "Plot of ratios and fitted values: tWLRLM") +
        theme_bw()
      
      wlrm_alpha_estimate <- alpha_estimate(estimate = diversity,
                                            error = diversity_se,
                                            estimand = "richness",
                                            name = "wlrm_transformed",
                                            interval = c(n + f0/d, n + f0*d),
                                            type = "parametric",
                                            model = "Negative Binomial",
                                            frequentist = TRUE,
                                            parametric = TRUE,
                                            reasonable = FALSE,
                                            interval_type = "Approximate: log-normal",
                                            plot = my_plot, 
                                            warnings = my_warning,
                                            other = list(para = summary(lm1)$coef[,1:2],
                                                         full = lm1,
                                                         cutoff = cutoff))
    }
    
  }
  wlrm_alpha_estimate
}