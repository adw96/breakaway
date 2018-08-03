
#' The untransformed weighted linear regression estimator for species richness estimation
#' 
#' This function implements the untransformed version of the species richness
#' estimation procedure outlined in Rocchetti, Bunge and Bohning (2011).
#' 
#' @param input_data An input type that can be processed by \code{convert()} or a \code{phyloseq} object
#' @param cutoff Maximum frequency count to use
#' @param print Deprecated; only for backwards compatibility
#' @param plot Deprecated; only for backwards compatibility
#' @param answers Deprecated; only for backwards compatibility
#'
#' @return An object of class \code{alpha_estimate}, or \code{alpha_estimates} for \code{phyloseq} objects
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
#' @export
wlrm_untransformed  <- function(input_data, 
                                cutoff=NA,
                                print=NULL, 
                                plot=NULL, 
                                answers=NULL) {
  
  if (class(input_data) == "phyloseq") {
    if (input_data %>% otu_table %>% taxa_are_rows) {
      return(input_data %>% 
               get_taxa %>%
               apply(2, function(x) wlrm_untransformed(make_frequency_count_table(x))) %>%
               alpha_estimates)
    } else {
      return(input_data %>% 
               otu_table %>%
               apply(1, function(x) wlrm_untransformed(make_frequency_count_table(x))) %>%
               alpha_estimates)
    }
  }
  
  my_data <- convert(input_data)
  
  n <- sum(my_data[,2])
  
  if (my_data[1,1] != 1 || my_data[1,2] == 0) {
    wlrm_alpha_estimate <- alpha_estimate(estimand = "richness",
                                          estimate = NA,
                                          error = NA,
                                          name = "wlrm_untransformed",
                                          interval = c(NA, NA),
                                          type = "parametric",
                                          model = "Negative Binomial",
                                          frequentist = TRUE, 
                                          parametric = TRUE,
                                          reasonable = FALSE,
                                          warnings = "no singleton count")
    
  } else {
    
    # TODO fix 
    f1 <- my_data[1,2]
    
    if (is.na(cutoff)) {
      iss <- my_data$index
      fis <- my_data$frequency
      length_fis <- length(fis)
      
      breaks <- which(iss[-1] - iss[-length_fis] > 1)
      cutoff <- ifelse(is.na(breaks[1]), length_fis, breaks[1])
    }
    
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
      
      my_data <- my_data[1:cutoff, ]
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
        
        if (is.nan(f0_se)) {
          my_warning <- "infinite std error"
          diversity_se <- NaN
          d <- NaN
        } else {
          diversity_se <- sqrt(f0_se^2+n*f0/(n+f0))
          d <- exp(1.96*sqrt(log(1+diversity_se^2/f0)))
          my_warning <- NULL
        }
        
        yhats <- fitted(lm2)/(my_data[1:(cutoff-1),1]+1)
        
        plot_data <- rbind(data.frame("x" = xs, "y" = lhs$y, 
                                      "type" = "Observed"),
                           data.frame("x" = xs, "y" = yhats, 
                                      "type" = "Fitted"),
                           data.frame("x" = 0, "y" = b0_hat, 
                                      "type" = "Prediction"))
        
        my_plot <- ggplot(plot_data, 
                          aes_string(x = "x", 
                                     y = "y",
                                     col = "type", 
                                     pch = "type")) +
          geom_point() +
          labs(x = "x", y = "f(x+1)/f(x)", title = "Plot of ratios and fitted values: tWLRLM") +
          theme_bw()
        
      } else {
        diversity <- NA
        diversity_se <- NA
        my_warning <- "negative richness estimate"
        my_plot <- NULL
        f0 <- NA; d <- NA
        # if(print) cat("The utWLRM failed to produce an estimate.\n")
      }
      wlrm_alpha_estimate <- alpha_estimate(estimate = diversity,
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
                                                         cutoff = cutoff))
      
    }
  }
  wlrm_alpha_estimate
}