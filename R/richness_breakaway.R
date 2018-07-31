#' Species richness estimation with breakaway
#' 
#' breakaway is a wrapper for modern species richness estimation
#' for modern datasets
#' 
#' 
#' @param input_data An input type that can be processed by \code{convert()}
#' @param output Deprecated; only for backwards compatibility
#' @param answers Deprecated; only for backwards compatibility
#' @param plot Deprecated; only for backwards compatibility
#' @param print Deprecated; only for backwards compatibility
#' @return An object of class \code{alpha_estimate}
#' @note \samp{breakaway} presents an estimator of species richness that is
#' well-suited to the high-diversity/microbial setting. However, many microbial
#' datasets display more diversity than the Kemp-type models can permit. In
#' this case, the log-transformed WLRM diversity estimator of Rocchetti et. al.
#' (2011) is returned. The authors' experience suggests that some datasets that
#' require the log-transformed WLRM contain ``false'' diversity, that is,
#' diversity attributable to sequencing errors (via an inflated singleton
#' count). The authors encourage judicious use of diversity estimators when the
#' dataset may contain these errors, and recommend the use of
#' \code{\link{breakaway_nof1}} as an exploratory tool in this case.
#' 
#' @author Amy Willis
#' 
#' @import magrittr
#' @import ggplot2
#' @import phyloseq
#' 
#' @seealso \code{\link{breakaway_nof1}}; \code{\link{kemp}}; \code{\link{apples}}
#' 
#' @references Willis, A. and Bunge, J. (2015). Estimating diversity via
#' frequency ratios. \emph{Biometrics}, \bold{71}(4), 1042--1049.
#' 
#' @keywords diversity microbial models nonlinear
#' @examples
#' 
#' breakaway(apples)
#' breakaway(apples, plot = FALSE, output = FALSE, answers = TRUE)
#' 
#' @export 
breakaway <- function(input_data, 
                      output = NULL, plot = NULL, 
                      answers = NULL, print = NULL) {
  UseMethod("breakaway", input_data)
}

#' @export
breakaway.phyloseq <- function(input_data, 
                               output = NULL, plot = NULL, 
                               answers = NULL, print = NULL) {
  
  if (input_data %>% otu_table %>% taxa_are_rows) {
    input_data %>% 
      otu_table %>%
      apply(2, breakaway) %>%
      alpha_estimates
  } else {
    input_data %>% 
      otu_table %>%
      apply(1, breakaway) %>%
      alpha_estimates
  }
}

#' @export
breakaway.matrix <- function(input_data, 
                             output = NULL, plot = NULL, 
                             answers = NULL, print = NULL) {
  
  input_data %>%
    as.data.frame %>%
    breakaway(output, plot, answers, print)
  
}

#' @export
breakaway.data.frame <- function(input_data, 
                                 output = NULL, plot = NULL, 
                                 answers = NULL, print = NULL) {
  
  ## if a frequency count matrix...
  if (dim(input_data)[2] == 2 & 
      all(input_data[,1] %>% sort == input_data[,1])) {
    breakaway.default(input_data)
  } else {
    
    warning("Assuming taxa are rows")
    
    input_data %>% 
      apply(2, breakaway) %>%
      alpha_estimates
    
  }
}

#' @export
breakaway.default <- function(input_data, 
                              output = NULL, plot = NULL, 
                              answers = NULL, print = NULL) {
  
  my_data <- convert(input_data)
  
  breakaway_alpha_estimate <- kemp(input_data)
  
  if (!is.null(breakaway_alpha_estimate$warnings) |
      breakaway_alpha_estimate$error >  breakaway_alpha_estimate$estimate) {
      
    breakaway_alpha_estimate <- wlrm_untransformed(input_data)
    
    if (!is.null(breakaway_alpha_estimate$warnings)) {
      
      breakaway_alpha_estimate <- chao_bunge(input_data)
      
      if (!is.null(breakaway_alpha_estimate$warnings)) {
        
        breakaway_alpha_estimate <- wlrm_transformed(input_data)
        
        if (!is.null(breakaway_alpha_estimate$warnings)) {
          
          breakaway_alpha_estimate <- wlrm_transformed(input_data)
          
          if (!is.null(breakaway_alpha_estimate$warnings)) {
            
            breakaway_alpha_estimate <- chao1(input_data)
            
            if (!is.null(breakaway_alpha_estimate$warnings)) {
              
              warning("Could not estimate missing taxa; outputting sample richness")
              breakaway_alpha_estimate <- sample_richness(input_data)
              
            }
            

          }
          
        }
        
      }
    }
    
  }
  breakaway_alpha_estimate$name <- "breakaway"
  breakaway_alpha_estimate
  
}
