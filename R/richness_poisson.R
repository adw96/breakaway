#' PoissonModel
#' 
#' A model to estimate the number of missing taxa under a Poisson Model
#' 
#' @param input_data A frequency count table
#' @param cutoff The largest frequency to use for predicting f0. Defaults to 10.
#' 
#' @return An object of class \code{alpha_estimate}.
#' 
#' @import magrittr 
#' @importFrom stats uniroot
#' 
#' @export
poisson_model <- function(input_data, 
                          cutoff = 10) {
  
  my_warnings <- NULL
  
  if ((intersect(class(input_data), 
                 c("phyloseq", "otu_table")) %>% length) > 0) {
    return(physeq_wrap(fn = poisson_model, physeq = input_data,
                       cutoff))
  }
  
  
  my_data <- convert(input_data)
  
  if (is.na(cutoff)) cutoff <- 10
  # cutoff <- cutoff_wrap(my_data, requested = cutoff) 
  
  input_data <- my_data
  
  # print(input_data)
  included <- input_data[input_data$index <= cutoff, ]
  excluded <- input_data[input_data$index > cutoff, ]
  
  if (nrow(included) == 0) {
    included <- input_data
    excluded <- list("index" = Inf, "frequency" = 0)
    warning("Cut-off was too low: no data available for estimation")
  }
  
  # s = sum f_j
  cc <- included$frequency %>% sum
  cc_excluded  <- excluded$frequency %>% sum
  
  # n = sum j f_j
  nn <- crossprod(included$index, included$frequency)
  
  ## MLE for lambda is solution to 
  ## (1-exp(-lambda))/lambda = c/n
  poisson_fn <- function(lambda) {
    (1-exp(-lambda))/lambda - cc/nn
  }
  
  
  if (cc == nn) {
    
    warning("Only one frequency count was observed below cutoff.\nConsider increasing cutoff.")
    
    lambda_hat <- 0 # maximiser of poisson_fn is zero
    ccc_hat <- cc + cc_excluded
    ccc_se <- 0
    f0 <- 0
    d <- 1
    
  } else {
    # eqn nonlinear; so use Newton's method
    # for large lambda, lambda/(1-exp(-lambda)) \approx lambda, so choose upper bound above n/c
    # for large lambda, lambda/(1-exp(-lambda)) \approx lambda, so choose lower bound below 
    
    lambda_hat <- uniroot(poisson_fn, c(cc / (100*nn), nn/cc + 100))$root
    
    ccc_subset <- cc / (1-exp(-lambda_hat)) 
    ccc_hat <- ccc_subset + cc_excluded
    
    ccc_se <- sqrt(ccc_subset/(exp(lambda_hat)-1-lambda_hat))
    
    f0 <- ccc_subset - cc # previously ccc_hat - (cc + cc_excluded)
    d <- exp(1.96*sqrt(log(1+ccc_se^2/f0)))
  }
  
  
  alpha_estimate(estimate = ccc_hat,
                 error = ccc_se,
                 estimand = "richness",
                 name = "PoissonModel",
                 interval = c(cc + cc_excluded + f0/d, cc + cc_excluded + f0*d), 
                 type = "parametric",
                 model = "Poisson",
                 frequentist = TRUE,
                 parametric = TRUE,
                 reasonable = FALSE,
                 interval_type = "Approximate: log-normal",
                 other = list("lambda_hat" = lambda_hat,
                              "cutoff" = cutoff))
}
