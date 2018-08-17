#' PoissonModel
#' 
#' A model to estimate the number of missing taxa under a Poisson Model
#' 
#' @param input_data A frequency count table
#' @param cutoff The largest frequency to use for predicting f0. Defaults to 10.
#' 
#' @import magrittr 
#' @importFrom stats uniroot
#' 
#' @export
poisson_model <- function(input_data, 
                          cutoff = 10) {
  
  my_data <- convert(input_data)
  
  input_data <- convert(input_data)
  included <- input_data[input_data$index <= cutoff, ]
  excluded <- input_data[input_data$index > cutoff, ]
  
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
  
  # eqn nonlinear; so use Newton's method
  lambda_hat <- uniroot(poisson_fn, c(0.0001, 1000000))$root
  
  ccc_subset <- cc / (1-exp(-lambda_hat)) 
  ccc_hat <- ccc_subset + cc_excluded
  
  ccc_se <- sqrt(ccc_subset/(exp(lambda_hat)-1-lambda_hat))
  
  f0 <- ccc_hat - (cc + cc_excluded)
  d <- ifelse(f0 == 0,
              1,
              exp(1.96*sqrt(log(1+ccc_se^2/f0))))
  
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
