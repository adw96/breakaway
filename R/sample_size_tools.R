#' Estimate the sample size needed to do an unpaired one-way test using betta
#' 
#' Estimate the sample size needed to do an unpaired one-way test using betta
#' 
#' 
#' @param control_group_est An estimate of the alpha diversity parameter for
#' the control group
#' @param se_est An estimate of the (common) standard deviation
#' @param diff An estimate of the difference between the control and
#' non-control groups
#' @param alpha Minimum significance level desired
#' @param prop Desired power
#' @param samples Number of bootstrap resamples used to estimate the sample
#' size. Increase for a more accurate estimate.
#' @param precision How much to increment the sample size as we try to increase
#' the power
#' @return An estimate of the necessary sample size and some details
#' @export sample_size_estimate
sample_size_estimate <- function(control_group_est, se_est, diff = 5, 
                                 alpha = 0.05, prop = 0.8, samples = 20, 
                                 precision = 5) {
  n <- 5
  pvalues <- rep(1, samples)
  
  while (mean(pvalues < alpha) < prop) {
    iter <- 1
    while (iter <= samples) {
      create_ests <- c(rnorm(n, control_group_est, se_est), rnorm(n, control_group_est + diff, se_est))
      ses <- rep(se_est, 2*n)
      design_matrix <- cbind(1, c(rep(0, n), rep(1, n)))
      pp <- try(betta(create_ests, ses, design_matrix)$table[2,3], silent = T)
      if (class(pp) == "numeric") {
        pvalues[iter] <- pp
        iter <- iter + 1
      }
    }
    cat("Power at a sample size of ", n, ": ", 100*mean(pvalues < alpha), "%\n", sep="")
    n <- n + precision
  }
  cat("Sample size needed: ", n - precision, "\n") ## bc incremented
  cat("Note: This is number of subjects *per group* i.e.", n-precision, "from the control", 
      "group and", n-precision, "from the test group. It should be considered a lower bound",
      "at best! Please consult Amy with any questions. ")
}

















#' Plot the power obtained with sample size
#' 
#' Plot the power obtained with sample size
#' 
#' 
#' @param control_group_est An estimate of the alpha diversity parameter for
#' the control group
#' @param se_est An estimate of the (common) standard deviation
#' @param diff An estimate of the difference between the control and
#' non-control groups
#' @param samples Number of bootstrap resamples used to estimate the sample
#' size. Increase for a more accurate estimate.
#' @return A plot of the power with the sample size
#' @export sample_size_figure
sample_size_figure <- function(control_group_est, se_est, diff = 5, samples = 20) {
  nn <- seq(5, 50, length.out = 10)
  pvalues <- rep(NA, samples)
  
  plot(0, 0, type="n", xlim = c(0, max(nn)), ylim = c(0, 0.3), xlab = "Sample Size", ylab = "p-value")
  for (n in nn) {
    iter <- 1
    while (iter <= samples) {
      create_ests <- c(rnorm(n, control_group_est, se_est), rnorm(n, control_group_est + diff, se_est))
      ses <- rep(se_est, 2*n)
      design_matrix <- cbind(1, c(rep(0, n), rep(1, n)))
      pp <- try(betta(create_ests, ses, design_matrix)$table[2,3])
      if (class(pp) == "numeric") {
        pvalues[iter] <- pp
        iter <- iter + 1
      } 
    }
    points(rep(n, samples), pvalues)
  }
}
