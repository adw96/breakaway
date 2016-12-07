sample_size_est <- function(control_group_est, se_est, diff = 5, alpha = 0.05, prop = 0.8, samples = 20, precision = 5) {
  n <- 5
  pvalues <- rep(1, samples)
  
  while (mean(pvalues < alpha) < prop) {
    iter <- 1
    while (iter <= samples) {
      create_ests <- c(rnorm(n, control_group_est, se_est), rnorm(n, control_group_est + diff, se_est))
      ses <- rep(se_est, 2*n)
      design_matrix <- cbind(1, c(rep(0, n), rep(1, n)))
      pp <- try(breakaway::betta(create_ests, ses, design_matrix)$table[2,3], silent = T)
      if (class(pp) == "numeric") {
        pvalues[iter] <- pp
        iter <- iter + 1
      }
    }
    cat("Power at a sample size of ", n, ": ", 100*mean(pvalues < alpha), "%\n", sep="")
    n <- n + precision
  }
  cat("Sample size needed: ", n - precision, "\n") ## bc incremented
}