#' @export
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
      pp <- try(breakaway::betta(create_ests, ses, design_matrix)$table[2,3])
      if (class(pp) == "numeric") {
        pvalues[iter] <- pp
        iter <- iter + 1
      } 
    }
    points(rep(n, samples), pvalues)
  }
}