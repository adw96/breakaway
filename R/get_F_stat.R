#' Calculate F statistic under null hypothesis LB = 0 using output from betta()
#' or betta_random()
#'
#' This function calculates an F statistic for a test of null hypothesis LB = 0
#' (against an unrestricted alternative) where B is a
#' vector of p fixed effects returned by betta() or betta_random() and L is an
#' m x p matrix with linearly independent rows.
#'
#' @param fitted_betta A fitted betta object -- i.e., the output of either betta() or
#' betta_random() -- containing fixed effect estimates of interest.
#' @param L An m x p matrix defining the null LB = 0. L must have full row rank.
#' @return A list containing
#' \item{F_stat}{The calculated F statistic}
get_F_stat <- function(fitted_betta,
                            L){
  #store estimated covariance matrix for beta in C
  C <- fitted_betta$cov
  #pull estimate for beta out of betta output, format as column matrix
  beta_hat <- matrix(fitted_betta$table[,"Estimates"],ncol = 1)
  #compute LB
  LB <- L%*%beta_hat
  #determine numerator degrees of freedom
  q <- nrow(L)
  #determine denominator degrees of freedom
  v <- length(fitted_betta$blups) # is this always equal to n?
  v <- v - q

  Fstat <-  t(LB)%*%solve(L%*%C%*%t(L))%*%LB/q
  Fstat <- as.numeric(Fstat)

  return(Fstat)
}
