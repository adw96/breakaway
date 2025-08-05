#' Confidence intervals and testing for linear combinations of fixed effects estimated via betta() or betta_random()
#'
#' This function provides point estimates, standard errors, and equal-tailed confidence intervals for
#' linear combinations of fixed effects estimated via betta() or betta_random(). A p-value for a Wald
#' test of the null that the linear combination of effects is equal to zero (against a two-sided alternative)
#' is also returned.
#'
#'
#' @param fitted_betta A fitted betta object -- i.e., the output of either betta() or
#' betta_random() -- containing fixed effect estimates of interest.
#' @param linear_com The linear combination of fixed effects for which a point estimate,
#' confidence interval, and hypothesis test are to be produced.
#' @param signif_cutoff The type-I significance threshold for confidence intervals. Defaults
#' to 0.05.
#
#' @return \item{table}{ A table containing a point estimate, standard error, lower and upper confidence bounds,
#' and a p-value for the linear combination of fixed effects specified in input. The p-value is generated via a two-sided
#' Wald test of the null that the linear combination of fixed effects is equal to zero.}
#' @author David Clausen
#'
#' @seealso \code{\link{betta}};
#' @references Willis, A., Bunge, J., and Whitman, T. (2015). Inference for
#' changes in biodiversity. \emph{arXiv preprint.}
#' @keywords diversity
#' @examples
#'
#' # generate example data
#' df <- data.frame(chats = c(2000, 3000, 4000, 3000), ses = c(100, 200, 150, 180),
#'                  Cont_var = c(100, 150, 100, 50))
#'
#' # fit betta()
#' example_fit <- betta(formula = chats ~ Cont_var, ses = ses, data = df)
#'
#' # generate point estimate and 95% CI for mean richness at Cont_var = 125
#'
#' betta_lincom(fitted_betta = example_fit,
#' linear_com = c(1, 125)) # this tells betta_lincom to estimate value of beta_0 + 125*beta_1,
#' # where beta_0 is the intercept, and beta_1 is the (true value of the) coefficient on Cont_var
#'
#'
#' @export
betta_lincom <- function(fitted_betta,
                         linear_com,
                         signif_cutoff = 0.05){

  cov_mat <- fitted_betta$cov

  est <- sum(linear_com*fitted_betta$table[,1])

  se <- sqrt(matrix(linear_com, nrow = 1) %*%
               cov_mat %*%
               matrix(linear_com, ncol = 1))

  upper <- est + qnorm(1 - signif_cutoff/2)*se
  lower <- est - qnorm(1 - signif_cutoff/2)*se
  pval <- 2*(1-pnorm(abs(est/se)))


  results <- matrix( c(est, se, lower, upper, pval), nrow = 1)

  colnames(results) <- c("Estimates", "Standard Errors", "Lower CIs", "Upper CIs", "p-values")

  results <- as.data.frame(results)

  if(sum(pval ==0)>0){
    results$`p-values` <- as.character(results$`p-values`)
    results$`p-values`[pval == 0] <- " < 1e-20"
  } else{
    #so class of p-values column doesn't depend on value of p-values
    results$`p-values` <- as.character(signif(results$`p-values`,3))
  }
  return(results)
}
