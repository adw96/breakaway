
#' Conduct F test of null hypothesis LB = 0 using output from betta()
#' or betta_random()
#'
#' This function performs an F-test of a null hypothesis LB = 0 where B is a
#' vector of p fixed effects returned by betta() or betta_random() and L is an
#' m x p matrix with linearly independent rows.
#'
#'
#' @param fitted_betta A fitted betta object -- i.e., the output of either betta() or
#' betta_random() -- containing fixed effect estimates of interest.
#' @param L An m x p matrix defining the null LB = 0. L must have full row rank.

#
#' @return A list containing
#' \item{pval}{The p-value}
#' @author David Clausen
#'
#' @seealso \code{\link{betta}};
#' @seealso \code{\link{betta_random}};
#' @seealso \code{\link{betta}};
#' @references Willis, A., Bunge, J., and Whitman, T. (2015). Inference for
#' changes in biodiversity. \emph{arXiv preprint.}
#' @keywords diversity
#' @examples
#'
#' # generate example data
#' df <- data.frame(chats = c(2000, 3000, 4000, 3000,
#' 2000, 3000, 4000, 3000), ses = c(100, 200, 150, 180,
#' 100, 200, 150, 180),
#'                  Cont_var = c(100, 150, 100, 50,
#'                  100, 150, 100, 50),
#'                  Cont_var_2 = c(50,200,25,125,
#'                  50,200,25,125))
#'
#' # fit betta()
#' example_fit <- betta(formula = chats ~ Cont_var + Cont_var_2, ses = ses, data = df)
#'
#'
#' # construct L for hypothesis that B_cont_var = B_cont_var_2 = 0
#' L <- rbind(c(0,1,0),
#'            c(0,0,1))
#'
#' F_test_results <- F_test(example_fit, L)
#'
#' @export
F_test <- function(fitted_betta,
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

pval <- pf(Fstat,q,v,lower.tail = F)

return(list("pval" = pval,
            "Fstat" = Fstat))

}
