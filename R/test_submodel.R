#' Conduct F test of null hypothesis LB = 0 using output from betta()
#' or betta_random()
#'
#' This function performs an F-test of a null hypothesis LB = 0 where B is a
#' vector of p fixed effects returned by betta() or betta_random() and L is an
#' m x p matrix with linearly independent rows.
#'
#' @param fitted_betta A fitted betta object -- i.e., the output of either betta() or
#' betta_random() -- containing fixed effect estimates of interest.
#' @param submodel_formula A formula defining which submodel to treat as the null.
#' It is not necessary to include random effects in this formula (they will be ignored
#' if included -- the submodel will be fit with the same random effect structure
#' as the full model regardless of input.)
#' @param method A character variable indicating which method should be used to
#' estimate the distribution of the test statistic under the null.
#' @param nboot Number of bootstrap samples to use if method = "bootstrap".
#' Ignored if method = "asymptotic".
#' @return A list containing
#' \item{pval}{The p-value}
#' \item{F_stat}{The calculated F statistic}
#' \item{boot_F}{A vector of bootstrapped F statistics if bootstrap has been used.
#' Otherwise NULL.}
#'
#' @author David Clausen
#'
#' @seealso \code{\link{betta}};
#' @seealso \code{\link{betta_random}};
#' @seealso \code{\link{betta}};
#' @seealso \code{\link{F_test}}
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
#' example_fit <- betta(formula = chats ~ Cont_var + Cont_var_2, ses = ses,
#' data = df)
#'
#'
#' # construct L for hypothesis that B_cont_var = B_cont_var_2 = 0
#' L <- rbind(c(0,1,0),
#'            c(0,0,1))
#'
#' F_test_results <- F_test(example_fit,
#' L,
#' nboot = 10) #nboot = 10 for speed here; recommend >= 1000 in practice
#'
#' @export
test_submodel <- function(fitted_betta,
                          submodel_formula,
                          method = "bootstrap",
                          nboot = 1000){
  # use submodel_formula to derive L
  submodel_X <- stats::model.matrix(submodel_formula,
                           fitted_betta$function.args$data)

  L <- matrix(nrow = 0, ncol = ncol(fitted_betta$function.args$X))

  for(k in 1:ncol(L)){
    if( !(colnames(fitted_betta$function.args$X)[k] %in% colnames(submodel_X))
    ){
      new_row <- matrix(0, nrow = 1, ncol = ncol(L))
      new_row[1,k] <- 1
      L <- rbind(L,new_row)
    }
  }

  return(F_test(fitted_betta, L, method, nboot))

}
