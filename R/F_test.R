
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
#' @param method A character variable indicating which method should be used to
#' estimate the distribution of the test statistic under the null.
#' @param nboot Number of bootstrap samples to use if method = "bootstrap".
#' Ignored if method = "asymptotic".
#
#' @return A list containing
#' \item{pval}{The p-value}
#' \item{F_stat}{The calculated F statistic}
#' \item{boot_F}{A vector of bootstrapped F statistics if bootstrap has been used.
#' Otherwise NULL.}
#' @author David Clausen
#'
#' @seealso \code{\link{betta}};
#' @seealso \code{\link{betta_random}};
#' @seealso \code{\link{betta}};
#' @seealso \code{\link{get_F_stat}}
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
                   L,
                   method = "bootstrap",
                   nboot = 1000){

  if(!(method %in% c("bootstrap","asymptotic"))){
    stop("Value provided for argument `method` must be equal to `asymptotic` or `bootstrap.`")
  }

Fstat <- get_F_stat(fitted_betta,
                    L)
q <- nrow(L)
v <- length(fitted_betta$blups) - nrow(fitted_betta$table)

if(method == "asymptotic"){
pval <- pf(Fstat,q,v,lower.tail = F)
}

if(method == "bootstrap"){
if(fitted_betta$function.args$model_type == "fixed"){
  if(!all(sapply(as.numeric(unique(L)), function(x) x %in% c(0,1)))){
    stop("Parametric bootstrap only implemented for tests of submodels
         derived from full model by setting elements of beta to zero.
         Use method = 'asymptotic.'")
  }
  null_X = fitted_betta$function.args$X[,which(apply(L,2,sum) == 0),drop = F]
  null_fit <- betta(chats = fitted_betta$function.args$chats,
                    ses = fitted_betta$function.args$ses,
                    X = null_X,
                    data = fitted_betta$function.args$data
                    )
  sims <- simulate_betta(null_fit,
                         nboot)
  boot_F <- sapply(1:nboot,
                   function(i){
                    boot_fit <-  betta(chats = sims[[i]],
                                       ses = fitted_betta$function.args$ses,
                                       X = fitted_betta$function.args$X);
                    return(get_F_stat(boot_fit,L))})

}

  if(fitted_betta$function.args$model_type == "mixed"){
    if(!all(sapply(as.numeric(unique(L)), function(x) x %in% c(0,1)))){
      stop("Parametric bootstrap only implemented for tests of submodels
         derived from full model by setting elements of beta to zero.
         Use method = 'asymptotic.'")
    }
    null_X = fitted_betta$function.args$X[,which(apply(L,2,sum) == 0),drop = F]
    null_fit <- betta_random(chats = fitted_betta$function.args$chats,
                      ses = fitted_betta$function.args$ses,
                      X = null_X,
                      groups = fitted_betta$function.args$groups,
                      data = fitted_betta$function.args$data
    )
    sims <- simulate_betta_random(null_fit,
                           nboot)
    boot_F <- sapply(1:nboot,
                     function(i){
                       boot_fit <-  betta_random(chats = sims[[i]],
                                          ses = fitted_betta$function.args$ses,
                                          X = fitted_betta$function.args$X,
                                          groups = fitted_betta$function.args$groups);
                       return(get_F_stat(boot_fit,L))})

  }
  pval <- mean(boot_F>=Fstat)
} else{
  boot_F <- NULL
}

return(list("pval" = pval,
            "Fstat" = Fstat,
            "boot_F" = boot_F))

}
