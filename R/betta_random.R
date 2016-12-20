#' modelling total diversity with random effects
#' 
#' This function extends betta() to permit random effects modelling.
#' 
#' 
#' @param chats A vector of estimates of total diversity at different sampling
#' locations.
#' @param ses The standard errors in \samp{chats}, the diversity estimates.
#' @param X A numeric matrix of covariates corresponding to fixed effects. If
#' not supplied, an intercept-only model will be fit.
#' @param groups A categorical variable representing the random-effects groups
#' that each of the estimates belong to.
#' @return \item{table}{ A coefficient table for the model parameters. The
#' columns give the parameter estimates, standard errors, and p-values,
#' respectively. This model is only as effective as your diversity estimation
#' procedure; for this reason please confirm that your estimates are
#' appropriate and that your model is not misspecified. \samp{betta_pic} may be
#' useful for this purpose.  } \item{cov}{ Estimated covariance matrix of the
#' parameter estimates.  } \item{ssq_u}{ The estimate of the heterogeneity
#' variance.  } \item{ssq_g}{ Estimates of within-group variance. The estimate
#' will be zero for groups with only one observation.  } \item{homogeneity}{
#' The test statistic and p-value for the test of homogeneity.  }
#' \item{global}{ The test statistic and p-value for the test of model
#' explanatory power.  } \item{blups}{ The conditional expected values of the
#' diversity estimates (conditional on the random effects). Estimates of
#' variability for the random effects case are unavailable at this time; please
#' contact the maintainer if needed.  }
#' @author Amy Willis
#' @seealso \code{\link{betta}}; \code{\link{betta_pic}}
#' @references Willis, A., Bunge, J., and Whitman, T. (2015). Inference for
#' changes in biodiversity. \emph{arXiv preprint.}
#' @keywords diversity
#' @examples
#' 
#' betta_random(c(2000, 3000, 4000, 3000), c(100, 200, 150, 180),
#'              X = cbind("Int"=1, "Cont_var"=c(100, 150, 100, 50)),
#'              groups = c("a", "a", "a", "b"))
#' 
#' ## handles missing data
#' betta_random(c(2000, 3000, 4000, 3000), c(100, 200, 150, NA),
#'              groups= c("a", NA, "b", "b"))
#' 
#' @export
betta_random <- function(chats, ses, X=NA, groups) {
  if (isTRUE(is.na(X))) { X <- matrix(rep(1,length(chats)),ncol=1) }
  consider <- !(is.na(chats) | is.na(ses) | is.na(groups) | apply(is.na(X),1,sum))

  chats_effective <- chats[consider]; ses_effective <- ses[consider]; X_effective <- as.matrix(X[consider,])
  groups_effective <- groups[consider]
  n <- dim(X_effective)[1]; p <- dim(X_effective)[2]; gs <- length(unique(groups_effective))

  likelihood <- function(input) {
    ssq_u <- input[1]
    beta <- input[2:(p+1)]
    ssq_group <- input[(p+2):(p+gs+1)]
    group_variance <- ssq_group[groups_effective]
    W <- diag(1/(ssq_u+ses_effective^2+group_variance))
    -0.5*(sum(log(ssq_u+ses_effective^2+group_variance)+(chats_effective-X_effective%*%beta)^2/(ssq_u+ses_effective^2+group_variance)) + log(det(t(X_effective)%*%W%*%X_effective)))
  }

  within_groups_start <- c(by(chats_effective, groups_effective, var, simplify = T))
  within_groups_start[which(is.na(within_groups_start))]  <- 0

  mystart <- c(var(chats_effective),
               solve(t(X_effective)%*%X_effective)%*%t(X_effective)%*%chats_effective,
               within_groups_start)
  output <- optim(mystart, likelihood, hessian=FALSE, control=list(fnscale=-1),
                  lower=c(0, rep(-Inf, p), rep(0, gs)), method="L-BFGS-B")
  ssq_u <- output$par[1]
  beta <- output$par[2:(p+1)]
  ssq_group <- output$par[(p+2):(p+gs+1)]

  W <- diag(1/(ssq_u+ses_effective^2+ssq_group[groups_effective]))
  vars <- 1/diag(t(X_effective)%*%W%*%X_effective)

  global <- t(beta)%*%(t(X_effective)%*%W%*%X_effective)%*%beta ## global test

  Q <- sum((chats_effective-X_effective%*%beta)^2/ses_effective^2)

  mytable <- list()
  mytable$table <- cbind("Estimates"=beta,"Standard Errors"=sqrt(vars),"p-values"=round(2*(1-pnorm(abs(beta/sqrt(vars)))),3))
  try(rownames(mytable$table) <- colnames(X), silent = T)
  mytable$cov <- solve(t(X_effective)%*%W%*%X_effective)
  mytable$ssq_u <- ssq_u
  mytable$ssq_group <- ssq_group
  mytable$homogeneity <- c(Q,1-pchisq(Q,n-p))
  mytable$global <- c(global,1-pchisq(global,p-1))

  us <-  c(ssq_u*W%*%(chats_effective-X_effective%*%beta))
  blups <- rep(NA, length(chats))
  blups[consider] <- c(X_effective%*%beta + us)
  mytable$blups <- blups

  return(mytable)
}
