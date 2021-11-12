#' modelling total diversity with random effects
#'
#' This function extends betta() to permit random effects modelling.
#'
#'
#' @param chats A vector of estimates of total diversity at different sampling
#' locations. Required with the \code{groups} argument, and optionally with the
#' \code{X} argument.
#' @param ses The standard errors in \code{chats}, the diversity estimates. This
#' can either be a vector of standard errors (with the arguments \code{chats} and
#' \code{X}), or the name of the variable in the dataframe \code{data} that contains
#' the standard errors (with the arguments \code{formula} and \code{data}).
#' @param X A numeric matrix of covariates corresponding to fixed effects. If
#' not supplied, an intercept-only model will be fit. Optional with the \code{chats}
#' and \code{groups} arguments.
#' @param groups A categorical variable representing the random-effects groups
#' that each of the estimates belong to. Required with the \code{chats} argument and
#' optionally with the \code{X} argument.
#' @param formula A formula object of the form \eqn{y ~ x | group}. Required with
#' the \code{data} argument.
#' @param data A dataframe containing the response, response standard errors, covariates,
#' and grouping variable. Required with the \code{formula} argument.
#' @return \item{table}{ A coefficient table for the model parameters. The
#' columns give the parameter estimates, standard errors, and p-values,
#' respectively. This model is only as effective as your diversity estimation
#' procedure; for this reason please confirm that your estimates are
#' appropriate and that your model is not misspecified. betta_pic may be
#' useful for this purpose.  } \item{cov}{ Estimated covariance matrix of the
#' parameter estimates.  } \item{ssq_u}{ The estimate of the heterogeneity
#' variance.  } \item{ssq_g}{ Estimates of within-group variance. The estimate
#' will be zero for groups with only one observation.  } \item{homogeneity}{
#' The test statistic and p-value for the test of homogeneity.  }
#' \item{global}{ The test statistic and p-value for the test of model
#' explanatory power.  } \item{blups}{ The conditional expected values of the
#' diversity estimates (conditional on the random effects). Estimates of
#' variability for the random effects case are unavailable at this time; please
#' contact the maintainer if needed. \item{function.args}{A list containing
#' values initially passed to betta_random.}  }
#' @author Amy Willis
#'
#' @importFrom stats coef dexp dgeom dnbinom dpois fitted lm model.matrix nls optim
#' @importFrom stats pchisq pnorm predict quantile rbeta rbinom rnbinom rnorm runif sd var vcov
#' @import lme4
#'
#' @seealso \code{\link{betta}};
#' @references Willis, A., Bunge, J., and Whitman, T. (2015). Inference for
#' changes in biodiversity. \emph{arXiv preprint.}
#' @keywords diversity
#' @examples
#'
#' df <- data.frame(chats = c(2000, 3000, 4000, 3000), ses = c(100, 200, 150, 180), Cont_var = c(100, 150, 100, 50), groups = c("a", "a", "a", "b"))
#'
#' # formula notation
#' betta_random(formula = chats ~ Cont_var| groups, ses = ses,  data = df)
#'
#' # direct input
#' betta_random(c(2000, 3000, 4000, 3000), c(100, 200, 150, 180), X = cbind(Int = 1, Cont_var = c(100, 150, 100, 50)), groups = c("a", "a", "a", "b"))
#'
#' ## handles missing data
#' betta_random(c(2000, 3000, 4000, 3000), c(100, 200, 150, NA), groups = c("a", NA,
#'     "b", "b"))
#'
#' @export
betta_random <- function(chats = NULL, ses, X = NULL, groups = NULL, formula = NULL, data = NULL) {
  if (!is.null(formula)) {
    if (is.null(data)) {
      stop("Please include a dataframe that corresponds with your formula.")
    }
  } else {
    if (is.null(chats)) {
      stop("Please include 'ses' along with either the arguments 'chats' and 'groups'
           or the arguments 'formula' and 'data'.")
    }
    if (is.null(groups)) {
      stop("Please include a vector of group memberships as 'groups'.")
    }
  }
  if (!is.null(formula)) {
    ses <- data[,deparse(substitute(ses))]
    group_var <- lme4:::barnames(lme4::findbars(lme4:::RHSForm(formula)))
    groups <- data[,group_var]
    full_form <- lme4::subbars(formula)
    sm_form <- update(full_form, paste("~.-",group_var))
    X <- stats::model.matrix(sm_form, data)
    chats <- stats::model.response(stats::model.frame(formula = sm_form, data = data))
  }
  if (isTRUE(is.null(X))) { X <- matrix(rep(1,length(chats)),ncol=1) }
  consider <- !(is.na(chats) | is.na(ses) | is.na(groups) | apply(is.na(X),1,sum))

  chats_effective <- chats[consider]; ses_effective <- ses[consider]; X_effective <- as.matrix(X[consider,])
  groups_effective <- groups[consider]
  n <- dim(X_effective)[1]; p <- dim(X_effective)[2]; gs <- length(unique(groups_effective))


  # check for design matrix that isn't full rank and throw error
  rank <- qr(X)$rank
  if (rank < ncol(X)) {
    stop(    "Your design matrix is not full rank. We recommend that you
         examine your design matrix for linear dependency and remove
         redundant columns.")
  }

  likelihood <- function(input) {
    ssq_u <- input[1]
    beta <- input[2:(p+1)]
    ssq_group <- input[(p+2):(p+gs+1)]
    group_variance <- ssq_group[groups_effective]
    W <- diag(1/(ssq_u+ses_effective^2+group_variance))
    -0.5*(sum(log(ssq_u+ses_effective^2+group_variance)+(chats_effective - X_effective %*% beta)^2/(ssq_u+ses_effective^2+group_variance)) + log(det(t(X_effective) %*% W %*% X_effective)))
  }

  within_groups_start <- c(by(chats_effective, groups_effective, var, simplify = T))
  within_groups_start[which(is.na(within_groups_start))]  <- 0

  initial_est <- c(var(chats_effective),
                   solve(t(X_effective) %*% X_effective) %*% t(X_effective) %*% chats_effective,
                   within_groups_start)

  output <- optim(initial_est,
                  likelihood,
                  hessian=FALSE,
                  control=list(fnscale=-1),
                  lower=c(0, rep(-Inf, p), rep(0, gs)),
                  method="L-BFGS-B")
  ssq_u <- output$par[1]
  beta <- output$par[2:(p+1)]
  ssq_group <- output$par[(p+2):(p+gs+1)]

  W <- diag(1/(ssq_u+ses_effective^2+ssq_group[groups_effective]))
  vars <- 1/diag(t(X_effective) %*% W %*% X_effective)

  global <- t(beta) %*% (t(X_effective) %*% W %*% X_effective) %*% beta ## global test

  Q <- sum((chats_effective - X_effective %*% beta)^2/ses_effective^2)

  mytable <- list()
  mytable$table <- cbind("Estimates"=beta,"Standard Errors"=sqrt(vars),"p-values"=round(2*(1-pnorm(abs(beta/sqrt(vars)))),3))
  try(rownames(mytable$table) <- colnames(X), silent = T)
  mytable$cov <- solve(t(X_effective) %*% W %*% X_effective)
  mytable$ssq_u <- ssq_u
  mytable$ssq_group <- ssq_group
  mytable$homogeneity <- c(Q,1-pchisq(Q,n-p))
  mytable$global <- c(global,1-pchisq(global,p-1))

  us <-  c(ssq_u*W %*% (chats_effective - X_effective %*% beta))
  blups <- rep(NA, length(chats))
  blups[consider] <- c(X_effective %*% beta + us)
  mytable$blups <- blups

  function.args <- list("chats" = chats,
       "ses" = ses,
       "X" = X,
       "groups" = groups,
       "formula" = formula,
       "data" = data,
       "model_type" = "mixed")

  mytable$function.args <- function.args

  return(mytable)
}
