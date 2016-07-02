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
