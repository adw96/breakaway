betta <- function(chats,ses,X=NA) {
  if (isTRUE(is.na(X))) { X <- matrix(rep(1,length(chats)),ncol=1) }
  consider <- !(is.na(chats) | is.na(ses) | apply(is.na(X),1,sum))

  chats_effective <- chats[consider]; ses_effective <- ses[consider]; X_effective <- as.matrix(X[consider,])
  n <- dim(X_effective)[1]; p <- dim(X_effective)[2]

  likelihood <- function(input) {
    ssq_u <- input[1]
    beta <- input[2:length(input)]
    W <- diag(1/(ssq_u+ses_effective^2))
    -0.5*(sum(log(ssq_u+ses_effective^2)+(chats_effective-X_effective%*%beta)^2/(ssq_u+ses_effective^2)) + log(det(t(X_effective)%*%W%*%X_effective)))
  }

  mystart <- c(var(chats_effective),solve(t(X_effective)%*%X_effective)%*%t(X_effective)%*%chats_effective)
  output <- optim(mystart,likelihood,hessian=FALSE,control=list(fnscale=-1),lower=c(0,-Inf), method="L-BFGS-B") # fnscale => maximises if -ve
  ssq_u <- output$par[1]
  beta <- output$par[2:length(output$par)]

  W <- diag(1/(ssq_u+ses_effective^2))
  vars <- 1/diag(t(X_effective)%*%W%*%X_effective)

  global <- t(beta)%*%(t(X_effective)%*%W%*%X_effective)%*%beta ## global test

  Q <- sum((chats_effective-X_effective%*%beta)^2/ses_effective^2)
  R <- diag(ses_effective^2); G <- diag(ssq_u,n)

  getvar <- function() {
    C <- matrix(NA,nrow=n+p,ncol=n+p);
    C[1:p,1:p] <- t(X_effective)%*%solve(R)%*%X_effective;
    C[(p+1):(n+p),(p+1):(n+p)] <- solve(R)+solve(G);
    C[1:p,(p+1):(n+p)] <- t(X_effective)%*%solve(R);
    C[(p+1):(n+p),1:p] <- solve(R)%*%X_effective;
    return(solve(C))
  }

  mytable <- list()
  mytable$table <- cbind("Estimates"=beta,"Standard Errors"=sqrt(vars),"p-values"=round(2*(1-pnorm(abs(beta/sqrt(vars)))),3))
  mytable$cov <- solve(t(X_effective)%*%W%*%X_effective)
  mytable$ssq_u <- ssq_u
  mytable$homogeneity <- c(Q,1-pchisq(Q,n-p))
  mytable$global <- c(global,1-pchisq(global,p-1))

  us <-  c(ssq_u*W%*%(chats_effective-X_effective%*%beta))
  blups <- rep(NA, length(chats))
  blups[consider] <- c(X_effective%*%beta + us)
  mytable$blups <- blups

  if (class(try(getvar(), silent=T)) != "try-error") {
    C <- getvar()
    blupvars <- rep(NA, length(chats))
    blupvars[consider] <- c(sqrt(diag(cbind(X_effective,diag(1,n))%*%C%*%t(cbind(X_effective,diag(1,n))))))
    mytable$blupses <- blupvars
  }

  return(mytable)
}
