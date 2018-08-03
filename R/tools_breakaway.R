# #' @import stats
# minibreak <- function(lhs, xs, ys, my_data, myweights) {
#   structure_1_0 <- function(x,beta0,beta1) (beta0+beta1*x)/(1+x)
#   structure_1_1 <- function(x,beta0,beta1,alpha1) (beta0+beta1*x)/(1+alpha1*x)
#   structure_2_1 <- function(x,beta0,beta1,alpha1,beta2) (beta0+beta1*x+beta2*x^2)/(1+alpha1*x)
#   structure_2_2 <- function(x,beta0,beta1,alpha1,beta2,alpha2) (beta0+beta1*x+beta2*x^2)/(1+alpha1*x+alpha2*x^2)
#   structure_3_2 <- function(x,beta0,beta1,alpha1,beta2,alpha2,beta3) (beta0+beta1*x+beta2*x^2+beta3*x^3)/(1+alpha1*x+alpha2*x^2)
#   structure_3_3 <- function(x,beta0,beta1,alpha1,beta2,alpha2,beta3,alpha3) (beta0+beta1*x+beta2*x^2+beta3*x^3)/(1+alpha1*x+alpha2*x^2+alpha3*x^3)
#   structure_4_3 <- function(x,beta0,beta1,alpha1,beta2,alpha2,beta3,alpha3,beta4) (beta0+beta1*x+beta2*x^2+beta3*x^3+beta4*x^4)/(1+alpha1*x+alpha2*x^2+alpha3*x^3)
#   structure_4_4 <- function(x,beta0,beta1,alpha1,beta2,alpha2,beta3,alpha3,beta4,alpha4) (beta0+beta1*x+beta2*x^2+beta3*x^3+beta4*x^4)/(1+alpha1*x+alpha2*x^2+alpha3*x^3+alpha4*x^4)
#   
#   all_ests <- list()
#   
#   model_1_0 <-  lm(ys~xs,weights=myweights)
#   start_1_0 <- list(beta0 = model_1_0$coef[1], beta1 = model_1_0$coef[2])
#   
#   lhs2 <- lhs; lhs2$x <- xs # otherwise singularity
#   result_1_0 <- try( nls( lhs2$y ~ structure_1_0(x,beta0,beta1), data = lhs2, start_1_0, weights=myweights), silent=T )
#   if(class(result_1_0) == "try-error") {
#     coef_1_0 <- c(start_1_0$beta0,start_1_0$beta1)
#   } else {
#     coef_1_0 <- coef(result_1_0)
#   }
#   
#   beta0_tilde <- coef_1_0[1]
#   beta1_tilde <- coef_1_0[2]
#   xbar <- mean(xs)
#   start_1_1 <- list(beta0 = (beta0_tilde+beta1_tilde*xbar)/(1+xbar), beta1 = beta1_tilde/(1+xbar), alpha1 = 1/(1+xbar))
#   result_1_1 <- try( nls( lhs$y ~ structure_1_1(x,beta0,beta1,alpha1), data = lhs, start_1_1, weights=myweights), silent=T )
#   
#   if(class(result_1_1) == "try-error") {
#     temp_lm <- lm(lhs$y~lhs$x+I((lhs$x)^2),weights=myweights)$coef
#     start_2_1 <- list(beta0 = temp_lm[1], beta1 = temp_lm[2], alpha1 = 0, beta2 = temp_lm[3])
#   } else {
#     start_2_1 <- as.list(c(summary(result_1_1)$coef[,1],"beta2"=0))
#   }
#   result_2_1 <- try(nls( lhs$y ~ structure_2_1(x,beta0,beta1,alpha1,beta2), data = lhs, start_2_1, weights=myweights),silent=T)
#   
#   if(class(result_2_1) == "try-error") {
#     temp_lm <- lm(lhs$y~lhs$x+I((lhs$x)^2),weights=myweights)$coef
#     start_2_2 <- list(beta0 = temp_lm[1], temp_lm[2], alpha1 = 0, beta2 = temp_lm[3], alpha2 = 0)
#   } else start_2_2 <- as.list(c(summary(result_2_1)$coef[,1],"alpha2"=0))
#   result_2_2 <- try(nls( lhs$y ~ structure_2_2(x,beta0,beta1,alpha1,beta2,alpha2), data = lhs, start_2_2, weights=myweights), silent=T)
#   
#   if(class(result_2_2) == "try-error") {
#     temp_lm <- lm(lhs$y~lhs$x+I((lhs$x)^2)+I((lhs$x)^3),weights=myweights)$coef
#     start_3_2 <- list(beta0 = temp_lm[1], beta1 = temp_lm[2], alpha1 = 0, beta2 = temp_lm[3], alpha2 = 0, beta3 = temp_lm[4])
#   } else start_3_2 <- as.list(c(summary(result_2_2)$coef[,1],"beta3"=0))
#   result_3_2 <- try( nls( lhs$y ~ structure_3_2(x,beta0,beta1,alpha1,beta2,alpha2,beta3), data = lhs, start_3_2, weights=myweights), silent=T)
#   
#   if(class(result_3_2) == "try-error") {
#     temp_lm <- lm(lhs$y~lhs$x+I((lhs$x)^2)+I((lhs$x)^3), weights=myweights)$coef
#     start_3_3 <- list(beta0 = temp_lm[1], beta1 = temp_lm[2], alpha1 = 0, beta2 = temp_lm[3], alpha2 = 0, beta3 = temp_lm[4], alpha3 = 0)
#   } else start_3_3 <- as.list(c(summary(result_3_2)$coef[,1],"alpha3"=0))
#   result_3_3 <- try( nls( lhs$y ~ structure_3_3(x,beta0,beta1,alpha1,beta2,alpha2,beta3,alpha3), data = lhs, start_3_3, weights=myweights), silent=T)
#   
#   if(class(result_3_3) == "try-error") {
#     temp_lm <- lm(lhs$y~lhs$x+I((lhs$x)^2)+I((lhs$x)^3)+I((lhs$x)^4),weights=myweights)$coef
#     start_4_3 <- list(beta0 = temp_lm[1], beta1 = temp_lm[2], alpha1 = 0, beta2 = temp_lm[3], alpha2 = 0,beta3 = temp_lm[4], alpha3 = 0,beta4 = temp_lm[5])
#   } else start_4_3 <- as.list(c(summary(result_3_3)$coef[,1],"beta4"=0))
#   result_4_3 <- try( nls( lhs$y ~ structure_4_3(x,beta0,beta1,alpha1,beta2,alpha2,beta3,alpha3,beta4), data = lhs, start_4_3, weights=myweights), silent=T)
#   
#   if(class(result_4_3) == "try-error") {
#     temp_lm <- lm(lhs$y~lhs$x+I((lhs$x)^2)+I((lhs$x)^3)+I((lhs$x)^4),weights=myweights)$coef
#     start_4_4 <- list(beta0 = temp_lm[1], beta1 = temp_lm[2], alpha1 = 0, beta2 = temp_lm[3], alpha2 = 0,
#                       beta3 = temp_lm[4], alpha3 = 0, beta4 = temp_lm[5], alpha4 = 0)
#   } else start_4_4 <- as.list(c(summary(result_4_3)$coef[,1],"alpha4"=0))
#   result_4_4 <- try( nls( lhs$y ~ structure_4_4(x,beta0,beta1,alpha1,beta2,alpha2,beta3,alpha3,beta4,alpha4), data = lhs, start_4_4, weights=myweights), silent=T)
#   
#   info <- list("model_1_0"=result_1_0,"model_1_1"=result_1_1,
#                "model_2_1"=result_2_1,"model_2_2"=result_2_2,"model_3_2"=result_3_2,
#                "model_3_3"=result_3_3,"model_4_3"=result_4_3,"model_4_4"=result_4_4)
#   
#   converged <- lapply(info,class)=="nls"
#   
#   workinginfo <- info[converged]
#   
#   b0est <- lapply(workinginfo,predict,list(x=-xbar))
#   if (length(b0est$model_1_0)!=0) b0est$model_1_0 <- predict(workinginfo$model_1_0,list(x=0))
#   f0est <- lapply(b0est,function(x) my_data[1,2]/x)
#   
#   rootsokay <- as.logical(lapply(workinginfo,rootcheck,lhs,nof1=FALSE))
#   
#   sqerrors <- lapply(workinginfo,sqerror,lhs)
#   residses <- lapply(workinginfo,residse)
#   
#   useable <- rootsokay & (f0est>0)
#   
#   workinginfo$useful <- cbind(f0est,rootsokay,sqerrors,residses,useable)
#   return(workinginfo)
# }

rootcheck <- function(model,lhs,nof1=FALSE) {
  if(length(coef(model))==2)  {
    root <- -coef(model)[1]/coef(model)[2]
  } else {
    root <- Re(polyroot(c(1,coef(model)[substring(names(coef(model)),1,1)=="a"]))[
      abs(Im(polyroot(c(1,coef(model)[substring(names(coef(model)),1,1)=="a"]))))<10^-5])
  }
  if (length(root)==0) {
    outcome <- 1
  } else if (sum((root > ifelse(nof1,min(lhs$x-1),min(lhs$x))) & (root < max(lhs$x)))==0 ) { #uses sum of pointwise booleans to ensure any roots are caught
    outcome <- 1
  } else {
    outcome <- 0
  }
  return(outcome)
}

sqerror <- function(model, lhs) {
  fits <- fitted(model)
  return(sum(((lhs$y-fits)^2)/fits))
}

residse <- function(model) {
  return(summary(model)$sigma)
}


#' #' @import stats
#' minibreak_nof1 <- function(lhs, xs, ys, my_data, myweights) {
#'   structure_1_0 <- function(x,beta0,beta1) (beta0+beta1*x)/(1+x)
#'   structure_1_1 <- function(x,beta0,beta1,alpha1) (beta0+beta1*x)/(1+alpha1*x)
#'   structure_2_1 <- function(x,beta0,beta1,alpha1,beta2) (beta0+beta1*x+beta2*x^2)/(1+alpha1*x)
#'   structure_2_2 <- function(x,beta0,beta1,alpha1,beta2,alpha2) (beta0+beta1*x+beta2*x^2)/(1+alpha1*x+alpha2*x^2)
#'   structure_3_2 <- function(x,beta0,beta1,alpha1,beta2,alpha2,beta3) (beta0+beta1*x+beta2*x^2+beta3*x^3)/(1+alpha1*x+alpha2*x^2)
#'   structure_3_3 <- function(x,beta0,beta1,alpha1,beta2,alpha2,beta3,alpha3) (beta0+beta1*x+beta2*x^2+beta3*x^3)/(1+alpha1*x+alpha2*x^2+alpha3*x^3)
#'   structure_4_3 <- function(x,beta0,beta1,alpha1,beta2,alpha2,beta3,alpha3,beta4) (beta0+beta1*x+beta2*x^2+beta3*x^3+beta4*x^4)/(1+alpha1*x+alpha2*x^2+alpha3*x^3)
#'   structure_4_4 <- function(x,beta0,beta1,alpha1,beta2,alpha2,beta3,alpha3,beta4,alpha4) (beta0+beta1*x+beta2*x^2+beta3*x^3+beta4*x^4)/(1+alpha1*x+alpha2*x^2+alpha3*x^3+alpha4*x^4)
#'   
#'   all_ests <- list()
#'   
#'   model_1_0 <-  lm(ys~xs,weights=myweights)
#'   start_1_0 <- list(beta0 = model_1_0$coef[1], beta1 = model_1_0$coef[2])
#'   
#'   lhs2 <- lhs; lhs2$x <- xs # otherwise singularity
#'   result_1_0 <- try( nls( lhs2$y ~ structure_1_0(x,beta0,beta1), data = lhs2, start_1_0, weights=myweights), silent=T )
#'   if(class(result_1_0) == "try-error") {
#'     coef_1_0 <- c(start_1_0$beta0,start_1_0$beta1)
#'   } else {
#'     coef_1_0 <- coef(result_1_0)
#'   }
#'   
#'   beta0_tilde <- coef_1_0[1]
#'   beta1_tilde <- coef_1_0[2]
#'   xbar <- mean(c(1,xs))
#'   start_1_1 <- list(beta0 = (beta0_tilde+beta1_tilde*xbar)/(1+xbar), beta1 = beta1_tilde/(1+xbar), alpha1 = 1/(1+xbar))
#'   result_1_1 <- try( nls( lhs$y ~ structure_1_1(x,beta0,beta1,alpha1), data = lhs, start_1_1, weights=myweights), silent=T )
#'   
#'   if(class(result_1_1) == "try-error") {
#'     temp_lm <- lm(lhs$y~lhs$x+I((lhs$x)^2),weights=myweights)$coef
#'     start_2_1 <- list(beta0 = temp_lm[1], beta1 = temp_lm[2], alpha1 = 0, beta2 = temp_lm[3])
#'   } else {
#'     start_2_1 <- as.list(c(summary(result_1_1)$coef[,1],"beta2"=0))
#'   }
#'   result_2_1 <- try(nls( lhs$y ~ structure_2_1(x,beta0,beta1,alpha1,beta2), data = lhs, start_2_1, weights=myweights),silent=T)
#'   
#'   if(class(result_2_1) == "try-error") {
#'     temp_lm <- lm(lhs$y~lhs$x+I((lhs$x)^2),weights=myweights)$coef
#'     start_2_2 <- list(beta0 = temp_lm[1], temp_lm[2], alpha1 = 0, beta2 = temp_lm[3], alpha2 = 0)
#'   } else start_2_2 <- as.list(c(summary(result_2_1)$coef[,1],"alpha2"=0))
#'   result_2_2 <- try(nls( lhs$y ~ structure_2_2(x,beta0,beta1,alpha1,beta2,alpha2), data = lhs, start_2_2, weights=myweights), silent=T)
#'   
#'   if(class(result_2_2) == "try-error") {
#'     temp_lm <- lm(lhs$y~lhs$x+I((lhs$x)^2)+I((lhs$x)^3),weights=myweights)$coef
#'     start_3_2 <- list(beta0 = temp_lm[1], beta1 = temp_lm[2], alpha1 = 0, beta2 = temp_lm[3], alpha2 = 0, beta3 = temp_lm[4])
#'   } else start_3_2 <- as.list(c(summary(result_2_2)$coef[,1],"beta3"=0))
#'   result_3_2 <- try( nls( lhs$y ~ structure_3_2(x,beta0,beta1,alpha1,beta2,alpha2,beta3), data = lhs, start_3_2, weights=myweights), silent=T)
#'   
#'   if(class(result_3_2) == "try-error") {
#'     temp_lm <- lm(lhs$y~lhs$x+I((lhs$x)^2)+I((lhs$x)^3), weights=myweights)$coef
#'     start_3_3 <- list(beta0 = temp_lm[1], beta1 = temp_lm[2], alpha1 = 0, beta2 = temp_lm[3], alpha2 = 0, beta3 = temp_lm[4], alpha3 = 0)
#'   } else start_3_3 <- as.list(c(summary(result_3_2)$coef[,1],"alpha3"=0))
#'   result_3_3 <- try( nls( lhs$y ~ structure_3_3(x,beta0,beta1,alpha1,beta2,alpha2,beta3,alpha3), data = lhs, start_3_3, weights=myweights), silent=T)
#'   
#'   if(class(result_3_3) == "try-error") {
#'     temp_lm <- lm(lhs$y~lhs$x+I((lhs$x)^2)+I((lhs$x)^3)+I((lhs$x)^4),weights=myweights)$coef
#'     start_4_3 <- list(beta0 = temp_lm[1], beta1 = temp_lm[2], alpha1 = 0, beta2 = temp_lm[3], alpha2 = 0,beta3 = temp_lm[4], alpha3 = 0,beta4 = temp_lm[5])
#'   } else start_4_3 <- as.list(c(summary(result_3_3)$coef[,1],"beta4"=0))
#'   result_4_3 <- try( nls( lhs$y ~ structure_4_3(x,beta0,beta1,alpha1,beta2,alpha2,beta3,alpha3,beta4), data = lhs, start_4_3, weights=myweights), silent=T)
#'   
#'   if(class(result_4_3) == "try-error") {
#'     temp_lm <- lm(lhs$y~lhs$x+I((lhs$x)^2)+I((lhs$x)^3)+I((lhs$x)^4),weights=myweights)$coef
#'     start_4_4 <- list(beta0 = temp_lm[1], beta1 = temp_lm[2], alpha1 = 0, beta2 = temp_lm[3], alpha2 = 0,
#'                       beta3 = temp_lm[4], alpha3 = 0, beta4 = temp_lm[5], alpha4 = 0)
#'   } else start_4_4 <- as.list(c(summary(result_4_3)$coef[,1],"alpha4"=0))
#'   result_4_4 <- try( nls( lhs$y ~ structure_4_4(x,beta0,beta1,alpha1,beta2,alpha2,beta3,alpha3,beta4,alpha4), data = lhs, start_4_4, weights=myweights), silent=T)
#'   
#'   info <- list("model_1_0"=result_1_0,"model_1_1"=result_1_1,
#'                "model_2_1"=result_2_1,"model_2_2"=result_2_2,"model_3_2"=result_3_2,
#'                "model_3_3"=result_3_3,"model_4_3"=result_4_3,"model_4_4"=result_4_4)
#'   
#'   converged <- lapply(info,class)=="nls"
#'   
#'   workinginfo <- info[converged]
#'   
#'   b0est <- lapply(workinginfo,predict,list(x=-xbar))
#'   if (length(b0est$model_1_0)!=0) b0est$model_1_0 <- predict(workinginfo$model_1_0,list(x=0))
#'   
#'   firstratioest <- lapply(workinginfo,predict,list(x=-xbar+1))
#'   if (length(firstratioest$model_1_0)!=0) firstratioest$model_1_0 <- predict(workinginfo$model_1_0,list(x=1))
#'   
#'   f1est <- lapply(firstratioest,function(x) my_data[1,2]/x)
#'   
#'   f0est <- as.numeric(f1est)/as.numeric(b0est)
#'   
#'   rootsokay <- as.logical(lapply(workinginfo,rootcheck,lhs,nof1=TRUE))
#'   
#'   sqerrors <- lapply(workinginfo,sqerror,lhs)
#'   residses <- lapply(workinginfo,residse)
#'   
#'   useable <- rootsokay & (f0est>0) & (f1est>0)
#'   
#'   workinginfo$useful <- cbind(f0est,rootsokay,sqerrors,residses,useable)
#'   return(workinginfo)
#' }
