#' function for species richness estimation
#' 
#' This function implements the species richness estimation procedure outlined
#' in Willis & Bunge (2015). The diversity estimate, standard error, estimated
#' model coefficients, model details and plot of the fitted model are returned.
#' 
#' 
#' @param data The sample frequency count table for the population of interest.
#' The first row must correspond to the singletons. Acceptable formats include
#' a matrix, data frame, or file path (csv or txt). The standard frequency
#' count table format is used: two columns, the first of which contains the
#' frequency of interest (eg. 1 for singletons, species observed once, 2 for
#' doubletons, species observed twice, etc.) and the second of which contains
#' the number of species observed this many times. Frequencies (first column)
#' should be ordered least to greatest. At least 6 contiguous frequencies are
#' necessary. Do not concatenate large frequencies. See dataset apples for
#' sample formatting.
#' @param print Logical: whether the results should be printed to screen. If
#' \samp{FALSE}, answers should be set to \samp{TRUE} so that results will be
#' returned.
#' @param plot Logical: whether the data and model fit should be plotted.
#' @param answers Logical: whether the function should return an argument. If
#' \samp{FALSE}, print should be set to \samp{TRUE}.
#' @param force Logical: force \samp{breakaway} to run in the presence of
#' frequency count concatenation, i.e. combining all species observed 27 or 28
#' or 29 or more into frequency count 27. \samp{force=TRUE} will force
#' \samp{breakaway} to fit models in the presence of this. This will generally
#' result in the model being misspecified, and \samp{breakaway}'s diversity
#' estimates cannot be considered reliable in this case. The option
#' \samp{force} issue does not address censoring or small data on the rare end
#' of the frequency count tables, i.e. insufficiently many contiguous
#' frequencies.
#' @return \item{code}{ A category representing algorithm behaviour.
#' \samp{code=1} indicates no nonlinear models converged and the transformed
#' WLRM diversity estimate of Rocchetti et. al. (2011) is returned.
#' \samp{code=2} indicates that the iteratively reweighted model converged and
#' was returned. \samp{code=3} indicates that iterative reweighting did not
#' converge but a model based on a simplified variance structure was returned
#' (in this case, the variance of the frequency ratios is assumed to be
#' proportional to the denominator frequency index). Please peruse your fitted
#' model before using your diversity estimate.  } \item{name}{ The ``name'' of
#' the selected model. The first integer represents the numerator polynomial
#' degree and the second integer represents the denominator polynomial degree
#' of the model for the frequency ratios. See Willis & Bunge (2015) for
#' details.  } \item{para}{ Estimated model parameters and standard errors.  }
#' \item{est}{ The estimate of total (observed plus unobserved) diversity.  }
#' \item{seest}{ The standard error in the diversity estimate.  } \item{full}{
#' The chosen nonlinear model for frequency ratios.  } \item{ci}{ An asymmetric
#' 95\% confidence interval for diversity.  }
#' @note \samp{breakaway} presents an estimator of species richness that is
#' well-suited to the high-diversity/microbial setting. However, many microbial
#' datasets display more diversity than the Kemp-type models can permit. In
#' this case, the log-transformed WLRM diversity estimator of Rocchetti et. al.
#' (2011) is returned. The authors' experience suggests that some datasets that
#' require the log-transformed WLRM contain ``false'' diversity, that is,
#' diversity attributable to sequencing errors (via an inflated singleton
#' count). The authors encourage judicious use of diversity estimators when the
#' dataset may contain these errors, and recommend the use of
#' \code{\link{breakaway_nof1}} as an exploratory tool in this case.
#' @author Amy Willis
#' @seealso \code{\link{breakaway_nof1}}; \code{\link{apples}}
#' @references Willis, A. and Bunge, J. (2015). Estimating diversity via
#' frequency ratios. \emph{Biometrics}, \bold{71}(4), 1042--1049.
#' 
#' Rocchetti, I., Bunge, J. and Bohning, D. (2011). Population size estimation
#' based upon ratios of recapture probabilities. \emph{Annals of Applied
#' Statistics}, \bold{5}.
#' @keywords diversity microbial models nonlinear
#' @examples
#' breakaway(apples)
#' breakaway(apples,plot=FALSE,print=FALSE,answers=TRUE) 
#' 
breakaway <- function(my_data, print=TRUE, plot=TRUE, answers=FALSE, force=FALSE, useAll=FALSE) {
  
  ## read in data
  if( !(is.matrix(my_data) || is.data.frame(my_data))) {
    filename <- my_data
    ext <- substr(filename, nchar(filename)-2, nchar(filename))
    if (ext == "csv") {
      my_data <- read.table(file=filename, header=0,sep=",")
      if( my_data[1,1] !=1) my_data <- read.table(filename, header=1,sep=",")
    } else if (ext == "txt") {
      my_data <- read.table(file=filename, header=0)
    } else stop("Please input your data as a txt or csv file,
               or as an R dataframe or matrix.")
  }
  
  ## if tables() is used to create the frequency tables, the frequency index column is usually a factor, so fix this here
  if ( is.factor(my_data[,1]) ) {
    fs <- as.numeric(as.character(my_data[,1]))
    my_data <- cbind(fs,my_data[,2])
    my_data <- my_data[my_data[,1]!=0,]
  }
  
  if(length(my_data) <= 1) {
    stop("Input data is of length 1 or 0. Huh?")
  }
  
  if (my_data[1,1]!=1 || my_data[1,2]==0) {
    stop("You don't have an observed singleton count.\n breakaway isn't built for that data structure.\n")
  } 
  
  my_data <- my_data[!(my_data[,2]==0 | is.na(my_data[,2])),]
  orig_my_data <- my_data
  n <- sum(orig_my_data[,2])
  f1 <- my_data[1,2]
  
  if (useAll) {
    warning("The useAll option is in beta! Use at your own risk.\nPlease report any issues via github or to Amy directly.")
    ## Use all contiguous frequencies. Default is to cut off at the first break in frequencies.
    
    frequency_vector <- my_data[ , 1]
    frequency_counts <- my_data[ , 2]
    
    frequency_ratios <- frequency_counts[-1] / frequency_counts[-length(frequency_counts)]
    frequency_vector_differences <- frequency_vector[-1] - frequency_vector[-length(frequency_vector)]
    
    xs <- frequency_vector[frequency_vector_differences == 1]
    xs <- xs[-length(xs)] # need to take off the last 1 because we lose 1 with ratios
    xbar <- mean(xs)
    
    y <- frequency_ratios[frequency_vector_differences == 1]
    ys <- (xs+1)*y
    
    lhs <- list("x"=xs-xbar,"y"=y)
    cutoff <- length(y)
  } else {
    ## breakaway's default is to cut off at the first break in frequencies
    
    ## Horrible way to do this. What was 22 y.o. Amy thinking?!
    ## finds the break in contiguity
    cutoff <- ifelse(is.na(which(my_data[-1,1]-my_data[-length(my_data[,1]),1]>1)[1]),length(my_data[,1]),which(my_data[-1,1]-my_data[-length(my_data[,1]),1]>1)[1])
    
    my_data <- my_data[1:cutoff,]
    ys <- (my_data[1:(cutoff-1),1]+1)*my_data[2:cutoff,2]/my_data[1:(cutoff-1),2]
    xs <- 1:(cutoff-1)
    xbar <- mean(xs)
    lhs <- list("x"=xs-xbar,"y"=my_data[2:cutoff,2]/my_data[1:(cutoff-1),2])
  }
  
  if (cutoff < 6) { ## check for unusual data structures
    stop("You don't have enough contiguous frequencies.\nbreakaway needs at least 6!\n")
  } 
  
  if ( !force && ( (my_data[cutoff,2]/my_data[cutoff-1,2])>10 )) {
    stop("\tIt looks like you've concatenated some of your data! Please truncate and try again.\n")
  } 
  
  weights_inv <- 1/xs
  run <- .minibreak(lhs,xs,ys,my_data,weights_inv)
  result <- list()
  choice <- list()
  
  ### If no models converged, use the WLRM
  if (sum(as.numeric(run$useful[,5]))==0) {
    choice$outcome <- 0
  
    if (useAll) {
      stop("No models converged. Your dataset may be too small for breakaway's model set, or too irregular.")
    } else {
      if(print) cat("No breakaway models converged.")
      weights_trans <- (1/my_data[-1,2]+1/my_data[-cutoff,2])^-1
      lm1 <- lm(log(ys)~xs,weights=weights_trans)
      b0_hat <- summary(lm1)$coef[1,1]; b0_se <- summary(lm1)$coef[1,2]
      f0 <- f1*exp(-b0_hat)
      diversity <- f0 + n
      f0_se <- sqrt( (exp(-b0_hat))^2*f1*(b0_se^2*f1+1) )  #consistent with rbbo
      diversity_se <- sqrt(f0_se^2+n*f0/(n+f0))
      
      if(print) {
        cat("################## breakaway ##################\n")
        cat("\tThe best estimate of total diversity is", round(diversity),
            "\n \t with std error",round(diversity_se),"\n")
        cat("\tThe model employed was the WLRM\n")
      }
      if(answers) {
        result$name <- "WLRM"
        result$para <- summary(lm1)$coef[,1:2]
        result$est <- diversity
        result$seest <- as.vector(diversity_se)
        result$full <- lm1
        d <- exp(1.96*sqrt(log(1+result$seest^2/f0)))
        result$ci <- c(n+f0/d,n+f0*d)
        return(result)
      }
      
      result$code <- 1
    }
    
  } else { 
    ## Otherwise, YAY! Something worked for 1/x weighting
    choice$outcome <- 1
    choice$model <- rownames(run$useful)[min(which(run$useful[,5]==1))] #pick smallest
    choice$full <-  run[[noquote(choice$model)]]
    
    oldest <- run$useful[run$useful[,5]==1,1][[1]]
    est <- 0
    its <- 0
    while ( choice$outcome & abs(oldest-est)>1 & its < 30) {
      oldest <- est
      C <- round(n+oldest,0)
      unscaledprobs <- c(1,cumprod(fitted(choice$full)))
      p <- unscaledprobs/sum(unscaledprobs)
      as <- p[-1]^2/p[-cutoff]^3 * (1-exp(-C*p[-cutoff]))^3/(1-exp(-C*p[-1]))^2 * (1-C*p[-cutoff]/(exp(C*p[-cutoff])-1))
      bs <- p[-1]/p[-cutoff]^2 * (1-exp(-C*p[-cutoff]))^2/(1-exp(-C*p[-1])) * (1-C*p[-1]/(exp(C*p[-1])-1))
      ratiovars <- (as + bs)/C
      
      if(its==0) {
        weights1 <- 1/ratiovars
      }
      
      run <- try ( .minibreak(lhs,xs,ys,my_data,1/ratiovars), silent = 1)
      
      if ( class(run) == "try-error") {
        ratiovars <- (p[-1]^2/p[-cutoff]^3 + p[-1]/p[-cutoff]^2)/C
        run <- try ( .minibreak(lhs,xs,ys,my_data,1/ratiovars), silent = 1)
        if ( class(run) == "try-error") {
          if(print) {print("Numerical errors result in non-convergence") }
        }
      }
      
      choice <- list()
      if ( class(run)=="try-error" ) {
        choice$outcome <- 0
      } else if (sum(as.numeric(run$useful[,5]))==0) {
        choice$outcome <- 0
      } else {
        choice$outcome <- 1
        choice$model <- rownames(run$useful)[min(which(run$useful[,5]==1))]
        choice$full <-  run[[noquote(choice$model)]]
        
        est <- run$useful[run$useful[,5]==1,1][[1]]
        its <- its + 1
        result$code <- 2
      }
    }
    if( !choice$outcome) {
      if(print) cat("We used 1/x weighting. \n")
      run <- .minibreak(lhs,xs,ys,my_data,weights_inv)
      choice$outcome <- 1
      choice$model <- rownames(run$useful)[min(which(run$useful[,5]==1))]
      choice$full <-  run[[noquote(choice$model)]]
      result$code <- 3
    }
    
    if(choice$model=="model_1_0") {
      b0_hat <- coef(choice$full)[1]
      b0_var <- vcov(choice$full)[1,1]
    } else {
      effective_coef <- c(coef(choice$full),rep(0,9-length(coef(choice$full))))
      b <- effective_coef[1]-effective_coef[2]*xbar+effective_coef[4]*xbar^2-effective_coef[6]*xbar^3+effective_coef[8]*xbar^4
      a <- 1-effective_coef[3]*xbar+effective_coef[5]*xbar^2-effective_coef[7]*xbar^3+effective_coef[9]*xbar^4
      
      nabla <- c(1/a, -xbar/a, b*xbar/a^2, xbar^2/a, -b*xbar^2/a^2, -xbar^3/a,b*xbar^3/a^2, xbar^4/a, -b*xbar^4/a^2)
      nabla <- nabla[1:length(coef(choice$full))]
      
      b0_hat <- b/a
      b0_var <- t(nabla)%*%vcov(choice$full)%*%nabla
    }
    
    f0 <- run$useful[rownames(run$useful)==choice$model,1][[1]]
    f0_var <- f1*b0_hat^-2*(1-f1/n+f1*b0_hat^-2*b0_var) #1st order d.m.
    
    diversity <- f0 + n
    diversity_se <- sqrt(n*f0/diversity + f0_var)
    
    if (choice$model == "model_1_0") the_function <- "f_{x+1}/f_{x} ~ (beta0+beta1*x)/(1+x)"
    if (choice$model == "model_1_1") the_function <- "f_{x+1}/f_{x} ~ (beta0+beta1*(x-xbar))/(1+alpha1*(x-xbar))"
    if (choice$model == "model_2_1") the_function <- "f_{x+1}/f_{x} ~ (beta0+beta1*(x-xbar)+beta2*(x-xbar)^2)/(1+alpha1*(x-xbar))"
    if (choice$model == "model_2_2") the_function <- "f_{x+1}/f_{x} ~ (beta0+beta1*(x-xbar)+beta2*(x-xbar)^2)/(1+alpha1*(x-xbar)+alpha2)"
    if (choice$model == "model_3_2") the_function <- "f_{x+1}/f_{x} ~ (beta0+beta1*(x-xbar)+beta2*(x-xbar)^2+beta3*(x-xbar)^3)/(1+alpha1*(x-xbar)+alpha2*(x-xbar)^2)"
    if (choice$model == "model_3_3") the_function <- "f_{x+1}/f_{x} ~ (beta0+beta1*(x-xbar)+beta2*(x-xbar)^2+beta3*(x-xbar)^3)/(1+alpha1*(x-xbar)+alpha2*(x-xbar)^2+alpha3*(x-xbar)^3)"
    if (choice$model == "model_4_3") the_function <- "f_{x+1}/f_{x} ~ (beta0+beta1*(x-xbar)+beta2*(x-xbar)^2+beta3*(x-xbar)^3+beta4*(x-xbar)^4)/(1+alpha1*(x-xbar)+alpha2*(x-xbar)^2+alpha3*(x-xbar)^3)"
    if (choice$model == "model_4_4") the_function <- "f_{x+1}/f_{x} ~ (beta0+beta1*(x-xbar)+beta2*(x-xbar)^2+beta3*(x-xbar)^3+beta4*(x-xbar)^4)/(1+alpha1*(x-xbar)+alpha2*(x-xbar)^2+alpha3*(x-xbar)^3+alpha4*(x-xbar)^4)"
    
    parameter_table <-  coef(summary(choice$full))[,1:2]
    colnames(parameter_table) <- c("Coef estimates","Coef std errors")
    
    if(print) {
      cat("################## breakaway ##################\n")
      cat("\tThe best estimate of total diversity is", round(diversity),
          "\n \t with std error",round(diversity_se),"\n")
      cat("\tThe model employed was", choice$model,"\n")
      cat("\tThe function selected was\n\t ", the_function,"\n")
      print(parameter_table)
      cat("xbar\t\t\t",xbar)
    }
    
    if(plot)  {
      yhats <- fitted(choice$full)
      par(new=FALSE)
      plot(xs,lhs$y,xlim=c(0,max(xs)+1),ylim=c(min(lhs$y,yhats),max(lhs$y)*1.05),
           ylab="f(x+1)/f(x)",xlab="x",main="Plot of ratios and fitted values under model");
      points(xs,yhats,pch=18)
      points(0,b0_hat,col="red",pch=18)
      legend(0,max(lhs$y),c("Fitted values", "Prediction"),pch=c(18,18),col=c("black", "red"),cex=0.8,bty = "n")
    }
    
    if(answers) {
      result$name <- choice$model
      result$para <- parameter_table
      result$est <- diversity
      result$seest <- as.vector(diversity_se)
      result$full <- choice$full
      d <- exp(1.96*sqrt(log(1+result$seest^2/f0)))
      result$ci <- c(n+f0/d,n+f0*d)
      return(result)
    }
  }
}

.minibreak <- function(lhs, xs, ys, my_data, myweights) {
  structure_1_0 <- function(x,beta0,beta1) (beta0+beta1*x)/(1+x)
  structure_1_1 <- function(x,beta0,beta1,alpha1) (beta0+beta1*x)/(1+alpha1*x)
  structure_2_1 <- function(x,beta0,beta1,alpha1,beta2) (beta0+beta1*x+beta2*x^2)/(1+alpha1*x)
  structure_2_2 <- function(x,beta0,beta1,alpha1,beta2,alpha2) (beta0+beta1*x+beta2*x^2)/(1+alpha1*x+alpha2*x^2)
  structure_3_2 <- function(x,beta0,beta1,alpha1,beta2,alpha2,beta3) (beta0+beta1*x+beta2*x^2+beta3*x^3)/(1+alpha1*x+alpha2*x^2)
  structure_3_3 <- function(x,beta0,beta1,alpha1,beta2,alpha2,beta3,alpha3) (beta0+beta1*x+beta2*x^2+beta3*x^3)/(1+alpha1*x+alpha2*x^2+alpha3*x^3)
  structure_4_3 <- function(x,beta0,beta1,alpha1,beta2,alpha2,beta3,alpha3,beta4) (beta0+beta1*x+beta2*x^2+beta3*x^3+beta4*x^4)/(1+alpha1*x+alpha2*x^2+alpha3*x^3)
  structure_4_4 <- function(x,beta0,beta1,alpha1,beta2,alpha2,beta3,alpha3,beta4,alpha4) (beta0+beta1*x+beta2*x^2+beta3*x^3+beta4*x^4)/(1+alpha1*x+alpha2*x^2+alpha3*x^3+alpha4*x^4)
  
  all_ests <- list()
  
  model_1_0 <-  lm(ys~xs,weights=myweights)
  start_1_0 <- list(beta0 = model_1_0$coef[1], beta1 = model_1_0$coef[2])
  
  lhs2 <- lhs; lhs2$x <- xs # otherwise singularity
  result_1_0 <- try( nls( lhs2$y ~ structure_1_0(x,beta0,beta1), data = lhs2, start_1_0, weights=myweights), silent=T )
  if(class(result_1_0) == "try-error") {
    coef_1_0 <- c(start_1_0$beta0,start_1_0$beta1)
  } else {
    coef_1_0 <- coef(result_1_0)
  }
  
  beta0_tilde <- coef_1_0[1]
  beta1_tilde <- coef_1_0[2]
  xbar <- mean(xs)
  start_1_1 <- list(beta0 = (beta0_tilde+beta1_tilde*xbar)/(1+xbar), beta1 = beta1_tilde/(1+xbar), alpha1 = 1/(1+xbar))
  result_1_1 <- try( nls( lhs$y ~ structure_1_1(x,beta0,beta1,alpha1), data = lhs, start_1_1, weights=myweights), silent=T )
  
  if(class(result_1_1) == "try-error") {
    temp_lm <- lm(lhs$y~lhs$x+I((lhs$x)^2),weights=myweights)$coef
    start_2_1 <- list(beta0 = temp_lm[1], beta1 = temp_lm[2], alpha1 = 0, beta2 = temp_lm[3])
  } else {
    start_2_1 <- as.list(c(summary(result_1_1)$coef[,1],"beta2"=0))
  }
  result_2_1 <- try(nls( lhs$y ~ structure_2_1(x,beta0,beta1,alpha1,beta2), data = lhs, start_2_1, weights=myweights),silent=T)
  
  if(class(result_2_1) == "try-error") {
    temp_lm <- lm(lhs$y~lhs$x+I((lhs$x)^2),weights=myweights)$coef
    start_2_2 <- list(beta0 = temp_lm[1], temp_lm[2], alpha1 = 0, beta2 = temp_lm[3], alpha2 = 0)
  } else start_2_2 <- as.list(c(summary(result_2_1)$coef[,1],"alpha2"=0))
  result_2_2 <- try(nls( lhs$y ~ structure_2_2(x,beta0,beta1,alpha1,beta2,alpha2), data = lhs, start_2_2, weights=myweights), silent=T)
  
  if(class(result_2_2) == "try-error") {
    temp_lm <- lm(lhs$y~lhs$x+I((lhs$x)^2)+I((lhs$x)^3),weights=myweights)$coef
    start_3_2 <- list(beta0 = temp_lm[1], beta1 = temp_lm[2], alpha1 = 0, beta2 = temp_lm[3], alpha2 = 0, beta3 = temp_lm[4])
  } else start_3_2 <- as.list(c(summary(result_2_2)$coef[,1],"beta3"=0))
  result_3_2 <- try( nls( lhs$y ~ structure_3_2(x,beta0,beta1,alpha1,beta2,alpha2,beta3), data = lhs, start_3_2, weights=myweights), silent=T)
  
  if(class(result_3_2) == "try-error") {
    temp_lm <- lm(lhs$y~lhs$x+I((lhs$x)^2)+I((lhs$x)^3), weights=myweights)$coef
    start_3_3 <- list(beta0 = temp_lm[1], beta1 = temp_lm[2], alpha1 = 0, beta2 = temp_lm[3], alpha2 = 0, beta3 = temp_lm[4], alpha3 = 0)
  } else start_3_3 <- as.list(c(summary(result_3_2)$coef[,1],"alpha3"=0))
  result_3_3 <- try( nls( lhs$y ~ structure_3_3(x,beta0,beta1,alpha1,beta2,alpha2,beta3,alpha3), data = lhs, start_3_3, weights=myweights), silent=T)
  
  if(class(result_3_3) == "try-error") {
    temp_lm <- lm(lhs$y~lhs$x+I((lhs$x)^2)+I((lhs$x)^3)+I((lhs$x)^4),weights=myweights)$coef
    start_4_3 <- list(beta0 = temp_lm[1], beta1 = temp_lm[2], alpha1 = 0, beta2 = temp_lm[3], alpha2 = 0,beta3 = temp_lm[4], alpha3 = 0,beta4 = temp_lm[5])
  } else start_4_3 <- as.list(c(summary(result_3_3)$coef[,1],"beta4"=0))
  result_4_3 <- try( nls( lhs$y ~ structure_4_3(x,beta0,beta1,alpha1,beta2,alpha2,beta3,alpha3,beta4), data = lhs, start_4_3, weights=myweights), silent=T)
  
  if(class(result_4_3) == "try-error") {
    temp_lm <- lm(lhs$y~lhs$x+I((lhs$x)^2)+I((lhs$x)^3)+I((lhs$x)^4),weights=myweights)$coef
    start_4_4 <- list(beta0 = temp_lm[1], beta1 = temp_lm[2], alpha1 = 0, beta2 = temp_lm[3], alpha2 = 0,
                      beta3 = temp_lm[4], alpha3 = 0, beta4 = temp_lm[5], alpha4 = 0)
  } else start_4_4 <- as.list(c(summary(result_4_3)$coef[,1],"alpha4"=0))
  result_4_4 <- try( nls( lhs$y ~ structure_4_4(x,beta0,beta1,alpha1,beta2,alpha2,beta3,alpha3,beta4,alpha4), data = lhs, start_4_4, weights=myweights), silent=T)
  
  info <- list("model_1_0"=result_1_0,"model_1_1"=result_1_1,
               "model_2_1"=result_2_1,"model_2_2"=result_2_2,"model_3_2"=result_3_2,
               "model_3_3"=result_3_3,"model_4_3"=result_4_3,"model_4_4"=result_4_4)
  
  converged <- lapply(info,class)=="nls"
  
  workinginfo <- info[converged]
  
  b0est <- lapply(workinginfo,predict,list(x=-xbar))
  if (length(b0est$model_1_0)!=0) b0est$model_1_0 <- predict(workinginfo$model_1_0,list(x=0))
  f0est <- lapply(b0est,function(x) my_data[1,2]/x)
  
  rootsokay <- as.logical(lapply(workinginfo,.rootcheck,lhs,nof1=FALSE))
  
  sqerrors <- lapply(workinginfo,.sqerror,lhs)
  residses <- lapply(workinginfo,.residse)
  
  useable <- rootsokay & (f0est>0)
  
  workinginfo$useful <- cbind(f0est,rootsokay,sqerrors,residses,useable)
  return(workinginfo)
}

.rootcheck <- function(model,lhs,nof1=FALSE) {
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

.sqerror <- function(model, lhs) {
  fits <- fitted(model)
  return(sum(((lhs$y-fits)^2)/fits))
}

.residse <- function(model) {
  return(summary(model)$sigma)
}
