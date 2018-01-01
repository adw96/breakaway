#' Species richness estimation with breakaway
#' 
#' This function implements the species richness estimation procedure outlined
#' in Willis & Bunge (2015). The diversity estimate, standard error, estimated
#' model coefficients, model details and plot of the fitted model are returned.
#' 
#' 
#' @param my_data The sample frequency count table for the population of
#' interest. The first row must correspond to the singletons. Acceptable
#' formats include a matrix, data frame, or file path (csv or txt). The
#' standard frequency count table format is used: two columns, the first of
#' which contains the frequency of interest (eg. 1 for singletons, species
#' observed once, 2 for doubletons, species observed twice, etc.) and the
#' second of which contains the number of species observed this many times.
#' Frequencies (first column) should be ordered least to greatest. At least 6
#' contiguous frequencies are necessary. Do not concatenate large frequencies.
#' See dataset apples for sample formatting.
#' @param output Logical: whether the results should be printed to screen. If
#' \samp{FALSE}, answers should be set to \samp{TRUE} so that results will be
#' returned.
#' @param plot Logical: whether the data and model fit should be plotted.
#' @param answers Logical: whether the function should return an argument. If
#' \samp{FALSE}, output should be set to \samp{TRUE}.
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
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' breakaway(apples)
#' breakaway(apples, plot = FALSE, output = FALSE, answers = TRUE)
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' @export breakaway
breakaway <- function(my_data, output=TRUE, plot=FALSE, answers=FALSE, print = NULL) {
  
  if (!is.null(print)) {
    warning("Argument 'print' deprecated")
  }
  orig_my_data <- check_format(my_data)
  
  if (my_data[1,1]!=1 || my_data[1,2]==0) {
    stop("You don't have an observed singleton count.\n breakaway isn't built for that data structure.\n")
  } 
  
  n <- sum(orig_my_data[,2])
  f1 <- my_data[1,2]
  
  ## breakaway's default is to cut off at the first break in frequencies
  
  ## finds the break in contiguity
  iss <- my_data[, 1]
  fis <- my_data[, 2]
  length_fis <- length(fis)
  
  breaks <- which(iss[-1] - iss[-length_fis] > 1)
  cutoff <- ifelse(is.na(breaks[1]), length_fis, breaks[1])
  
  ## set up structures
  my_data <- my_data[1:cutoff,]
  ratios <- my_data[-1, 2]/my_data[-cutoff, 2]
  xs <- 1:(cutoff-1)
  ys <- (xs+1)*ratios
  xbar <- mean(xs)
  lhs <- list("x"=xs-xbar,"y"=ratios)
  
  
  if (cutoff < 6) { ## check for unusual data structures
    message("You don't have enough contiguous frequencies
            to use breakaway. Using Chao-Bunge instead...\n")
    
    result_temporary <- chao_bunge(my_data, answers = TRUE, output = FALSE)
    
    if(output) {
      cat("################## breakaway ##################\n")
      cat("\tThe best estimate of total diversity is", round(result_temporary$est),
          "\n \t with std error",round(result_temporary$seest),"\n")
      cat("\tThe model employed was Chao-Bunge (Negative binomial) \n")
    }
    
    if(answers) {
      return(result_temporary)
    }
    
  } else {
    
    
    weights_inv <- 1/xs
    run <- minibreak_all(lhs,xs,ys,my_data,weights_inv, withf1=TRUE)
    result <- list()
    choice <- list()
    
    ### If no models converged, use the WLRM
    if (sum(as.numeric(run$useful[,5]))==0) {
      choice$outcome <- 0
      
      if(output) cat("No breakaway models converged.")
      weights_trans <- (1/my_data[-1,2]+1/my_data[-cutoff,2])^-1
      lm1 <- lm(log(ys)~xs,weights=weights_trans)
      b0_hat <- summary(lm1)$coef[1,1]; b0_se <- summary(lm1)$coef[1,2]
      f0 <- f1*exp(-b0_hat)
      diversity <- f0 + n
      f0_se <- sqrt( (exp(-b0_hat))^2*f1*(b0_se^2*f1+1) )  #consistent with rbbo
      diversity_se <- sqrt(f0_se^2+n*f0/(n+f0))
      
      if(output) {
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
        
        run <- try ( minibreak_all(lhs,xs,ys,my_data,1/ratiovars, withf1 = TRUE), silent = 1)
        
        if ( class(run) == "try-error") {
          ratiovars <- (p[-1]^2/p[-cutoff]^3 + p[-1]/p[-cutoff]^2)/C
          run <- try ( minibreak_all(lhs,xs,ys,my_data,1/ratiovars, withf1 = TRUE), silent = 1)
          if ( class(run) == "try-error") {
            if(output) {print("Numerical errors result in non-convergence") }
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
        if(output) cat("We used 1/x weighting. \n")
        run <- minibreak_all(lhs,xs,ys,my_data,weights_inv, withf1 = TRUE)
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
      
      if(output) {
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
}

#' @export richness_clean
richness_clean <- function(count_table, FUN = breakaway) {
  richness_output <- try(FUN(count_table, answers = T, output = F, plot = F))
  if (class(richness_output) != "try-error") {
    data.frame("index" = "Richness", 
               "estimate" = richness_output$est, 
               "standard_error" = richness_output$seest, 
               "lower" = richness_output$ci[1], 
               "upper" = richness_output$ci[2])
  } else {
    data.frame("index" = "Richness", 
               "estimate" = NA, 
               "standard_error" = NA, 
               "lower" = sum(count_table[,2]), 
               "upper" = Inf)
  }
}

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















































#' species richness estimation without singletons
#' 
#' This function permits estimation of total diversity based on a sample
#' frequency count table. Unlike \code{\link{breakaway}}, it does not require
#' an input for the number of species observed once, making it an excellent
#' exploratory tool for microbial ecologists who believe that their sample may
#' contain spurious singletons. The underlying estimation procedure is similar
#' to that of \code{\link{breakaway}} and is outlined in Willis & Bunge (2014).
#' The diversity estimate, standard error, estimated model coefficients and
#' plot of the fitted model are returned.
#' 
#' 
#' @param my_data The sample frequency count table for the population of
#' interest. The first row must correspond to the doubletons. Acceptable
#' formats include a matrix, data frame, or file path (csv or txt). The
#' standard frequency count table format is used: two columns, the first of
#' which contains the frequency of interest (eg. 1 for singletons, species
#' observed once; 2 for doubletons, species observed twice, etc.) and the
#' second of which contains the number of species observed this many times.
#' Frequencies (first column) should be ordered least to greatest. At least 6
#' contiguous frequencies are necessary. Do not concatenate large frequencies.
#' See dataset \code{\link{apples}} for sample formatting.
#' @param output Logical: whether the results should be printed to screen. If
#' \samp{FALSE}, \samp{answers} should be set to \samp{TRUE} so that results
#' will be returned.
#' @param plot Logical: whether the data and model fit should be plotted.
#' @param answers Logical: whether the function should return an argument. If
#' \samp{FALSE}, \samp{output} should be set to \samp{TRUE}.
#' @param force Logical: force the procedure to run in the presence of
#' frequency count concatenation. A basic check procedure confirms that the
#' user has not appeared to concatenate multiple upper frequencies.
#' \samp{force=TRUE} will force the procedure to fit models in the presence of
#' this. \samp{breakaway_nof1}'s diversity estimates cannot be considered
#' reliable in this case.
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
#' degree and the second integer represents the denominator polynomial degree.
#' See Willis & Bunge (2014) for details.  } \item{para}{ Estimated model
#' parameters and standard errors.  } \item{est}{ The estimate of total
#' (observed plus unobserved) diversity.  } \item{seest}{ The standard error in
#' the diversity estimate.  } \item{full}{ The chosen nonlinear model for
#' frequency ratios.  }
#' @note It is common for microbial ecologists to believe that their dataset
#' contains false diversity. This often arises because sequencing errors result
#' in classifying already observed organisms as new organisms.
#' \samp{breakaway_nof1} was developed as an exploratory tool in this case.
#' Practitioners can run \samp{breakaway} on their dataset including the
#' singletons, and \samp{breakaway_nof1} on their dataset excluding the
#' singletons, and assess if the estimated levels of diversity are very
#' different. Great disparity may provide evidence of an inflated singleton
#' count, or at the very least, that \samp{breakaway} is especially sensitive
#' to the number of rare species observed. Note that \samp{breakaway_nof1} may
#' be less stable than \samp{breakaway} due to predicting based on a reduced
#' dataset, and have greater standard errors.
#' @author Amy Willis
#' @seealso \code{\link{breakaway}}; \code{\link{apples}}
#' @references Willis, A. (2015). Species richness estimation with high
#' diversity but spurious singletons. \emph{Under review.}
#' 
#' Willis, A. and Bunge, J. (2015). Estimating diversity via frequency ratios.
#' \emph{Biometrics.}
#' @keywords diversity error microbial models nonlinear
#' @examples
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' breakaway_nof1(apples[-1, ])
#' breakaway_nof1(apples[-1, ], plot = FALSE, output = FALSE, answers = TRUE)
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' @export breakaway_nof1
breakaway_nof1 <- function(my_data, output=TRUE, plot=TRUE, answers=FALSE, force=FALSE) {
  
  my_data <- check_format(my_data)
  
  
  if (my_data[1,1]!=2 || my_data[1,2]==0) {
    if(output) cat("breakaway_nof1 is for when you have no singleton count.\nYou need a leading doubleton count!\n")
  } else {
    my_data <- my_data[!(my_data[,2]==0 | is.na(my_data[,2])),]
    orig_my_data <- my_data
    n <- sum(orig_my_data[,2])
    f2 <- my_data[1,2]
    
    cutoff <- ifelse(is.na(which(my_data[-1,1]-my_data[-length(my_data[,1]),1]>1)[1]),length(my_data[,1]),which(my_data[-1,1]-my_data[-length(my_data[,1]),1]>1)[1])
    my_data <- my_data[1:cutoff,]
    ys <- (my_data[1:(cutoff-1),1]+1)*my_data[2:cutoff,2]/my_data[1:(cutoff-1),2]
    xs <- 2:cutoff
    xbar <- mean(c(1,xs))
    lhs <- list("x"=xs-xbar,"y"=my_data[2:cutoff,2]/my_data[1:(cutoff-1),2])
    
    if ( cutoff < 5) { ### check for unusual data structures
      if(output) cat("You don't have enough contiguous frequencies.\n breakaway needs at least 6!\n")
    } else if ((force==FALSE) && ( (my_data[cutoff,2]/my_data[cutoff-1,2])>10 )) {
      cat("\tIt looks like you've concatenated some of your my_data!\n Please truncate and try again.\n")
    } else {
      weights_inv <- 1/xs
      run <- minibreak_all(lhs,xs,ys,my_data,weights_inv, withf1 = FALSE)
      result <- list()
      choice <- list()
      
      if (sum(as.numeric(run$useful[,5]))==0) {
        choice$outcome <- 0
        if(output) cat("No breakaway models converged.")
        weights_trans <- (1/my_data[-1,2]+1/my_data[-cutoff,2])^-1
        lm1 <- lm(log(ys)~xs,weights=weights_trans)
        b0_hat <- summary(lm1)$coef[1,1]
        b0_se <- summary(lm1)$coef[1,2]
        b1_hat <- summary(lm1)$coef[2,1]
        b1_se <- summary(lm1)$coef[2,2]
        
        f1_pred <- 2*f2*exp(-(b0_hat+b1_hat))
        f0_pred <- f1_pred*exp(-b0_hat)
        diversity <- f0_pred + f1_pred + n
        
        covmatrix <- matrix(c(f2*(1-f2/n),rep(0,8)),nrow=3)
        covmatrix[2:3,2:3] <- vcov(lm1)
        
        derivs <- c(exp(-(b0_hat+b1_hat))*(1+exp(-b1_hat)),
                    f2*exp(-b1_hat)*(1+exp(-b1_hat))*exp(-b0_hat),
                    f2*exp(-b0_hat)*(-exp(-b1_hat)-2*exp(-b1_hat)))
        
        f0plusf1_var <- 4*t(derivs)%*%covmatrix%*%(derivs)
        
        diversity_se <- sqrt(f0plusf1_var+n*(f0_pred+f1_pred)/(n+f0_pred+f1_pred))
        
        if(output) {
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
          d <- exp(1.96*sqrt(log(1+result$seest^2/f0_pred)))
          result$ci <- c(n+f0_pred/d,n+f0_pred*d)
          return(result)
        }
        
        result$code <- 1
      } else { # something worked for 1/x weighting
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
          
          run <- try ( minibreak_all(lhs, xs, ys, my_data, 1/ratiovars, withf1=FALSE), silent = 1)
          
          if ( class(run) == "try-error") {
            ratiovars <- (p[-1]^2/p[-cutoff]^3 + p[-1]/p[-cutoff]^2)/C
            run <- try ( minibreak_all(lhs,xs,ys,my_data,1/ratiovars, withf1=FALSE), silent = 1)
            if ( class(run) == "try-error") {
              if(output) {print("Numerical errors result in non-convergence") }
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
          if(output) cat("Iterative reweighting didn't produce any outcomes after the first iteration, so we use 1/x\n")
          run <- minibreak_all(lhs,xs,ys,my_data,weights_inv, withf1 = FALSE)
          choice$outcome <- 1
          choice$model <- rownames(run$useful)[min(which(run$useful[,5]==1))]
          choice$full <-  run[[noquote(choice$model)]]
          result$code <- 3
        }
        
        effective_coef <- c(coef(choice$full),rep(0,9-length(coef(choice$full))))
        if (choice$model=="model_1_0") {
          effective_coef[3] <- 1
        }
        
        b <- effective_coef[1]-effective_coef[2]*xbar+effective_coef[4]*xbar^2-effective_coef[6]*xbar^3+effective_coef[8]*xbar^4
        a <- 1-effective_coef[3]*xbar+effective_coef[5]*xbar^2-effective_coef[7]*xbar^3+effective_coef[9]*xbar^4
        nabla_b <- c(1/a, -xbar/a, b*xbar/a^2, xbar^2/a, -b*xbar^2/a^2, -xbar^3/a, b*xbar^3/a^2, xbar^4/a, -b*xbar^4/a^2)
        nabla_b <- nabla_b[1: length(coef(choice$full))]
        
        b1 <- effective_coef[1]+effective_coef[2]*(1-xbar)+effective_coef[4]*(1-xbar)^2+effective_coef[6]*(1-xbar)^3+effective_coef[8]*(1-xbar)^4
        a1 <- 1+effective_coef[3]*(1-xbar)+effective_coef[5]*(1-xbar)^2+effective_coef[7]*(1-xbar)^3+effective_coef[9]*(1-xbar)^4
        nabla_c <- c(1/a1, (1-xbar)/a1, -b1*(1-xbar)/a1^2, (1-xbar)^2/a1, -b1*(1-xbar)^2/a1^2, (1-xbar)^3/a1,-b1*(1-xbar)^3/a1^2,(1-xbar)^4/a1,-b1*(1-xbar)^4/a1^2)
        nabla_c <- nabla_c[1: length(coef(choice$full))]
        
        b0_hat <- b/a
        c0_hat <- b1/a1
        
        if (choice$model=="model_1_0") {
          nabla_b <- c(nabla_b,0)
          new_cov_b <- matrix(0,nrow=3,ncol=3)
          new_cov_b[1:2,1:2] <- vcov(choice$full)
          b0_var <- t(nabla_b)%*%new_cov_b%*%nabla_b
          
          nabla_c <- c(nabla_c,0)
          new_cov_c <- matrix(0,nrow=3,ncol=3)
          new_cov_c[1:2,1:2] <- vcov(choice$full)
          c0_var <- t(nabla_c)%*%new_cov_c%*%nabla_c
        } else {
          b0_var <- t(nabla_b)%*%vcov(choice$full)%*%nabla_b
          c0_var <- t(nabla_c)%*%vcov(choice$full)%*%nabla_c
        }
        
        f1_pred <- f2/c0_hat
        f0_pred <- f1_pred/b0_hat
        diversity <- f0_pred + f1_pred + n
        
        covmatrix <- diag(c(f2*(1-f2/diversity),b0_var,c0_var))
        if (choice$model == "model_1_0") {
          cov_b0_c0 <- 0.5*sum(vcov(choice$full)[,1])
          covmatrix[2,3] <- cov_b0_c0
          covmatrix[3,2] <- cov_b0_c0
        } else {
          covmatrix[2,3] <- b0_var/a1
          covmatrix[3,2] <- b0_var/a1
        }
        derivs <- matrix(c(b0_hat^-1*c0_hat^-1,
                           -f2*b0_hat^-2*c0_hat^-1,
                           -f2*b0_hat^-1*c0_hat^-2,c0_hat^-1,0,-f2*b0_hat^-2),byrow=TRUE,ncol=3,nrow=2)
        f0plusf1_var <- derivs%*%covmatrix%*%t(derivs)
        
        var_est <- sum(f0plusf1_var) + n*(f0_pred+f1_pred)/diversity - 2*f0_pred*n/diversity - 2*f1_pred*n/diversity
        diversity_se <- ifelse(var_est<0, 0, sqrt(var_est)) ## Piece 6
        
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
        
        if(output) {
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
          points(c(0,1),c(b0_hat,c0_hat),col="red",pch=18)
          legend(0,max(lhs$y),c("Fitted values", "Prediction"),pch=c(18,18),col=c("black", "red"),cex=0.8,bty = "n")
        }
        
        if(answers) {
          result$name <- choice$model
          result$para <- parameter_table
          result$est <- diversity
          result$seest <- as.vector(diversity_se)
          result$full <- choice$full
          d <- exp(1.96*sqrt(log(1+result$seest^2/f0_pred)))
          result$ci <- c(n+f0_pred/d,n+f0_pred*d)
          return(result)
        }
      }
    }
  }
}

minibreak_all <- function(lhs, xs, ys, my_data, myweights, withf1 = NULL) {
  
  if (is.null(withf1)) stop("Empty argument `withf1`")
  
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
  
  xbar <- ifelse(withf1, mean(xs), mean(c(1,xs)))
  
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
  
  if (withf1) {
    
    f0est <- lapply(b0est,function(x) my_data[1,2]/x)
    
  } else {
    
    firstratioest <- lapply(workinginfo,predict,list(x=-xbar+1))
    if (length(firstratioest$model_1_0)!=0) firstratioest$model_1_0 <- predict(workinginfo$model_1_0,list(x=1))
    
    f1est <- lapply(firstratioest,function(x) my_data[1,2]/x)
    
    f0est <- as.numeric(f1est)/as.numeric(b0est)
    
    
  }
  
  rootsokay <- as.logical(lapply(workinginfo,rootcheck,lhs,nof1=!withf1))
  
  sqerrors <- lapply(workinginfo,sqerror,lhs)
  residses <- lapply(workinginfo,residse)
  
  if (withf1) {
    useable <- rootsokay & (f0est>0)
  } else {
    useable <- rootsokay & (f0est>0) & (f1est>0)
  }
  
  
  workinginfo$useful <- cbind(f0est,rootsokay,sqerrors,residses,useable)
  return(workinginfo)
}
