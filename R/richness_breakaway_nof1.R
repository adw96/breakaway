

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
#' @param input_data An input type that can be processed by \code{convert()}
#' @param output Deprecated; only for backwards compatibility
#' @param answers Deprecated; only for backwards compatibility
#' @param plot Deprecated; only for backwards compatibility
#' @param print Deprecated; only for backwards compatibility
#' @param force Deprecated; only for backwards compatibility
#' @return An object of class \code{alpha_estimate}  \item{code}{ A category representing algorithm behaviour.
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
#' diversity but spurious singletons. \emph{arXiv.}
#' 
#' Willis, A. and Bunge, J. (2015). Estimating diversity via frequency ratios.
#' \emph{Biometrics.}
#' @keywords diversity error microbial models nonlinear
#' @examples
#' 
#' breakaway_nof1(apples[-1, ])
#' breakaway_nof1(apples[-1, ], plot = FALSE, output = FALSE, answers = TRUE)
#' 
#' @export breakaway_nof1
breakaway_nof1 <- function(input_data, 
                           output = NULL, plot = NULL, 
                           answers = NULL, print = NULL, 
                           force = NULL) {
  
  my_data <- convert(input_data)
  
  if (my_data[1, 1] != 2 || my_data[1, 2]==0) {
    kemp_alpha_estimate <- alpha_estimate(estimand = "richness",
                                          estimate = NA,
                                          error = NA,
                                          model = "Kemp",
                                          name = "breakaway_nof1",
                                          frequentist = TRUE, 
                                          parametric = TRUE,
                                          reasonable = FALSE,
                                          warnings = "no doubleton count")
  } else {
    
    n <- sum(my_data$frequency)
    f2 <- input_data[1,2]
    
    cutoff <- ifelse(is.na(which(input_data[-1,1]-input_data[-length(input_data[,1]),1]>1)[1]),length(input_data[,1]),which(input_data[-1,1]-input_data[-length(input_data[,1]),1]>1)[1])
    input_data <- input_data[1:cutoff,]
    ys <- (input_data[1:(cutoff-1),1]+1)*input_data[2:cutoff,2]/input_data[1:(cutoff-1),2]
    xs <- 2:cutoff
    xbar <- mean(c(1,xs))
    lhs <- list("x"=xs-xbar,"y"=input_data[2:cutoff,2]/input_data[1:(cutoff-1),2])
    
    if (cutoff < 6) { ## check for unusual data structures
      
      kemp_alpha_estimate <- alpha_estimate(estimand = "richness",
                                            estimate = NA,
                                            error = NA,
                                            model = "Kemp",
                                            name = "breakaway_nof1",
                                            frequentist = TRUE, 
                                            parametric = TRUE,
                                            reasonable = FALSE,
                                            warnings = "insufficient contiguous frequencies")
      
    } else {
      weights_inv <- 1/xs
      run <- minibreak_all(lhs,xs,ys,input_data,weights_inv, withf1 = FALSE)
      result <- list()
      choice <- list()
      
      if (sum(as.numeric(run$useful[,5]))==0) {
        
        
        kemp_alpha_estimate <- alpha_estimate(estimand = "richness",
                                              estimate = NA,
                                              error = NA,
                                              model = "Kemp",
                                              name = "breakaway_nof1",
                                              frequentist = TRUE, 
                                              parametric = TRUE,
                                              reasonable = FALSE,
                                              warnings = "no kemp models converged",
                                              other = list(outcome = 0,
                                                           code = 1))
        # 
        # choice$outcome <- 0
        # if(output) cat("No breakaway models converged.")
        # weights_trans <- (1/input_data[-1,2]+1/input_data[-cutoff,2])^-1
        # lm1 <- lm(log(ys)~xs,weights=weights_trans)
        # b0_hat <- summary(lm1)$coef[1,1]
        # b0_se <- summary(lm1)$coef[1,2]
        # b1_hat <- summary(lm1)$coef[2,1]
        # b1_se <- summary(lm1)$coef[2,2]
        # 
        # f1_pred <- 2*f2*exp(-(b0_hat+b1_hat))
        # f0_pred <- f1_pred*exp(-b0_hat)
        # diversity <- f0_pred + f1_pred + n
        # 
        # covmatrix <- matrix(c(f2*(1-f2/n),rep(0,8)),nrow=3)
        # covmatrix[2:3,2:3] <- vcov(lm1)
        # 
        # derivs <- c(exp(-(b0_hat+b1_hat))*(1+exp(-b1_hat)),
        #             f2*exp(-b1_hat)*(1+exp(-b1_hat))*exp(-b0_hat),
        #             f2*exp(-b0_hat)*(-exp(-b1_hat)-2*exp(-b1_hat)))
        # 
        # f0plusf1_var <- 4*t(derivs)%*%covmatrix%*%(derivs)
        # 
        # diversity_se <- sqrt(f0plusf1_var+n*(f0_pred+f1_pred)/(n+f0_pred+f1_pred))
        # 
        # if(output) {
        #   cat("################## breakaway ##################\n")
        #   cat("\tThe best estimate of total diversity is", round(diversity),
        #       "\n \t with std error",round(diversity_se),"\n")
        #   cat("\tThe model employed was the WLRM\n")
        # }
        # if(answers) {
        #   result$name <- "WLRM"
        #   result$para <- summary(lm1)$coef[,1:2]
        #   result$est <- diversity
        #   result$seest <- as.vector(diversity_se)
        #   result$full <- lm1
        #   d <- exp(1.96*sqrt(log(1+result$seest^2/f0_pred)))
        #   result$ci <- c(n+f0_pred/d,n+f0_pred*d)
        #   return(result)
        # }
        # 
        # result$code <- 1
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
          
          run <- try ( minibreak_all(lhs, xs, ys, input_data, 1/ratiovars, withf1=FALSE), silent = 1)
          
          if ( class(run) == "try-error") {
            ratiovars <- (p[-1]^2/p[-cutoff]^3 + p[-1]/p[-cutoff]^2)/C
            run <- try ( minibreak_all(lhs,xs,ys,input_data,1/ratiovars, withf1=FALSE), silent = 1)
            if ( class(run) == "try-error") {
              # if(output) {print("Numerical errors result in non-convergence") }
              kemp_alpha_estimate <- alpha_estimate(estimand = "richness",
                                                    estimate = NA,
                                                    error = NA,
                                                    model = "Kemp",
                                                    name = "breakaway_nof1",
                                                    frequentist = TRUE, 
                                                    parametric = TRUE,
                                                    reasonable = FALSE,
                                                    warnings = "no kemp models converged (numerical errors)",
                                                    other = list(outcome = 0,
                                                                 code = 1))
              
            }
          }
          
          choice <- list()
          if (class(run) == "try-error" | sum(as.numeric(run$useful[,5]))==0) {
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
          # if(output) cat("Iterative reweighting didn't produce any outcomes after the first iteration, so we use 1/x\n")
          run <- minibreak_all(lhs,xs,ys,input_data,weights_inv, withf1 = FALSE)
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
        diversity <- unname(f0_pred + f1_pred + n)
        
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
        f0plusf1_var <- derivs %*% covmatrix %*% t(derivs)
        
        var_est <- sum(f0plusf1_var) + n*(f0_pred+f1_pred)/diversity - 2*f0_pred*n/diversity - 2*f1_pred*n/diversity
        diversity_se <- unname(ifelse(var_est<0, 0, c(sqrt(var_est)))) ## Piece 6
        
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
        
        d <- exp(1.96*sqrt(log(1 + diversity_se^2/f0_pred)))
        
        # if(output) {
        #   cat("################## breakaway ##################\n")
        #   cat("\tThe best estimate of total diversity is", round(diversity),
        #       "\n \t with std error",round(diversity_se),"\n")
        #   cat("\tThe model employed was", choice$model,"\n")
        #   cat("\tThe function selected was\n\t ", the_function,"\n")
        #   print(parameter_table)
        #   cat("xbar\t\t\t",xbar)
        # }
        
        # if(plot)  {
        #   yhats <- fitted(choice$full)
        #   par(new=FALSE)
        #   plot(xs,lhs$y,xlim=c(0,max(xs)+1),ylim=c(min(lhs$y,yhats),max(lhs$y)*1.05),
        #        ylab="f(x+1)/f(x)",xlab="x",main="Plot of ratios and fitted values under model");
        #   points(xs,yhats,pch=18)
        #   points(c(0,1),c(b0_hat,c0_hat),col="red",pch=18)
        #   legend(0,max(lhs$y),c("Fitted values", "Prediction"),pch=c(18,18),col=c("black", "red"),cex=0.8,bty = "n")
        # }
        
        yhats <- fitted(choice$full)
        
        plot_data <- rbind(data.frame("x" = xs, "y" = lhs$y, 
                                      "type" = "Observed"),
                           data.frame("x" = xs, "y" = yhats, 
                                      "type" = "Fitted"),
                           data.frame("x" = c(0, 1), 
                                      "y" = c(b0_hat,c0_hat), 
                                      "type" = "Prediction"))
        
        my_plot <- ggplot(plot_data, 
                          aes(x = x, 
                              y = y,
                              col = type, pch = type)) +
          geom_point() +
          labs(x = "x", y = "f(x+1)/f(x)", title = "Plot of ratios and fitted values: Kemp models (no f1)") +
          theme_bw()
        
        
        # if(answers) {
        #   result$name <- choice$model
        #   result$para <- parameter_table
        #   result$est <- diversity
        #   result$seest <- as.vector(diversity_se)
        #   result$full <- choice$full
        #   result$ci <- c(n+f0_pred/d,n+f0_pred*d)
        #   return(result)
        # }
        
        kemp_alpha_estimate <- alpha_estimate(estimate = diversity,
                                              error = diversity_se,
                                              model = "Kemp",
                                              name = "breakaway+nof1",
                                              estimand = "richness",
                                              interval = c(n + f0_pred/d, n + f0_pred*d),
                                              interval_type = "Approximate: log-normal", 
                                              frequentist = TRUE, 
                                              parametric = TRUE,
                                              reasonable = TRUE,
                                              warnings = NULL,
                                              # legacy arguments
                                              para = parameter_table,
                                              plot = my_plot,
                                              other = list(xbar = xbar,
                                                           code = 1,
                                                           the_function = the_function,
                                                           name = choice$model),
                                              full = choice$full)
        
        
      }
    }
  }
  kemp_alpha_estimate
}
