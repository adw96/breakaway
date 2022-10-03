#' Species richness estimation with Kemp-type models
#' 
#' This function implements the species richness estimation procedure outlined
#' in Willis & Bunge (2015). The diversity estimate, standard error, estimated
#' model coefficients, model details and plot of the fitted model are returned.
#' 
#' 
#' @param input_data An input type that can be processed by \code{convert()}
#' @param cutoff The maximum frequency count to use for model fitting
#' @param output Deprecated; only for backwards compatibility
#' @param answers Deprecated; only for backwards compatibility
#' @param plot Deprecated; only for backwards compatibility
#' @param print Deprecated; only for backwards compatibility
#' @param ... Additional arguments will be ignored; this is for backward compatibility
#' @return An object of class \code{alpha_estimate} \item{code}{ A category representing algorithm behaviour.
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
#' 
#' @author Amy Willis
#' 
#' @import magrittr
#' @import ggplot2
#' @import phyloseq
#' 
#' @seealso \code{\link{breakaway}}; \code{\link{breakaway_nof1}}; \code{\link{apples}}
#' @references Willis, A. and Bunge, J. (2015). Estimating diversity via
#' frequency ratios. \emph{Biometrics}, \bold{71}(4), 1042--1049.
#' 
#' Rocchetti, I., Bunge, J. and Bohning, D. (2011). Population size estimation
#' based upon ratios of recapture probabilities. \emph{Annals of Applied
#' Statistics}, \bold{5}.
#' @keywords diversity microbial models nonlinear
#' @examples
#' 
#' kemp(apples)
#' kemp(apples, plot = FALSE, output = FALSE, answers = TRUE)
#' 
#' @export 
kemp <- function(input_data,
                 cutoff = NA, 
                 output = NULL, plot = NULL, 
                 answers = NULL, print = NULL, ...) {
  UseMethod("kemp", input_data)
}

#' @export
kemp.phyloseq <- function(input_data,
                          cutoff = NA, ...) {
  
  physeq_wrap(fn = kemp, 
              physeq = input_data,
              cutoff = cutoff)
}

#' @export
kemp.matrix <- function(input_data,
                        cutoff = NA, ...) {
  
  input_data %>%
    as.data.frame %>%
    kemp(cutoff = cutoff)
  
}

#' @export
kemp.data.frame <- function(input_data,
                            cutoff = NA, ...) {
  
  ## if a frequency count matrix...
  if (dim(input_data)[2] == 2 & 
      all(input_data[,1] %>% sort == input_data[,1])) {
    kemp.default(input_data, cutoff = cutoff)
  } else {
    
    warning("Assuming taxa are rows")
    
    input_data %>% 
      apply(2, kemp, cutoff = cutoff) %>%
      alpha_estimates
    
  }
}

#' @export
kemp.default <- function(input_data,
                         cutoff = NA, ...)  {
  
  my_data <- convert(input_data)
  
  if (my_data[1, 1] != 1 || my_data[1, 2]==0) {
    kemp_alpha_estimate <- alpha_estimate(estimand = "richness",
                                          estimate = NA,
                                          error = NA,
                                          model = "Kemp",
                                          name = "kemp",
                                          frequentist = TRUE, 
                                          parametric = TRUE,
                                          reasonable = FALSE,
                                          warnings = "no singleton count")
  } else {
    n <- sum(my_data$frequency)
    f1 <- my_data[1, 2]
    
    cutoff <- cutoff_wrap(my_data, requested = cutoff) 
    
    if (cutoff < 6) { ## check for unusual data structures
      
      kemp_alpha_estimate <- alpha_estimate(estimand = "richness",
                                            estimate = NA,
                                            error = NA,
                                            model = "Kemp",
                                            name = "kemp",
                                            frequentist = TRUE, 
                                            parametric = TRUE,
                                            reasonable = FALSE,
                                            warnings = "insufficient contiguous frequencies")
      
    } else {
      
      ## set up structures
      my_data <- my_data[1:cutoff, ]
      ratios <- my_data[-1, 2]/my_data[-cutoff, 2]
      xs <- 1:(cutoff-1)
      ys <- (xs+1)*ratios
      xbar <- mean(xs)
      
      lhs <- list("x" = xs-xbar, "y" = ratios)
      
      stopifnot(length(lhs$x) == length(lhs$y))
      
      stopifnot(length(lhs$x) == length(xs))
      
      weights_inv <- 1/xs
      run <- minibreak_all(lhs, xs, ys, my_data, weights_inv, withf1 = TRUE)
      result <- list()
      choice <- list()
      
      ### If no models converged...
      if (sum(as.numeric(run$useful[,5])) == 0) {
        
        kemp_alpha_estimate <- alpha_estimate(estimand = "richness",
                                              estimate = NA,
                                              error = NA,
                                              model = "Kemp",
                                              name = "kemp",
                                              frequentist = TRUE, 
                                              parametric = TRUE,
                                              reasonable = FALSE,
                                              warnings = "no kemp models converged",
                                              other = list(outcome = 0,
                                                           code = 1))
        
      } else { 
        ## Otherwise, YAY! Something worked for 1/x weighting
        choice$outcome <- 1
        choice$model <- rownames(run$useful)[min(which(run$useful[,5]==1))] #pick smallest
        choice$full <-  run[[noquote(choice$model)]]
        
        oldest <- run$useful[run$useful[,5]==1,1][[1]]
        est <- 0
        its <- 0
        while ( choice$outcome & abs(oldest-est) > 1 & its < 30) {
          
          oldest <- est
          C <- round(n+oldest,0)
          unscaledprobs <- c(1,cumprod(fitted(choice$full)))
          p <- unscaledprobs/sum(unscaledprobs)
          as <- p[-1]^2/p[-cutoff]^3 * (1-exp(-C*p[-cutoff]))^3/(1-exp(-C*p[-1]))^2 * (1-C*p[-cutoff]/(exp(C*p[-cutoff])-1))
          bs <- p[-1]/p[-cutoff]^2 * (1-exp(-C*p[-cutoff]))^2/(1-exp(-C*p[-1])) * (1-C*p[-1]/(exp(C*p[-1])-1))
          ratiovars <- (as + bs)/C
          
          # if (its == 0) {
          #   weights1 <- 1/ratiovars
          # }
          
          run <- try ( minibreak_all(lhs, xs, ys, my_data,
                                     1/ratiovars, withf1 = TRUE),
                       silent = 1)
          
          # stopifnot(any(ratiovars < 0))
          # run <- minibreak_all(lhs, xs, ys, my_data,
          #                            1/ratiovars, withf1 = TRUE)
          
          if (inherits(run, "try-error")) {
            ratiovars <- (p[-1]^2/p[-cutoff]^3 + p[-1]/p[-cutoff]^2)/C
            # stopifnot(all(ratiovars > 0))
            run <- try ( minibreak_all(lhs, xs, ys, my_data, 1/ratiovars, withf1 = TRUE), silent = 1)
            if (inherits(run, "try-error")) {
              
              kemp_alpha_estimate <- alpha_estimate(estimand = "richness",
                                                    estimate = NA,
                                                    error = NA,
                                                    model = "Kemp",
                                                    name = "kemp",
                                                    frequentist = TRUE, 
                                                    parametric = TRUE,
                                                    reasonable = FALSE,
                                                    warnings = "no kemp models converged (numerical errors)",
                                                    other = list(outcome = 0,
                                                                 code = 1))
              
            }
          }
          
          choice <- list()
          if (inherits(run, "try-error")) {
            choice$outcome <- 0
          } else if (sum(as.numeric(run$useful[,5]))==0 | 
                     any(is.infinite(ratiovars)) |
                     ratiovars[1] > 1e20 |
                     any(ratiovars < 0)) {
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
        if (!choice$outcome) {
          # if(output) cat("We used 1/x weighting. \n")
          run <- minibreak_all(lhs, xs, ys, my_data, weights_inv, withf1 = TRUE)
          choice$outcome <- 1
          choice$model <- rownames(run$useful)[min(which(run$useful[,5]==1))]
          choice$full <-  run[[noquote(choice$model)]]
          result$code <- 3
          
        }
        
        if(choice$model=="model_1_0") {
          b0_hat <- coef(choice$full)[1]
          b0_var <- vcov(choice$full)[1,1]
        } else {
          effective_coef <- c(coef(choice$full), rep(0,9-length(coef(choice$full))))
          b <- effective_coef[1]-effective_coef[2]*xbar+effective_coef[4]*xbar^2-effective_coef[6]*xbar^3+effective_coef[8]*xbar^4
          a <- 1-effective_coef[3]*xbar+effective_coef[5]*xbar^2-effective_coef[7]*xbar^3+effective_coef[9]*xbar^4
          
          nabla <- c(1/a, -xbar/a, b*xbar/a^2, xbar^2/a, -b*xbar^2/a^2, -xbar^3/a, b*xbar^3/a^2, xbar^4/a, -b*xbar^4/a^2)
          nabla <- nabla[1:length(coef(choice$full))]
          
          b0_hat <- b/a
          b0_var <- t(nabla) %*% vcov(choice$full) %*% nabla
        }
        
        f0 <- run$useful[rownames(run$useful) == choice$model,1][[1]]
        f0_var <- f1*b0_hat^-2*(1 - f1/n + f1*b0_hat^-2 * b0_var) #1st order d.m.
        
        diversity <- unname(f0 + n)
        diversity_se <- unname(c(sqrt(n*f0/diversity + f0_var)))
        
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
        
        d <- exp(1.96*sqrt(log(1+diversity_se^2/f0)))
        
        yhats <- fitted(choice$full)
        
        plot_data <- rbind(data.frame("x" = xs, 
                                      "y" = lhs$y, 
                                      "type" = "Observed"),
                           data.frame("x" = xs, 
                                      "y" = yhats, 
                                      "type" = "Fitted"),
                           data.frame("x" = 0, 
                                      "y" = b0_hat, 
                                      "type" = "Prediction"))
        my_plot <- ggplot(plot_data, 
                          aes_string(x = "x", 
                                     y = "y",
                                     col = "type", 
                                     pch = "type")) +
          geom_point() +
          labs(x = "x", y = "f(x+1)/f(x)", title = "Plot of ratios and fitted values: Kemp models") +
          theme_bw()
        
        kemp_alpha_estimate <- alpha_estimate(estimate = diversity,
                                              error = diversity_se,
                                              model = "Kemp",
                                              name = "kemp",
                                              estimand = "richness",
                                              interval = c(n + f0/d, n + f0*d),
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
                                              full = choice$full,
                                              f0_hat = f0)
        
        
      }
    }
  }
  kemp_alpha_estimate
  
}

rootcheck <- function(model, lhs, nof1=FALSE) {
  # TODO rewrite to be clearer
  if(length(coef(model)) == 2)  {
    root <- -coef(model)[1]/coef(model)[2]
  } else {
    model_coefs_names <- model %>% coef %>% names %>% substring(1, 1)
    a_coefs <- model_coefs_names == "a"
    numerator_coefs <- c(1, coef(model)[a_coefs])
    numerator_roots <- numerator_coefs %>% polyroot
    root <- numerator_roots[(numerator_roots %>% Im %>% abs) < 10^-5] %>% Re
  }
  if (length(root) == 0 | 
      (sum((root > ifelse(nof1, min(lhs$x - 2), min(lhs$x - 1))) & (root < max(lhs$x))) == 0 )) { 
    # uses sum of pointwise booleans to ensure any roots are caught
    # edit 2019/11/04: want to make sure root is not between f1 = 1 and f1 = 0
    # This was issue with issue68
    outcome <- 1
  } else {
    outcome <- 0
  }
  return(outcome)
}

sqerror <- function(model, lhs) {
  fits <- fitted(model)
  
  stopifnot( length(lhs$y) == length(fits))
  
  return(sum(((lhs$y-fits)^2)/fits))
}

residse <- function(model) {
  return(summary(model)$sigma)
}

minibreak_all <- function(lhs, xs, ys, input_data, myweights, withf1 = NULL) {
  
  stopifnot( length(lhs[[1]]) == length(lhs[[2]]))
  
  stopifnot( length(lhs$x) == length(xs))
  stopifnot( length(lhs$x) == length(lhs$y))
  
  if (is.null(withf1)) stop("Empty argument `withf1`")
  
  structure_1_0 <- function(x, beta0, beta1) (beta0+beta1*x)/(1+x)
  structure_1_1 <- function(x, beta0, beta1,alpha1) (beta0+beta1*x)/(1+alpha1*x)
  structure_2_1 <- function(x, beta0, beta1,alpha1, beta2) (beta0+beta1*x+beta2*x^2)/(1+alpha1*x)
  structure_2_2 <- function(x, beta0, beta1,alpha1, beta2,alpha2) (beta0+beta1*x+beta2*x^2)/(1+alpha1*x+alpha2*x^2)
  structure_3_2 <- function(x, beta0, beta1,alpha1, beta2,alpha2, beta3) (beta0+beta1*x+beta2*x^2+beta3*x^3)/(1+alpha1*x+alpha2*x^2)
  structure_3_3 <- function(x, beta0, beta1,alpha1, beta2,alpha2, beta3,alpha3) (beta0+beta1*x+beta2*x^2+beta3*x^3)/(1+alpha1*x+alpha2*x^2+alpha3*x^3)
  structure_4_3 <- function(x, beta0, beta1,alpha1, beta2,alpha2, beta3,alpha3, beta4) (beta0+beta1*x+beta2*x^2+beta3*x^3+beta4*x^4)/(1+alpha1*x+alpha2*x^2+alpha3*x^3)
  structure_4_4 <- function(x, beta0, beta1,alpha1, beta2,alpha2, beta3,alpha3, beta4,alpha4) (beta0+beta1*x+beta2*x^2+beta3*x^3+beta4*x^4)/(1+alpha1*x+alpha2*x^2+alpha3*x^3+alpha4*x^4)
  
  all_ests <- list()
  
  model_1_0 <-  lm(ys~xs, weights=myweights)
  start_1_0 <- list(beta0 = model_1_0$coef[1], beta1 = model_1_0$coef[2])
  
  lhs2 <- lhs; lhs2$x <- xs # otherwise singularity
  stopifnot( length(lhs2[[1]]) == length(lhs2[[2]]))
  
  result_1_0 <- try( nls( y ~ structure_1_0(x, beta0, beta1),
                          data = lhs2, start_1_0, weights=myweights),
                     silent=T )
  if (inherits(result_1_0, "try-error")) {
    coef_1_0 <- c(start_1_0$beta0, start_1_0$beta1)
  } else {
    coef_1_0 <- coef(result_1_0)
  }
  
  beta0_tilde <- coef_1_0[1]
  beta1_tilde <- coef_1_0[2]
  
  xbar <- ifelse(withf1, mean(xs), mean(c(1, xs)))
  
  start_1_1 <- list(beta0 = (beta0_tilde+beta1_tilde*xbar)/(1+xbar), beta1 = beta1_tilde/(1+xbar), alpha1 = 1/(1+xbar))
  result_1_1 <- try( nls( y ~ structure_1_1(x, beta0, beta1,alpha1), data = lhs, start_1_1, weights=myweights), silent=T )
  
  if (inherits(result_1_1, "try-error")) {
    temp_lm <- lm(y ~ x + I(x^2), weights=myweights, data = lhs)$coef
    start_2_1 <- list(beta0 = temp_lm[1], beta1 = temp_lm[2], alpha1 = 0, beta2 = temp_lm[3])
  } else {
    start_2_1 <- as.list(c(summary(result_1_1)$coef[,1],"beta2"=0))
  }
  
  result_2_1 <- try(nls( y ~ structure_2_1(x, beta0, beta1,alpha1, beta2), data = lhs, start_2_1, weights=myweights), silent=T)
  
  if (inherits(result_2_1, "try-error")) {
    temp_lm <- lm(y ~ x + I(x^2), weights=myweights, data = lhs)$coef
    start_2_2 <- list(beta0 = temp_lm[1], temp_lm[2], alpha1 = 0, beta2 = temp_lm[3], alpha2 = 0)
  } else start_2_2 <- as.list(c(summary(result_2_1)$coef[,1],"alpha2"=0))
  result_2_2 <- try(nls( y ~ structure_2_2(x, beta0, beta1,alpha1, beta2,alpha2), data = lhs, start_2_2, weights=myweights), silent=T)
  
  if (inherits(result_2_2, "try-error")) {
    temp_lm <- lm(y ~ x + I(x^2)+I(x^3), weights=myweights, data = lhs)$coef
    start_3_2 <- list(beta0 = temp_lm[1], beta1 = temp_lm[2], alpha1 = 0, beta2 = temp_lm[3], alpha2 = 0, beta3 = temp_lm[4])
  } else start_3_2 <- as.list(c(summary(result_2_2)$coef[,1],"beta3"=0))
  result_3_2 <- try( nls( y ~ structure_3_2(x, beta0, beta1,alpha1, beta2,alpha2, beta3), data = lhs, start_3_2, weights=myweights), silent=T)
  
  if (inherits(result_3_2, "try-error")) {
    temp_lm <- lm(y ~ x + I(x^2)+I(x^3), weights=myweights, data = lhs)$coef
    start_3_3 <- list(beta0 = temp_lm[1], beta1 = temp_lm[2], alpha1 = 0, beta2 = temp_lm[3], alpha2 = 0, beta3 = temp_lm[4], alpha3 = 0)
  } else start_3_3 <- as.list(c(summary(result_3_2)$coef[,1],"alpha3"=0))
  result_3_3 <- try( nls( y ~ structure_3_3(x, beta0, beta1,alpha1, beta2,alpha2, beta3,alpha3), data = lhs, start_3_3, weights=myweights), silent=T)
  
  if (inherits(result_3_3, "try-error")) {
    temp_lm <- lm(y ~ x + I(x^2)+I(x^3)+I(x^4), weights=myweights, data = lhs)$coef
    start_4_3 <- list(beta0 = temp_lm[1], beta1 = temp_lm[2], alpha1 = 0, beta2 = temp_lm[3], alpha2 = 0, beta3 = temp_lm[4], alpha3 = 0, beta4 = temp_lm[5])
  } else start_4_3 <- as.list(c(summary(result_3_3)$coef[,1],"beta4"=0))
  result_4_3 <- try( nls( y ~ structure_4_3(x, beta0, beta1,alpha1, beta2,alpha2, beta3,alpha3, beta4), data = lhs, start_4_3, weights=myweights), silent=T)
  
  if (inherits(result_4_3, "try-error")) {
    temp_lm <- lm(y ~ x + I(x^2)+I(x^3)+I(x^4), weights=myweights, data = lhs)$coef
    start_4_4 <- list(beta0 = temp_lm[1], beta1 = temp_lm[2], alpha1 = 0, beta2 = temp_lm[3], alpha2 = 0,
                      beta3 = temp_lm[4], alpha3 = 0, beta4 = temp_lm[5], alpha4 = 0)
  } else start_4_4 <- as.list(c(summary(result_4_3)$coef[,1],"alpha4"=0))
  result_4_4 <- try( nls( y ~ structure_4_4(x, beta0, beta1,alpha1, beta2,alpha2, beta3,alpha3, beta4,alpha4), data = lhs, start_4_4, weights=myweights), silent=T)
  
  info <- list("model_1_0"=result_1_0,"model_1_1"=result_1_1,
               "model_2_1"=result_2_1,"model_2_2"=result_2_2,"model_3_2"=result_3_2,
               "model_3_3"=result_3_3,"model_4_3"=result_4_3,"model_4_4"=result_4_4)
  
  converged <- lapply(info, class) == "nls"
  
  workinginfo <- info[converged]
  
  b0est <- lapply(workinginfo, predict, list(x = -xbar))
  if (length(b0est$model_1_0) != 0) {
    b0est$model_1_0 <- predict(workinginfo$model_1_0, list(x=0))
  }
  
  if (withf1) {
    
    f0est <- lapply(b0est,function(x) input_data[1,2]/x)
    
  } else {
    
    firstratioest <- lapply(workinginfo,predict, list(x=-xbar+1))
    if (length(firstratioest$model_1_0)!=0) firstratioest$model_1_0 <- predict(workinginfo$model_1_0, list(x=1))
    
    f1est <- lapply(firstratioest,function(x) input_data[1,2]/x)
    
    f0est <- as.numeric(f1est)/as.numeric(b0est)
    
  }
  
  rootsokay <- as.logical(lapply(workinginfo, 
                                 rootcheck, 
                                 lhs, 
                                 nof1 =! withf1))
  
  ## issue: model 1_0 is different
  
  sqerrors <- lapply(workinginfo, sqerror, lhs = lhs)
  
  residses <- lapply(workinginfo, residse)
  
  if (withf1) {
    useable <- rootsokay & (f0est>0)
  } else {
    useable <- rootsokay & (f0est>0) & (f1est>0)
  }
  
  workinginfo$useful <- cbind(f0est, rootsokay, sqerrors, residses, useable)
  return(workinginfo)
}
