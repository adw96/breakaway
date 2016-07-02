breakaway_nof1 <- function(data, print=TRUE, plot=TRUE, answers=FALSE, force=FALSE) {

  if( !(is.matrix(data) || is.data.frame(data))) {
    filename <- data
    ext <- substr(filename, nchar(filename)-2, nchar(filename))
    if (ext == "csv") {
      data <- read.table(file=filename, header=0,sep=",")
      if( data[1,1] !=1) data <- read.table(filename, header=1,sep=",")
    } else if (ext == "txt") {
      data <- read.table(file=filename, header=0)
    } else cat("Please input your data as a txt or csv file,
               or as an R dataframe or matrix.")
  }
  if ( is.factor(data[,1]) ) {
    fs <- as.numeric(as.character(data[,1]))
    data <- cbind(fs,data[,2])
    data <- data[data[,1]!=0,]
  }

  if(length(data)>1) {
    if (data[1,1]!=2 || data[1,2]==0) {
      if(print) cat("breakaway_nof1 is for when you have no singleton count.\nYou need a leading doubleton count!\n")
    } else {
      data <- data[!(data[,2]==0 | is.na(data[,2])),]
      orig_data <- data
      n <- sum(orig_data[,2])
      f2 <- data[1,2]

      cutoff <- ifelse(is.na(which(data[-1,1]-data[-length(data[,1]),1]>1)[1]),length(data[,1]),which(data[-1,1]-data[-length(data[,1]),1]>1)[1])
      data <- data[1:cutoff,]
      ys <- (data[1:(cutoff-1),1]+1)*data[2:cutoff,2]/data[1:(cutoff-1),2]
      xs <- 2:cutoff
      xbar <- mean(c(1,xs))
      lhs <- list("x"=xs-xbar,"y"=data[2:cutoff,2]/data[1:(cutoff-1),2])

      if ( cutoff < 5) { ### check for unusual data structures
        if(print) cat("You don't have enough contiguous frequencies.\n breakaway needs at least 6!\n")
      } else if ((force==FALSE) && ( (data[cutoff,2]/data[cutoff-1,2])>10 )) {
        cat("\tIt looks like you've concatenated some of your data!\n Please truncate and try again.\n")
      } else {
        weights_inv <- 1/xs
        run <- .minibreak_nof1(lhs,xs,ys,data,weights_inv)
        result <- list()
        choice <- list()

        if (sum(as.numeric(run$useful[,5]))==0) {
          choice$outcome <- 0
          if(print) cat("No breakaway models converged.")
          weights_trans <- (1/data[-1,2]+1/data[-cutoff,2])^-1
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

            run <- try ( .minibreak_nof1(lhs,xs,ys,data,1/ratiovars), silent = 1)

            if ( class(run) == "try-error") {
              ratiovars <- (p[-1]^2/p[-cutoff]^3 + p[-1]/p[-cutoff]^2)/C
              run <- try ( .minibreak_nof1(lhs,xs,ys,data,1/ratiovars), silent = 1)
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
            if(print) cat("Iterative reweighting didn't produce any outcomes after the first iteration, so we use 1/x\n")
            run <- .minibreak_nof1(lhs,xs,ys,data,weights_inv)
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
}

.minibreak_nof1 <- function(lhs, xs, ys, data, myweights) {
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
  xbar <- mean(c(1,xs))
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

  firstratioest <- lapply(workinginfo,predict,list(x=-xbar+1))
  if (length(firstratioest$model_1_0)!=0) firstratioest$model_1_0 <- predict(workinginfo$model_1_0,list(x=1))

  f1est <- lapply(firstratioest,function(x) data[1,2]/x)

  f0est <- as.numeric(f1est)/as.numeric(b0est)

  rootsokay <- as.logical(lapply(workinginfo,.rootcheck,lhs,nof1=TRUE))

  sqerrors <- lapply(workinginfo,.sqerror,lhs)
  residses <- lapply(workinginfo,.residse)

  useable <- rootsokay & (f0est>0) & (f1est>0)

  workinginfo$useful <- cbind(f0est,rootsokay,sqerrors,residses,useable)
  return(workinginfo)
}
