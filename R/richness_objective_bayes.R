#' Objective Bayes species richness estimate with the Negative Binomial model
#'
#' @param data A frequency count table
#' @param output Print summary of fitted model. Default is TRUE.
#' @param plot Plot estimates? Default is TRUE.
#' @param answers Return dataframe containing model summaries? Default is FALSE.
#' @param tau Truncation cutoff for model fitting. (Counts above this cutoff
#' are not used in model fitting.) Default is 10.
#' @param burn.in Number of MCMC burn-in iterations. Default is 1000.
#' @param iterations Number of MCMC iterations after burn-in. Default is 5000.
#' @param Metropolis.stdev.N Standard deviation of proposal distribution
#' parameter N within MCMC. Default is 100.
#' @param Metropolis.start.T1 Starting value of parameter T1. Default is -0.8.
#' @param Metropolis.stdev.T1 Standard deviation of proposal distribution
#' parameter T1 within MCMC. Default is 0.01.
#' @param Metropolis.start.T2 Starting value of parameter T2. Default is 0.8.
#' @param Metropolis.stdev.T2 Standard deviation of proposal distribution
#' parameter T2 within MCMC. Default is 0.01.
#' @return If answers set to TRUE, a list containing posterior median and
#' 95% posterior credible intervals for richness, fitted values, and
#' MCMC diagnostics. If answers set to FALSE (default), no value is returned.
#' @importFrom stats acf
#' @importFrom graphics hist par plot
#'
objective_bayes_negbin <- function(data,
                                   output=TRUE,
                                   plot=TRUE,
                                   answers=FALSE,
                                   tau=10,
                                   burn.in=1000,
                                   iterations=5000,
                                   Metropolis.stdev.N=100,
                                   Metropolis.start.T1=-0.8,
                                   Metropolis.stdev.T1=0.01,
                                   Metropolis.start.T2=0.8,
                                   Metropolis.stdev.T2=0.01#,
                                   #bars=5 apparent an orphaned argument
                                   ) {

  data <- check_format(data)
  fullfreqdata  <- data

  if (tau > max(data[,1])) {
    tau <- max(data[,1])
  }

  # calculate summary statistics on full data
  w<-sum(fullfreqdata[,2])
  n<-sum(fullfreqdata[,1]*fullfreqdata[,2])

  # subset data up to tau
  freqdata<-fullfreqdata[1:tau,]

  # calculate summary statistics on data up to tau
  w.tau<-sum(freqdata[,2])
  n.tau<-sum(freqdata[,1]*freqdata[,2])
  # calculate NP estimate of n0
  NP.est.n0<-w.tau/(1-freqdata[1,2]/n.tau)-w.tau


  ### Step 3: calculate posterior

  ## initialization
  iterations<-iterations+burn.in
  N<-rep(0,iterations)
  T1T2<-matrix(rep(c(Metropolis.start.T1,Metropolis.start.T2),each=iterations),ncol=2)
  # to track acceptance rate of T1T2
  a1<-0
  # to track acceptance rate of N
  a2<-0
  # starting value based on nonparametric estimate of n0
  N[1]<-ceiling(NP.est.n0)+w.tau
  # storage for deviance replicates
  D.post<-rep(0,iterations)



  for (i in 2:iterations){

    # print every 500th iteration number
    if (i %in% seq(0,iterations-burn.in,by=500)) {message(paste("starting iteration ",i," of ",iterations,sep=""))}

    ## sample from p(T1T2|N,x)

    ## propose value T1T2 from a bivariate normal dist.; make sure T1T2.new > {-1,0}
    repeat {
      T1T2.new <- rmvnorm(1, c(T1T2[i-1,1],T1T2[i-1,2]),
                          matrix(c(Metropolis.stdev.T1,0,0,Metropolis.stdev.T2),nrow=2))
      if(T1T2.new[1]>(-1) & T1T2.new[2]>0)
        break
    }

    # calculate log of acceptance ratio
    logr1<-(-1)*log(T1T2.new[1]^2+2*T1T2.new[1]+2)-log(1+T1T2.new[2]^2)+n.tau*log(T1T2.new[2])-
      (N[i-1]*(1+T1T2.new[1])+n.tau)*log(1+T1T2.new[1]+T1T2.new[2])+
      N[i-1]*(1+T1T2.new[1])*log(1+T1T2.new[1])+
      sum(freqdata[,2]*lgamma(1+T1T2.new[1]+freqdata[,1]))-
      w.tau*lgamma(1+T1T2.new[1])+log(T1T2[i-1,1]^2+2*T1T2[i-1,1]+2)+
      log(1+T1T2[i-1,2]^2)-n.tau*log(T1T2[i-1,2])+
      (N[i-1]*(1+T1T2[i-1,1])+n.tau)*log(1+T1T2[i-1,1]+T1T2[i-1,2])-
      N[i-1]*(1+T1T2[i-1,1])*log(1+T1T2[i-1,1])-
      sum(freqdata[,2]*lgamma(1+T1T2[i-1,1]+freqdata[,1]))+
      w.tau*lgamma(1+T1T2[i-1,1])

    # calculate acceptance ratio
    r1<-exp(logr1)

    # accept or reject propsed value
    if (runif(1) < min(r1,1)) {
      T1T2[i,]<-T1T2.new
      a1<-a1+1
    } else {
      T1T2[i,]<-T1T2[i-1,]
    }

    ## sample from p(N|A,G,x)

    ## make sure N.new >=w.tau
    repeat {
      N.new<-rnbinom(1,mu=N[i-1],size=Metropolis.stdev.N)
      if(N.new>w.tau-1)
        break
    }

    ## calculate log(N.new!/(N.new-w.tau)!)
    N3.new<-rep(0,w.tau)
    for (j in 0:(w.tau-1)){
      N3.new[j+1]<-log(N.new-j)
    }
    N2.new<-sum(N3.new)

    ## calculate log(N[i-1]!/(N[i-1]-w.tau)!)
    N3<-rep(0,w.tau)
    for (j in 0:(w.tau-1)){
      N3[j+1]<-log(N[i-1]-j)
    }
    N2<-sum(N3)

    # calculate log of acceptance ratio
    logr2<-(-1/2)*log(N.new)+N2.new+N.new*(1+T1T2[i,1])*log((1+T1T2[i,1])/(1+T1T2[i,1]+T1T2[i,2]))+
      log(dnbinom(N[i-1],mu=N.new,size=Metropolis.stdev.N))+
      (1/2)*log(N[i-1])-N2-
      N[i-1]*(1+T1T2[i,1])*log((1+T1T2[i,1])/(1+T1T2[i,1]+T1T2[i,2]))-
      log(dnbinom(N.new,mu=N[i-1],size=Metropolis.stdev.N))

    # calculate acceptance ratio
    r2<-exp(logr2)

    # accept or reject propsed value
    if (runif(1)<min(r2,1)) {N[i]<-N.new ; a2<-a2+1} else {N[i]<-N[i-1]}

    ## calculate deviance from current sample

    # calculate log(N[i]!/(N[i]-w.tau)!)
    N3.curr<-rep(0,w.tau)
    for (j in 0:(w.tau-1)){
      N3.curr[j+1]<-log(N[i]-j)
    }
    N2.curr<-sum(N3.curr)

    # calculate deviance
    D.post[i]<-(-2)*(N2.curr+n.tau*log(T1T2[i,2])-(N[i]*(1+T1T2[i,1])+n.tau)*log(1+T1T2[i,1]+T1T2[i,2])+N[i]*(1+T1T2[i,1])*log(1+T1T2[i,1])+sum(freqdata[,2]*lgamma(1+T1T2[i,1]+freqdata[,1]))-w.tau*lgamma(1+T1T2[i,1])-sum(lfactorial(freqdata[,2]))-sum(freqdata[,2]*lgamma(freqdata[,1]+1)))

  }


  ### Step 4: model diagnostics

  ## 1) deviance at posterior mean
  mean.T1<-mean(T1T2[(burn.in+1):iterations,1])
  mean.T2<-mean(T1T2[(burn.in+1):iterations,2])
  mean.N<-mean(N[(burn.in+1):iterations])

  ## calculate log(mean.N!/(mean.N-w.tau)!)
  N3.mean<-rep(0,w.tau)
  for (j in 0:(w.tau-1)){
    N3.mean[j+1]<-log(mean.N-j)
  }
  N2.mean<-sum(N3.mean)

  loglik.post.mean<-N2.mean+n.tau*log(mean.T2)-(mean.N*(1+mean.T1)+n.tau)*log(1+mean.T1+mean.T2)+mean.N*(1+mean.T1)*log(1+mean.T1)+sum(freqdata[,2]*lgamma(1+mean.T1+freqdata[,1]))-w.tau*lgamma(1+mean.T1)-sum(lfactorial(freqdata[,2]))-sum(freqdata[,2]*lgamma(freqdata[,1]+1))
  D.mean<-(-2)*loglik.post.mean

  ## 2) posterior mean and median deviances
  mean.D<-mean(D.post[(burn.in+1):iterations])
  median.D<-quantile(D.post[(burn.in+1):iterations],probs=.5,names=F)

  ## 3) model complexity
  p.D<-mean.D-D.mean

  ## 4) Deviance information criterion
  DIC<-2*p.D+D.mean


  ### Step 5: fitted values based on medians of the marginal posteriors
  median.T1<-quantile(T1T2[(burn.in+1):iterations,1],probs=.5,names=F)
  median.T2<-quantile(T1T2[(burn.in+1):iterations,2],probs=.5,names=F)
  median.N<-quantile(N[(burn.in+1):iterations],probs=.5,names=F)

  fits<-rep(0,tau)
  for (k in 1:tau){
    fits[k]<-(median.N)*dnbinom(k,size=median.T1+1,prob=(median.T1+1)/(median.T1+median.T2+1))
  }
  fitted.values<-data.frame(cbind(j=seq(1,tau),fits,count=freqdata[,2]))


  ### Step 6: estimate thinning to reduce correlated posterior samples
  lags<-acf(N[(burn.in+1):iterations],type="correlation",main="Autocorr plot",ylab="ACF",xlab="Lag", plot=F)
  lag.thin<-suppressWarnings(min(which(lags$acf<0.1)))
  if (lag.thin==Inf) {lag.thin<-paste(">",length(lags$lag),sep="")
  }


  ### Step 7: results
  hist.points<-hist(N[(burn.in+1):iterations]+w-w.tau,breaks=seq(w,max(N)+w-w.tau+1)-0.5, plot = FALSE)

  results<-data.frame(w=w,
                      n=n,
                      NP.est.N=NP.est.n0+w,
                      tau=tau,
                      w.tau=w.tau,
                      n.tau=n.tau,
                      iterations=iterations-burn.in,
                      burn.in=burn.in,
                      acceptance.rate.T1T2=a1/iterations,
                      acceptance.rate.N=a2/iterations,
                      lag=lag.thin,
                      mode.N=hist.points$mids[which.max(hist.points$density)],
                      mean.N=mean(N[(burn.in+1):iterations])+w-w.tau,
                      median.N=quantile(N[(burn.in+1):iterations]+w-w.tau,probs=.5,names=F),
                      LCI.N=quantile(N[(burn.in+1):iterations]+w-w.tau,probs=.025,names=F),
                      UCI.N=quantile(N[(burn.in+1):iterations]+w-w.tau,probs=.975,names=F),
                      stddev.N=sd((N[(burn.in+1):iterations]+w-w.tau)),
                      mean.D=mean.D,
                      median.D=median.D,
                      DIC
  )

  final_results <- list()
  final_results$est <- results$median.N
  final_results$ci <- c("lower 95%"=results$LCI.N, "upper 95%"=results$UCI.N)
  final_results$mean <- results$mean.N
  final_results$semeanest <- results$stddev.N
  final_results$dic  <- DIC
  final_results$fits <- fitted.values
  final_results$diagnostics<-c("acceptance rate N"=results$acceptance.rate.N,
                               "acceptance rate T1T2"=results$acceptance.rate.T1T2,
                               "lag"=results$lag)

  if (output) {
    # output results and fitted values
    print(final_results)
  }

  if (plot) {
    oldpar <- par(no.readonly = TRUE)
    par(mfrow=c(1,2))

    ## posterior histogram
    hist(N[(burn.in+1):iterations]+w-w.tau,
         main="Posterior distribution",xlab="Total Number of Species",
         col='purple',freq=F,ylab="Density")

    # make trace plot
    plot((burn.in+1):iterations,N[(burn.in+1):iterations]+w-w.tau,type="l",xlab="Iteration Number",ylab="Total Number of Species", main="Trace plot")
    on.exit(par(oldpar))
  }

  if (answers) {
    return(final_results)
  }

}


#' Objective Bayes species richness estimate with the Poisson model
#'
#' @param data A frequency count table
#' @param output Print summary of fitted model. Default is TRUE.
#' @param plot Plot estimates? Default is TRUE.
#' @param answers Return dataframe containing model summaries? Default is FALSE.
#' @param tau Truncation cutoff for model fitting. (Counts above this cutoff
#' are not used in model fitting.) Default is 10.
#' @param burn.in Number of MCMC burn-in iterations. Default is 100.
#' @param iterations Number of MCMC iterations after burn-in. Default is 2500.
#' @param Metropolis.stdev.N Standard deviation of proposal distribution
#' parameter N within MCMC. Default is 75.
#' @param Metropolis.start.lambda Starting value of parameter lambda.
#' Default is 1
#' @param Metropolis.stdev.lambda Standard deviation of proposal distribution
#' parameter T1 within MCMC. Default is 0.3.
#' @return If answers set to TRUE, a dataframe containing posterior median and
#' 95% posterior credible intervals for richness, fitted values, and
#' MCMC diagnostics. If answers set to FALSE (default), no value is returned.
#' @importFrom graphics hist par plot
#'
#' @export
objective_bayes_poisson <- function(data,
                                    output=TRUE,
                                    plot=TRUE,
                                    answers=FALSE,
                                    tau=10, burn.in=100,
                                    iterations=2500,
                                    Metropolis.stdev.N=75,
                                    Metropolis.start.lambda=1,
                                    Metropolis.stdev.lambda=0.3#,
                                    #bars=5
                                    ) {

  data <- check_format(data)
  fullfreqdata  <- data

  if (tau > max(data[,1])) {
    tau <- max(data[,1])
  }

  # calculate summary statistics on full data
  w<-sum(fullfreqdata[,2])
  n<-sum(fullfreqdata[,1]*fullfreqdata[,2])

  # subset data up to tau
  freqdata<-fullfreqdata[1:tau,]

  # calculate summary statistics on data up to tau
  w.tau<-sum(freqdata[,2])
  n.tau<-sum(freqdata[,1]*freqdata[,2])
  # calculate NP estimate of n0
  NP.est.n0<-w.tau/(1-freqdata[1,2]/n.tau)-w.tau


  ### Step 3: calculate posterior

  ## initialization
  iterations<-iterations+burn.in
  N<-rep(0,iterations)
  L<-c(Metropolis.start.lambda,rep(1,iterations-1))
  # to track acceptance rate of lambda
  a1<-0
  # to track acceptance rate of N
  a2<-0
  # starting value based on nonparametric estimate of n0
  N[1]<-ceiling(NP.est.n0)+w.tau
  # storage for deviance replicates
  D.post<-rep(0,iterations)

  for (i in 2:iterations){

    # print every 500th iteration number
    if (i %in% seq(0,iterations-burn.in,by=500)) {message(paste("starting iteration ",i," of ",iterations,sep=""))}

    ## sample from p(lambda|x,C)

    # propose value for lambda
    L.new<-abs(rnorm(1,mean=L[i-1],sd=Metropolis.stdev.lambda))

    # calculate log of acceptance ratio
    logr1<-(n.tau-1/2)*log(L.new)-L.new*N[i-1]-(n.tau-1/2)*log(L[i-1])+L[i-1]*N[i-1]

    # calculate acceptance ratio
    r1<-exp(logr1)

    # accept or reject propsed value
    if (runif(1)<min(r1,1)) {L[i]<-L.new ; a1<-a1+1} else {L[i]<-L[i-1]}

    ## sample from p(N|lambda,x)

    ## make sure N.new >=w.tau
    repeat {
      N.new<-rnbinom(1,mu=N[i-1],size=Metropolis.stdev.N)
      if(N.new>w.tau-1)
        break
    }

    ## calculate log(N.new!/(N.new-w.tau)!)
    N3.new<-rep(0,w.tau)
    N3.new[1:w.tau]<-log(N.new-0:(w.tau-1))
    N2.new<-sum(N3.new)

    ## calculate log(N[i-1]!/(N[i-1]-w.tau)!)
    N3<-rep(0,w.tau)
    N3[1:w.tau]<-log(N[i-1]- 0:(w.tau-1))
    N2<-sum(N3)

    # calculate log of acceptance ratio
    logr2<-(N2.new-(1/2)*log(N.new)-N.new*L[i])-(N2-(1/2)*log(N[i-1])-N[i-1]*L[i])+(log(dnbinom(N[i-1],mu=N.new,size=Metropolis.stdev.N)))-(log(dnbinom(N.new,mu=N[i-1],size=Metropolis.stdev.N)))

    # calculate acceptance ratio
    r2<-exp(logr2)

    # accept or reject propsed value
    if (runif(1)<min(r2,1)) {N[i]<-N.new ; a2<-a2+1} else {N[i]<-N[i-1]}

    ## calculate deviance from current sample

    # calculate log(N[i]!/(N[i]-w.tau)!)
    N3.curr<-rep(0,w.tau)
    for (j in 0:(w.tau-1)){
      N3.curr[1:w.tau]<-log(N[i]-0:(w.tau-1))
    }
    N2.curr<-sum(N3.curr)

    # calculate deviance
    D.post[i]<-(-2)*(N2.curr-sum(lfactorial(freqdata[,2]))-L[i]*(N[i])-sum(freqdata[,2]*log(factorial(freqdata[,1])))+n.tau*log(L[i]))

  }


  ### Step 4: model diagnostics

  ## 1) deviance at posterior mean
  mean.L<-mean(L[(burn.in+1):iterations])
  mean.N<-mean(N[(burn.in+1):iterations])

  ## calculate log(mean.N!/(mean.N-w.tau)!)
  N3.mean<-rep(0,w.tau)
  N3.mean[1:w.tau]<-log(mean.N-0:(w.tau-1))
  N2.mean<-sum(N3.mean)

  loglik.post.mean<-N2.mean-sum(lfactorial(freqdata[,2]))-mean.L*mean.N+n.tau*log(mean.L)-sum(freqdata[,2]*lfactorial(freqdata[,1]))
  D.mean<-(-2)*loglik.post.mean

  ## 2) posterior mean and median deviances
  mean.D<-mean(D.post[(burn.in+1):iterations])
  median.D<-quantile(D.post[(burn.in+1):iterations],probs=.5,
                     names=F)

  ## 3) model complexity
  p.D<-mean.D-D.mean

  ## 4) Deviance information criterion
  DIC<-2*p.D+D.mean


  ### Step 5: fitted values based on medians of the marginal posteriors
  median.L<-quantile(L[(burn.in+1):iterations],probs=.5,names=F)
  median.N<-quantile(N[(burn.in+1):iterations],probs=.5,names=F)

  fits<-rep(0,tau)
  fits[1:tau]<-(median.N)*dpois(1:tau,median.L)
  fitted.values<-data.frame(cbind(j=seq(1,tau),fits,count=freqdata[,2]))


  ### Step 6: estimate thinning to reduce correlated posterior samples
  lags<-acf(N[(burn.in+1):iterations],type="correlation",main="Autocorr plot",ylab="ACF",xlab="Lag", plot=F)
  lag.thin<-suppressWarnings(min(which(lags$acf<0.1)))
  if (lag.thin==Inf) {lag.thin<-paste(">",length(lags$lag),sep="")
  }


  ### Step 7: results
  hist.points<-hist(N[(burn.in+1):iterations]+w-w.tau,breaks=seq(w,max(N)+w-w.tau+1)-0.5, plot = FALSE)

  results<-data.frame(w=w,
                      n=n,
                      NP.est.N=NP.est.n0+w,
                      tau=tau,
                      w.tau=w.tau,
                      n.tau=n.tau,
                      iterations=iterations-burn.in,
                      burn.in=burn.in,
                      acceptance.rate.lambda=a1/iterations,
                      acceptance.rate.N=a2/iterations,
                      lag=lag.thin,
                      mode.N=hist.points$mids[which.max(hist.points$density)],
                      mean.N=mean(N[(burn.in+1):iterations])+w-w.tau,
                      median.N=quantile(N[(burn.in+1):iterations]+w-w.tau,probs=.5,names=F),
                      LCI.N=quantile(N[(burn.in+1):iterations]+w-w.tau,probs=.025,names=F),
                      UCI.N=quantile(N[(burn.in+1):iterations]+w-w.tau,probs=.975,names=F),
                      stddev.N=sd((N[(burn.in+1):iterations]+w-w.tau)),
                      mean.D=mean.D,
                      median.D=median.D,
                      DIC)

  final_results <- list()
  final_results$est <- results$median.N
  final_results$ci <- c("lower 95%"=results$LCI.N, "upper 95%"=results$UCI.N)
  final_results$mean <- results$mean.N
  final_results$semeanest <- results$stddev.N
  final_results$dic  <- DIC
  final_results$fits <- fitted.values
  final_results$diagnostics<-c("acceptance rate N"=results$acceptance.rate.N,
                               "acceptance rate lambda"=results$acceptance.rate.lambda,
                               "lag"=results$lag)

  if (output) {
    # output results and fitted values
    print(final_results)
  }

  if (plot) {
    oldpar <- par(no.readonly = TRUE)
    par(mfrow=c(1,2))

    ## posterior histogram
    hist(N[(burn.in+1):iterations]+w-w.tau,
         main="Posterior distribution",xlab="Total Number of Species",
         col='purple',freq=F,ylab="Density")

    # make trace plot
    plot((burn.in+1):iterations,N[(burn.in+1):iterations]+w-w.tau,type="l",xlab="Iteration Number",ylab="Total Number of Species", main="Trace plot")
    on.exit(par(oldpar))
    }

  if (answers) {
    return(final_results)
  }

}

#' Objective Bayes species richness estimate with the mixed-geometric model
#'
#' @param data A frequency count table
#' @param output Print summary of fitted model. Default is TRUE.
#' @param plot Plot estimates? Default is TRUE.
#' @param answers Return dataframe containing model summaries? Default is FALSE.
#' @param tau Truncation cutoff for model fitting. (Counts above this cutoff
#' are not used in model fitting.) Default is 10.
#' @param burn.in Number of MCMC burn-in iterations. Default is 1000.
#' @param iterations Number of MCMC iterations after burn-in. Default is 5000.
#' @param Metropolis.stdev.N Standard deviation of proposal distribution
#' parameter N within MCMC. Default is 100.
#' @param Metropolis.start.T1 Starting value of parameter T1. Default is 1.
#' @param Metropolis.stdev.T1 Standard deviation of proposal distribution
#' parameter T1 within MCMC. Default is 2.
#' @param Metropolis.start.T2 Starting value of parameter T2. Default is 3.
#' @param Metropolis.stdev.T2 Standard deviation of proposal distribution
#' parameter T2 within MCMC. Default is 2.
#' @return If answers set to TRUE, a list containing posterior median and
#' 95% posterior credible intervals for richness, fitted values, and
#' MCMC diagnostics. If answers set to FALSE (default), no value is returned.
#' @importFrom graphics hist par plot
#'
#' @export
objective_bayes_mixedgeo <- function(data, output=TRUE, plot=TRUE, answers=FALSE,
                                     tau=10, burn.in=100, iterations=2500, Metropolis.stdev.N=100,
                                     Metropolis.start.T1=1, Metropolis.stdev.T1=2,
                                     Metropolis.start.T2=3, Metropolis.stdev.T2=2#, bars=3
                                     ) {

  data <- check_format(data)
  fullfreqdata  <- data

  if (tau > max(data[,1])) {
    tau <- max(data[,1])
  }

  # calculate summary statistics on full data
  w<-sum(fullfreqdata[,2])
  n<-sum(fullfreqdata[,1]*fullfreqdata[,2])

  # subset data up to tau
  freqdata<-fullfreqdata[1:tau,]

  # calculate summary statistics on data up to tau
  w.tau<-sum(freqdata[,2])
  n.tau<-sum(freqdata[,1]*freqdata[,2])
  # calculate NP estimate of n0
  NP.est.n0<-w.tau/(1-freqdata[1,2]/n.tau)-w.tau


  ### Step 3: calculate posterior

  ## initialization
  iterations<-iterations+burn.in
  A<-rep(0,iterations)
  T1<-rep(0,iterations)
  T2<-rep(0,iterations)
  N<-rep(0,iterations)
  # to track acceptance rate of N
  a1<-0
  # starting value, nonparametric estimate of n0
  N[1]<-ceiling(NP.est.n0)+w.tau
  # starting value, MLE of T1
  T1[1]<-Metropolis.start.T1
  # starting value, MLE of T2
  T2[1]<-Metropolis.start.T2
  A[1]<-0.5
  # storage for deviance replicates
  D.post<-rep(0,iterations)

  for (i in 2:iterations){

    # print every 500th iteration number
    if (i %in% seq(0,iterations-burn.in,by=500)) {message(paste("starting iteration ",i," of ",iterations,sep=""))}

    ## sample from p(Z|A,T1,T2,X,N)

    ## create a new vector of length N[i-1]
    Z<-rep(0,length=N[i-1])

    ## create a full data vector
    X<-c(rep(0,N[i-1]-w.tau),rep(freqdata[,1],times=freqdata[,2]))

    ## sample random bernoulli with appropriate success prob for each Z[k]; do not allow for Z all zeros or ones
    for (k in 1:N[i-1]){
      Z[k]<-rbinom(1,1,prob=A[i-1]*(1/(1+T1[i-1]))*(T1[i-1]/(1+T1[i-1]))^X[k]/((A[i-1]*(1/(1+T1[i-1]))*(T1[i-1]/(1+T1[i-1]))^X[k])+((1-A[i-1])
                                                                                                                                    *(1/(1+T2[i-1]))*(T2[i-1]/(1+T2[i-1]))^X[k])))
    }

    ## sample from p(A|Z,T1,T2,X,N)

    ## sample from beta dist
    A[i]<-rbeta(1,shape1=sum(Z)+1,shape2=N[i-1]-sum(Z)+1)

    ## sample from p(T1|A,Z,T2,X,N) and p(T2|A,Z,T1,X,N)

    repeat{
      ## sample T1/(1+T1) and T2/(1+T2) from beta dists
      T1.trans<-rbeta(1,shape1=sum(X*Z)+0.5,shape2=sum(Z)+0.5)
      T2.trans<-rbeta(1,shape1=sum(X)-sum(X*Z)+0.5,shape2=N[i-1]-sum(Z)+0.5)

      ## back transform to T1 and T2
      T1[i]<-T1.trans/(1-T1.trans)
      T2[i]<-T2.trans/(1-T2.trans)

      ## keep T1<T2
      if(T1[i]<T2[i])
        break
    }

    ## sample from p(N|A,T1,T2,Z,X)

    ## make sure N.new >=w.tau
    repeat {
      N.new<-rnbinom(1,mu=N[i-1],size=Metropolis.stdev.N)
      if(N.new>w.tau-1)
        break
    }

    ## calculate log(N.new!/(N.new-w.tau)!)
    N3.new<-rep(0,w.tau)
    for (j in 0:(w.tau-1)){
      N3.new[j+1]<-log(N.new-j)
    }
    N2.new<-sum(N3.new)

    ## calculate log(N[i-1]!/(N[i-1]-w.tau)!)
    N3<-rep(0,w.tau)
    for (j in 0:(w.tau-1)){
      N3[j+1]<-log(N[i-1]-j)
    }
    N2<-sum(N3)

    ## calculate log of acceptance ratio
    logr1<-(-1/2)*log(N.new)+N2.new+N.new*log(A[i]*(1/(1+T1[i]))+(1-A[i])*(1/(1+T2[i])))+(1/2)*log(N[i-1])-N2-N[i-1]*log(A[i]*(1/(1+T1[i]))+(1-A[i])*(1/(1+T2[i])))+log(dnbinom(N[i-1],mu=N.new,size=Metropolis.stdev.N))-log(dnbinom(N.new,mu=N[i-1],size=Metropolis.stdev.N))

    ## calculate acceptance ratio
    r1<-exp(logr1)

    ## accept or reject the proposed value
    if (runif(1)<min(r1,1)) {N[i]<-N.new ; a1<-a1+1} else {N[i]<-N[i-1]}

    ## calculate deviance from current sample

    ## calculate log(N[i]!/(N[i]-w.tau)!)
    N3.curr<-rep(0,w.tau)
    for (j in 0:(w.tau-1)){
      N3.curr[j+1]<-log(N[i]-j)
    }
    N2.curr<-sum(N3.curr)

    # calculate deviance
    D.post[i]<-(-2)*(N2.curr+(N[i]-w.tau)*log(A[i]*(1/(1+T1[i]))+(1-A[i])*(1/(1+T2[i])))+sum(freqdata[,2]*log(A[i]*(1/(1+T1[i]))*(T1[i]/(1+T1[i]))^freqdata[,1]+(1-A[i])*(1/(1+T2[i]))*(T2[i]/(1+T2[i]))^freqdata[,1]))-sum(lfactorial(freqdata[,2])))

  }



  ### Step 4: model diagnostics

  ## 1) deviance at posterior mean
  mean.A<-mean(A[(burn.in+1):iterations])
  mean.T1<-mean(T1[(burn.in+1):iterations])
  mean.T2<-mean(T2[(burn.in+1):iterations])
  mean.N<-mean(N[(burn.in+1):iterations])

  ## calculate log(mean.N!/(mean.N-w.tau)!)
  N3.mean<-rep(0,w.tau)
  for (j in 0:(w.tau-1)){
    N3.mean[j+1]<-log(mean.N-j)
  }
  N2.mean<-sum(N3.mean)

  loglik.post.mean<-N2.curr+(mean.N-w.tau)*log(mean.A*(1/(1+mean.T1))+(1-mean.A)*(1/(1+mean.T2)))+sum(freqdata[,2]*log(mean.A*(1/(1+mean.T1))*(mean.T1/(1+mean.T1))^freqdata[,1]+(1-mean.A)*(1/(1+mean.T2))*(mean.T2/(1+mean.T2))^freqdata[,1]))-sum(lfactorial(freqdata[,2]))
  D.mean<-(-2)*loglik.post.mean

  ## 2) posterior mean and median deviances
  mean.D<-mean(D.post[(burn.in+1):iterations])
  median.D<-quantile(D.post[(burn.in+1):iterations],probs=.5,names=F)

  ## 3) model complexity
  p.D<-mean.D-D.mean

  ## 4) Deviance information criterion
  DIC<-2*p.D+D.mean


  ### Step 5: fitted values based on medians of the marginal posteriors
  median.A<-quantile(A[(burn.in+1):iterations],probs=.5,names=F)
  median.T1<-quantile(T1[(burn.in+1):iterations],probs=.5,names=F)
  median.T2<-quantile(T2[(burn.in+1):iterations],probs=.5,names=F)
  median.N<-quantile(N[(burn.in+1):iterations],probs=.5,names=F)

  fits<-rep(0,tau)
  for (k in 1:tau){
    fits[k]<-(median.N)*(median.A*dgeom(k,prob=1/(1+median.T1))+(1-median.A)*dgeom(k,prob=1/(1+median.T2)))
  }
  fitted.values<-data.frame(cbind(j=seq(1,tau),fits,count=freqdata[,2]))


  ### Step 6: estimate thinning to reduce correlated posterior samples
  lags<-acf(N[(burn.in+1):iterations],type="correlation",main="Autocorr plot",ylab="ACF",xlab="Lag", plot=F)
  lag.thin<-suppressWarnings(min(which(lags$acf<0.1)))
  if (lag.thin==Inf) {lag.thin<-paste(">",length(lags$lag),sep="")
  }

  ### Step 7: results
  hist.points<-hist(N[(burn.in+1):iterations]+w-w.tau,breaks=seq(w,max(N)+w-w.tau+1)-0.5, plot = FALSE)

  results<-data.frame(w=w,
                      n=n,
                      NP.est.N=NP.est.n0+w,
                      tau=tau,
                      w.tau=w.tau,
                      n.tau=n.tau,
                      iterations=iterations-burn.in,
                      burn.in=burn.in,
                      acceptance.rate.T1T2=1,
                      acceptance.rate.N=a1/iterations,
                      lag=lag.thin,
                      mode.N=hist.points$mids[which.max(hist.points$density)],
                      mean.N=mean(N[(burn.in+1):iterations])+w-w.tau,
                      median.N=quantile(N[(burn.in+1):iterations]+w-w.tau,probs=.5,names=F),
                      LCI.N=quantile(N[(burn.in+1):iterations]+w-w.tau,probs=.025,names=F),
                      UCI.N=quantile(N[(burn.in+1):iterations]+w-w.tau,probs=.975,names=F),
                      stddev.N=sd((N[(burn.in+1):iterations]+w-w.tau)),
                      mean.D=mean.D,
                      median.D=median.D,
                      DIC
  )

  final_results <- list()
  final_results$est <- results$median.N
  final_results$ci <- c("lower 95%"=results$LCI.N, "upper 95%"=results$UCI.N)
  final_results$mean <- results$mean.N
  final_results$semeanest <- results$stddev.N
  final_results$dic  <- DIC
  final_results$fits <- fitted.values
  final_results$diagnostics<-c("acceptance rate N"=results$acceptance.rate.N,
                               "acceptance rate T1T2"=results$acceptance.rate.T1T2,
                               "lag"=results$lag)

  if (output) {
    # output results and fitted values
    print(final_results)
  }

  if (plot) {
    oldpar <- par(no.readonly = TRUE)
    par(mfrow=c(1,2))

    ## posterior histogram
    hist(N[(burn.in+1):iterations]+w-w.tau,
         main="Posterior distribution",xlab="Total Number of Species",
         col='purple',freq=F,ylab="Density")

    # make trace plot
    plot((burn.in+1):iterations,N[(burn.in+1):iterations]+w-w.tau,type="l",xlab="Iteration Number",ylab="Total Number of Species", main="Trace plot")
    on.exit(par(oldpar))
  }

  if (answers) {
    return(final_results)
  }

}

#' Estimate species richness with an objective Bayes method using a geometric model
#'
#' @param data A frequency count table
#' @param output Print summary of fitted model. Default is TRUE.
#' @param plot Plot estimates? Default is TRUE.
#' @param answers Return dataframe containing model summaries? Default is FALSE.
#' @param tau Truncation cutoff for model fitting. (Counts above this cutoff
#' are not used in model fitting.) Default is 10.
#' @param burn.in Number of MCMC burn-in iterations. Default is 100.
#' @param iterations Number of MCMC iterations after burn-in. Default is 2500.
#' @param Metropolis.stdev.N Standard deviation of proposal distribution
#' parameter N within MCMC. Default is 75.
#' @param Metropolis.start.theta Starting value of parameter theta.
#' Default is 1
#' @param Metropolis.stdev.theta Standard deviation of proposal distribution
#' parameter theta within MCMC. Default is 0.3.
#' @return If answers set to TRUE, a list containing posterior median and
#' 95% posterior credible intervals for richness, fitted values, and
#' MCMC diagnostics. If answers set to FALSE (default), no value is returned.
#' @importFrom graphics hist par plot
#'
#' @export
objective_bayes_geometric <- function(data,
                                      output=TRUE,
                                      plot=TRUE, answers=FALSE,
                                      tau=10, burn.in=100,
                                      iterations=2500,
                                      Metropolis.stdev.N=75,
                                      Metropolis.start.theta=1,
                                      Metropolis.stdev.theta=0.3) {

  data <- check_format(data)

  if (tau > max(data[,1])) {
    tau <- max(data[,1])
  }

  fullfreqdata  <- data
  # calculate NP estimate of n0
  w<-sum(fullfreqdata[,2])
  n<-sum(fullfreqdata[,1]*fullfreqdata[,2])

  # subset data below tau
  freqdata<-fullfreqdata[1:tau,]

  # calculate summary statistics and MLE estimate of n0 and C
  w.tau<-sum(freqdata[,2])
  n.tau<-sum(freqdata[,1]*freqdata[,2])
  R.hat<-(n.tau/w.tau-1)
  MLE.est.n0<-w.tau/R.hat
  MLE.est.N<-MLE.est.n0+w

  ### Step 3: calculate posterior

  ## initialization
  iterations <- iterations+burn.in
  N<-rep(0,iterations)
  R<-c(Metropolis.start.theta,rep(1,iterations-1))
  # to track acceptance rate of theta
  a1<-0
  # to track acceptance rate of N
  a2<-0
  # starting value based on MLE of C
  N[1]<-ceiling(MLE.est.N)
  D.post<-rep(0,iterations)

  for (i in 2:iterations){

    ## sample from p(theta|C,x)

    # propose value for theta
    R.new <- abs(rnorm(1, mean=R[i-1],
                       sd=Metropolis.stdev.theta))


    # calculate log of acceptance ratio
    logr1 <- (-N[i-1]-n.tau-1/2) * log(1+R.new) +
      (n.tau-1/2)*log(R.new) -
      (-N[i-1]-n.tau-1/2) * log(1+R[i-1]) -
      (n.tau-1/2)*log(R[i-1])


    # calculate acceptance ratio
    r1<-exp(logr1)

    # accept or reject propsed value
    if (runif(1)<min(r1,1)) {
      R[i]<-R.new ; a1<-a1+1
    }
    else {
      R[i]<-R[i-1]
    }
    ## sample from p(C|theta,x)

    ## make sure N.new >=w.tau
    repeat {
      N.new<-rnbinom(1,mu=N[i-1],size=Metropolis.stdev.N)
      if(N.new>w.tau-1)
        break
    }

    ## calculate log(N.new!/(N.new-w.tau)!)
    N3.new<-rep(0,w.tau)
    for (j in 0:(w.tau-1)){
      N3.new[j+1]<-log(N.new-j)
    }
    N2.new<-sum(N3.new)

    ## calculate log(N[i-1]!/(N[i-1]-w.tau)!)
    N3<-rep(0,w.tau)
    for (j in 0:(w.tau-1)){
      N3[j+1]<-log(N[i-1]-j)
    }
    N2<-sum(N3)

    # calculate log of acceptance ratio
    logr2<-(N2.new-(1/2)*log(N.new)-N.new*log(1+R[i]))-(N2-(1/2)*log(N[i-1])-N[i-1]*log(1+R[i]))+(log(dnbinom(N[i-1],mu=N.new,size=Metropolis.stdev.N)))-(log(dnbinom(N.new,mu=N[i-1],size=Metropolis.stdev.N)))

    # calculate acceptance ratio
    r2<-exp(logr2)


    # accept or reject propsed value
    if (runif(1)<min(r2,1)) {
      N[i]<-N.new ; a2<-a2+1
    }
    else {
      N[i]<-N[i-1]
    }
    ## calculate deviance from current sample

    # calculate log(N[i]!/(N[i]-w.tau)!)
    N3.curr<-rep(0,w.tau)
    for (j in 0:(w.tau-1)){
      N3.curr[j+1]<-log(N[i]-j)
    }
    N2.curr<-sum(N3.curr)

    # calculate deviance
    D.post[i]<-(-2)*(N2.curr-sum(log(factorial(freqdata[,2])))+n.tau*log(R[i])+(-N[i]-n.tau)*log(1+R[i]))

  }

  ### Step 4: model diagnostics

  ## 1) deviance at posterior mean
  mean.R<-mean(R[(burn.in+1):iterations])
  mean.N<-mean(N[(burn.in+1):iterations])

  ## calculate log(mean.N!/(mean.N-w.tau)!)
  N3.mean<-rep(0,w.tau)
  for (j in 0:(w.tau-1)){
    N3.mean[j+1]<-log(mean.N-j)
  }
  N2.mean<-sum(N3.mean)

  loglik.post.mean<-N2.mean-sum(log(factorial(freqdata[,2])))+n.tau*log(mean.R)-(mean.N+n.tau)*log(1+mean.R)
  D.mean<-(-2)*loglik.post.mean

  ## 2) posterior mean and median deviances
  mean.D<-mean(D.post[(burn.in+1):iterations])
  median.D<-quantile(D.post[(burn.in+1):iterations],probs=.5,names=F)

  ## 3) model complexity
  p.D<-mean.D-D.mean

  ## 4) Deviance information criterion
  DIC<-2*p.D+D.mean

  ### Step 5: fitted values based on medians of the marginal posteriors
  median.R<-quantile(R[(burn.in+1):iterations],probs=.5,names=F)
  median.N<-quantile(N[(burn.in+1):iterations],probs=.5,names=F)

  fits<-rep(0,tau)
  for (k in 1:tau){
    fits[k]<-(median.N)*dexp(k,1/(1+median.R))
  }
  fitted.values<-data.frame(cbind(j=seq(1,tau),fits,count=freqdata[,2]))




  ### Step 6: results
  hist.points<-hist(N[(burn.in+1):iterations]+w-w.tau,breaks=seq(w,max(N)+w-w.tau+1)-0.5, plot = plot)


  results<-data.frame(w=w,
                      n=n,
                      MLE.est.C=MLE.est.N,
                      tau=tau,
                      w.tau=w.tau,
                      n.tau=n.tau,
                      iterations=iterations,
                      burn.in=burn.in,
                      acceptance.rate.theta=a1/iterations,
                      acceptance.rate.N=a2/iterations,
                      mode.C=hist.points$mids[which.max(hist.points$density)],
                      mean.C=mean(N[(burn.in+1):iterations])+w-w.tau,
                      median.C=quantile(N[(burn.in+1):iterations]+w-w.tau,probs=.5,names=F),
                      LCI.C=quantile(N[(burn.in+1):iterations]+w-w.tau,probs=.025,names=F),
                      UCI.C=quantile(N[(burn.in+1):iterations]+w-w.tau,probs=.975,names=F),
                      stddev.C=sqrt(var((N[(burn.in+1):iterations]+w-w.tau))),
                      mean.D=mean.D,
                      median.D=median.D,
                      DIC
  )

  final_results <- list()
  final_results$results <- t(results)
  final_results$fits <- fitted.values

  if (output) {
    # output results and fitted values
    print(final_results)
  }

  if (plot) {
    oldpar <- par(no.readonly = TRUE)
    par(mfrow=c(2,2))

    # trace plot for C

    # first thin values of C if there are more than 10,000 iterations
    # must be a divisor of (iterations-burn.in)
    iterations.trace<-min(10000,iterations-burn.in)
    N.thin<-rep(0,iterations.trace)
    for (k in 1:iterations.trace){
      N.thin[k]<-N[k*((iterations-burn.in)/iterations.trace)]
    }

    # make trace plot
    plot(1:iterations.trace,N.thin,xlab="Iteration Number",ylab="Total Number of Species", main="Trace plot")

    # autocorrelation plot for C
    acf(N[(burn.in+1):iterations],type="correlation",main="Autocorr plot",ylab="ACF",xlab="Lag")

    # histogram of C with a bar for each discrete value
    hist(N[(burn.in+1):iterations]+w-w.tau,breaks=seq(w,max(N[
      (burn.in+1):iterations])+w-w.tau+1)-0.5,main="Posterior distribution",xlab="Total Number of Species",
      col='purple',freq=F,ylab="Density")
    on.exit(par(oldpar))
  }

  if (answers) {
    return(final_results)
  }

}

rmvnorm <- function (n, mean = rep(0, nrow(sigma)), sigma = diag(length(mean)),
                     method = c("eigen", "svd", "chol"), pre0.9_9994 = FALSE) {
  if (!isSymmetric(sigma, tol = sqrt(.Machine$double.eps),
                   check.attributes = FALSE)) {
    stop("sigma must be a symmetric matrix")
  }
  if (length(mean) != nrow(sigma))
    stop("mean and sigma have non-conforming size")
  method <- match.arg(method)
  R <- if (method == "eigen") {
    ev <- eigen(sigma, symmetric = TRUE)
    if (!all(ev$values >= -sqrt(.Machine$double.eps) * abs(ev$values[1]))) {
      warning("sigma is numerically not positive semidefinite")
    }
    t(ev$vectors %*% (t(ev$vectors) * sqrt(pmax(ev$values,
                                                0))))
  }
  else if (method == "svd") {
    s. <- svd(sigma)
    if (!all(s.$d >= -sqrt(.Machine$double.eps) * abs(s.$d[1]))) {
      warning("sigma is numerically not positive semidefinite")
    }
    t(s.$v %*% (t(s.$u) * sqrt(pmax(s.$d, 0))))
  }
  else if (method == "chol") {
    R <- chol(sigma, pivot = TRUE)
    R[, order(attr(R, "pivot"))]
  }
  retval <- matrix(rnorm(n * ncol(sigma)), nrow = n, byrow = !pre0.9_9994) %*%
    R
  retval <- sweep(retval, 2, mean, "+")
  colnames(retval) <- names(mean)
  retval
}
