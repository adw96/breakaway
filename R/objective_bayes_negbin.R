#' @export
objective_bayes_negbin <- function(data, print=TRUE, plot=TRUE, answers=FALSE, write=FALSE,
                          tau=10, burn.in=1000, iterations=5000, Metropolis.stdev.N=100,
                          Metropolis.start.T1=-0.8, Metropolis.stdev.T1=0.05,
                          Metropolis.start.T2=0.8, Metropolis.stdev.T2=0.05, bars=3) {


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

  fullfreqdata  <- data
  # calculate NP estimate of n0
  w<-sum(fullfreqdata[,2])
  n<-sum(fullfreqdata[,1]*fullfreqdata[,2])
  NP.est.n0<-w/(1-fullfreqdata[1,2]/n)-w

  # subset data below tau
  freqdata<-fullfreqdata[1:tau,]

  # calculate summary statistics
  w.tau<-sum(freqdata[,2])
  n.tau<-sum(freqdata[,1]*freqdata[,2])
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

    ## sample from p(T1T2|N,x)

    ## propose value T1T2 from a bivariate normal dist.; make sure T1T2.new > {-1,0}
    repeat {
      T1T2.new <- MASS::mvrnorm(1, c(T1T2[i-1,1],T1T2[i-1,2]),matrix(c(Metropolis.stdev.T1,0,0,Metropolis.stdev.T2),nrow=2))
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
    if (runif(1)<min(r1,1)) {T1T2[i,]<-T1T2.new ; a1<-a1+1}
    else T1T2[i,]<-T1T2[i-1,]

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
    if (runif(1)<min(r2,1)) {N[i]<-N.new ; a2<-a2+1}
    else N[i]<-N[i-1]

    ## calculate deviance from current sample

    # calculate log(N[i]!/(N[i]-w.tau)!)
    N3.curr<-rep(0,w.tau)
    for (j in 0:(w.tau-1)){
      N3.curr[j+1]<-log(N[i]-j)
    }
    N2.curr<-sum(N3.curr)

    # calculate deviance
    D.post[i]<-(-2)*(N2.curr+n.tau*log(T1T2[i,2])-(N[i]*(1+T1T2[i,1])+n.tau)*log(1+T1T2[i,1]+T1T2[i,2])+N[i]*(1+T1T2[i,1])*log(1+T1T2[i,1])+sum(freqdata[,2]*lgamma(1+T1T2[i,1]+freqdata[,1]))-w.tau*lgamma(1+T1T2[i,1])-sum(log(factorial(freqdata[,2])))-sum(freqdata[,2]*lgamma(freqdata[,1]+1)))

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

  D.mean<-(-2)*(N2.mean+n.tau*log(mean.T2)-(mean.N*(1+mean.T1)+n.tau)*log(1+mean.T1+mean.T2)+mean.N*(1+mean.T1)*log(1+mean.T1)+sum(freqdata[,2]*lgamma(1+mean.T1+freqdata[,1]))-w.tau*lgamma(1+mean.T1)-sum(log(factorial(freqdata[,2])))-sum(freqdata[,2]*lgamma(freqdata[,1]+1)))

  ## 2) posterior mean and median deviances
  mean.D<-mean(D.post[(burn.in+1):iterations])
  median.D<-quantile(D.post[(burn.in+1):iterations],probs=.5,names=F)

  ## 3) model complexity
  p.D<-mean.D-D.mean

  ## 4) Deviance information criterion
  DIC<-2*mean.D-D.mean



  ### Step 5: fitted values based on medians of the marginal posteriors
  median.T1<-quantile(T1T2[(burn.in+1):iterations,1],probs=.5,names=F)
  median.T2<-quantile(T1T2[(burn.in+1):iterations,2],probs=.5,names=F)
  median.N<-quantile(N[(burn.in+1):iterations],probs=.5,names=F)

  fits<-rep(0,tau)
  for (k in 1:tau){
    fits[k]<-(median.N)*dnbinom(k,size=median.T1+1,prob=(median.T1+1)/(median.T1+median.T2+1))
  }
  fitted.values<-data.frame(cbind(j=seq(1,tau),fits,count=freqdata[,2]))


  ### Step 6: results
  if (plot) hist.points<-hist(N[(burn.in+1):iterations]+w-w.tau,breaks=seq(w,max(N)+w-w.tau+1)-0.5)

  results<-data.frame(w=w,
                      n=n,
                      NP.est.C=NP.est.n0+w,
                      tau=tau,
                      w.tau=w.tau,
                      n.tau=n.tau,
                      iterations=iterations,
                      burn.in=burn.in,
                      acceptance.rate.T1T2=a1/iterations,
                      acceptance.rate.N=a2/iterations,
                      mode.N=hist.points$mids[which.max(hist.points$density)],
                      mean.N=mean(N[(burn.in+1):iterations])+w-w.tau,
                      median.N=quantile(N[(burn.in+1):iterations]+w-w.tau,probs=.5,names=F),
                      LCI.N=quantile(N[(burn.in+1):iterations]+w-w.tau,probs=.025,names=F),
                      UCI.N=quantile(N[(burn.in+1):iterations]+w-w.tau,probs=.975,names=F),
                      stddev.C=sqrt(var((N[(burn.in+1):iterations]+w-w.tau))),
                      mean.D=mean.D,
                      median.D=median.D,
                      DIC
  )

  final_results <- list()
  final_results$results <- t(results)
  final_results$fits <- fitted.values

  if (print) {
    # print results and fitted values
    print(final_results)
  }

  if (write) {
    # output results and fitted values
    write.csv(t(results),filename.analysis,quote=F)
    write.csv(fitted.values,filename.fits,quote=F)
  }

  if (plot) {
    plot.new()
    mat <- matrix(c(1,2,3,3), 2, byrow=T)
    layout(mat, c(1,1), c(1,1))

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
    hist(N[(burn.in+1):iterations]+w-w.tau,
         breaks=seq(w,max(N[(burn.in+1):iterations])+w-w.tau+1)-0.5,
         main="Posterior distribution",xlab="Total Number of Species",
         col='purple',freq=F,ylab="Density")

    # # a histogram with # bars for each discrete value
    # hist(N[(burn.in+1):iterations]+w-w.tau,
    #      breaks=seq(w,max(N[(burn.in+1):iterations])+w-w.tau+1,
    #      length=(max(N[(burn.in+1):iterations])-w.tau+1)/bars+1)-0.5,
    #      main="Posterior distribution",xlab="Total Number of Species",
    #      col='purple',freq=F,ylab="Density",
    #      xlim=c(w,quantile(N[(burn.in+1):iterations]+w-w.tau,probs=.975,names=F)))
  }

  if (answers) {
    return(final_results)
  }

}
