objective_bayes_geometric <- function(data, print=TRUE, plot=TRUE, answers=FALSE, write=FALSE,
                           tau=10, burn.in=100, iterations=2500, Metropolis.stdev.N=75,
                           Metropolis.start.theta=1, Metropolis.stdev.theta=0.3) {

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
  iterations<-iterations+burn.in
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
    R.new<-abs(rnorm(1,mean=R[i-1],sd=Metropolis.stdev.theta))

    # calculate log of acceptance ratio
    logr1<-(-N[i-1]-n.tau-1/2)*log(1+R.new)+(n.tau-1/2)*log(R.new)-(-N[i-1]-n.tau-1/2)*log(1+R[i-1])-(n.tau-1/2)*log(R[i-1])

    # calculate acceptance ratio
    r1<-exp(logr1)

    # accept or reject propsed value
    if (runif(1)<min(r1,1)) {R[i]<-R.new ; a1<-a1+1}
    else R[i]<-R[i-1]

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
  DIC<-2*mean.D-D.mean

  ### Step 5: fitted values based on medians of the marginal posteriors
  median.R<-quantile(R[(burn.in+1):iterations],probs=.5,names=F)
  median.N<-quantile(N[(burn.in+1):iterations],probs=.5,names=F)

  fits<-rep(0,tau)
  for (k in 1:tau){
    fits[k]<-(median.N)*dexp(k,1/(1+median.R))
  }
  fitted.values<-data.frame(cbind(j=seq(1,tau),fits,count=freqdata[,2]))




  ### Step 6: results
  if (plot) hist.points<-hist(N[(burn.in+1):iterations]+w-w.tau,breaks=seq(w,max(N)+w-w.tau+1)-0.5)


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
  }

  if (answers) {
    return(final_results)
  }

}
