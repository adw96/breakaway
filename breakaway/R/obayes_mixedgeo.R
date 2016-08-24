obayes_mixedgeo <- function(data, print=TRUE, plot=TRUE, answers=FALSE, write=FALSE,
                             tau=10, burn.in=100, iterations=2500, Metropolis.stdev.N=100,
                             Metropolis.start.T1=1, Metropolis.stdev.T1=2,
                             Metropolis.start.T2=3, Metropolis.stdev.T2=2, bars=3) {

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
    if (runif(1)<min(r1,1)) {N[i]<-N.new ; a1<-a1+1}
    else N[i]<-N[i-1]

    ## calculate deviance from current sample

    ## calculate log(N[i]!/(N[i]-w.tau)!)
    N3.curr<-rep(0,w.tau)
    for (j in 0:(w.tau-1)){
      N3.curr[j+1]<-log(N[i]-j)
    }
    N2.curr<-sum(N3.curr)

    # calculate deviance
    D.post[i]<-(-2)*(N2.curr+(N[i]-w.tau)*log(A[i]*(1/(1+T1[i]))+(1-A[i])*(1/(1+T2[i])))+sum(freqdata[,2]*log(A[i]*(1/(1+T1[i]))*(T1[i]/(1+T1[i]))^freqdata[,1]+(1-A[i])*(1/(1+T2[i]))*(T2[i]/(1+T2[i]))^freqdata[,1]))-sum(log(factorial(freqdata[,2]))))

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

  D.mean<-(-2)*(N2.curr+(mean.N-w.tau)*log(mean.A*(1/(1+mean.T1))+(1-mean.A)*(1/(1+mean.T2)))+sum(freqdata[,2]*log(mean.A*(1/(1+mean.T1))*(mean.T1/(1+mean.T1))^freqdata[,1]+(1-mean.A)*(1/(1+mean.T2))*(mean.T2/(1+mean.T2))^freqdata[,1]))-sum(log(factorial(freqdata[,2]))))

  ## 2) posterior mean and median deviances
  mean.D<-mean(D.post[(burn.in+1):iterations])
  median.D<-quantile(D.post[(burn.in+1):iterations],probs=.5,names=F)

  ## 3) model complexity
  p.D<-mean.D-D.mean

  ## 4) Deviance information criterion
  DIC<-2*mean.D-D.mean



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

  ### Step 6: results
  hist.points<-hist(N[(burn.in+1):iterations]+w-w.tau,breaks=seq(w,max(N)+w-w.tau+1)-0.5)

  results<-data.frame(w=w,
                      n=n,
                      NP.est.C=NP.est.n0+w,
                      tau=tau,
                      w.tau=w.tau,
                      n.tau=n.tau,
                      iterations=iterations,
                      burn.in=burn.in,
                      acceptance.rate.T1T2=1,
                      acceptance.rate.N=a1/iterations,
                      mode.N=hist.points$mids[which.max(hist.points$density)],
                      mean.N=mean(N[(burn.in+1):iterations])+w-w.tau,
                      median.N=quantile(N[(burn.in+1):iterations]+w-w.tau,probs=.5,names=F),
                      LCI.N=quantile(N[(burn.in+1):iterations]+w-w.tau,probs=.025,names=F),
                      UCI.N=quantile(N[(burn.in+1):iterations]+w-w.tau,probs=.975,names=F),
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
    hist(N[(burn.in+1):iterations]+w-w.tau,
         breaks=seq(w,max(N[(burn.in+1):iterations])+w-w.tau+1)-0.5,
         main="Posterior distribution",xlab="Total Number of Species",
         col='purple',freq=F,ylab="Density")

    # a histogram with # bars for each discrete value
    hist(N[(burn.in+1):iterations]+w-w.tau,
         breaks=seq(w,max(N[(burn.in+1):iterations])+w-w.tau+1,
                    length=(max(N[(burn.in+1):iterations])-w.tau+1)/bars+1)-0.5,
         main="Posterior distribution",xlab="Total Number of Species",col='purple',
         freq=F,ylab="Density",xlim=c(w,quantile(N[(burn.in+1):iterations]+w-w.tau,probs=.975,names=F)))
  }

  if (answers) {
    return(final_results)
  }

}
