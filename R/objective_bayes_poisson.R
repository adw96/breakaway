#' @export
objective_bayes_poisson <- function(data, print=TRUE, plot=TRUE, answers=FALSE, write=NULL,
                                    tau=10, burn.in=100, iterations=2500, Metropolis.stdev.N=75,
                                    Metropolis.start.lambda=1, Metropolis.stdev.lambda=0.3, bars=3) {

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
  L<-c(Metropolis.start.lambda,rep(1,iterations-1))
  # to track acceptance rate of lambda
  a1<-0
  # to track acceptance rate of C
  a2<-0
  # starting value based on nonparametric estimate of n0
  N[1]<-ceiling(NP.est.n0)+w.tau
  D.post<-rep(0,iterations)

  for (i in 2:iterations){

    ## sample from p(lambda|x,C)

    # propose value for lambda
    L.new<-abs(rnorm(1,mean=L[i-1],sd=Metropolis.stdev.lambda))

    # calculate log of acceptance ratio
    logr1<-(n.tau-1/2)*log(L.new)-L.new*N[i-1]-(n.tau-1/2)*log(L[i-1])+L[i-1]*N[i-1]

    # calculate acceptance ratio
    r1<-exp(logr1)

    # accept or reject propsed value
    if (runif(1)<min(r1,1)) {L[i]<-L.new ; a1<-a1+1}
    else L[i]<-L[i-1]

    ## sample from p(C|lambda,x)

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
    if (runif(1)<min(r2,1)) {N[i]<-N.new ; a2<-a2+1}
    else N[i]<-N[i-1]

    ## calculate deviance from current sample

    # calculate log(N[i]!/(N[i]-w.tau)!)
    N3.curr<-rep(0,w.tau)
    for (j in 0:(w.tau-1)){
      N3.curr[1:w.tau]<-log(N[i]-0:(w.tau-1))
    }
    N2.curr<-sum(N3.curr)

    # calculate deviance
    D.post[i]<-(-2)*(N2.curr-sum(log(factorial(freqdata[,2])))-L[i]*(N[i])-sum(freqdata[,2]*log(factorial(freqdata[,1])))+n.tau*log(L[i]))

  }

  ### Step 4: model diagnostics

  ## 1) deviance at posterior mean
  mean.L<-mean(L[(burn.in+1):iterations])
  mean.N<-mean(N[(burn.in+1):iterations])

  ## calculate log(mean.N!/(mean.N-w.tau)!)
  N3.mean<-rep(0,w.tau)
  N3.mean[1:w.tau]<-log(mean.N-0:(w.tau-1))
  N2.mean<-sum(N3.mean)

  loglik.post.mean<-N2.mean-sum(log(factorial(freqdata[,2])))-mean.L*mean.N+n.tau*log(mean.L)-sum(freqdata[,2]*log(factorial(freqdata[,1])))
  D.mean<-(-2)*loglik.post.mean

  ## 2) posterior mean and median deviances
  mean.D<-mean(D.post[(burn.in+1):iterations])
  median.D<-quantile(D.post[(burn.in+1):iterations],probs=.5,
                     names=F)

  ## 3) model complexity
  p.D<-mean.D-D.mean

  ## 4) Deviance information criterion
  DIC<-2*mean.D-D.mean

  ### Step 5: fitted values based on medians of the marginal posteriors
  median.L<-quantile(L[(burn.in+1):iterations],probs=.5,names=F)
  median.N<-quantile(N[(burn.in+1):iterations],probs=.5,names=F)

  fits<-rep(0,tau)
  fits[1:tau]<-(median.N)*dpois(1:tau,median.L)
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
                      acceptance.rate.lambda=a1/iterations,
                      acceptance.rate.N=a2/iterations,
                      mode.C=hist.points$mids[which.max(hist.points$density)],
                      mean.C=mean(N[(burn.in+1):iterations])+w-w.tau,
                      median.C=quantile(N[(burn.in+1):iterations]+w-w.tau,probs=.5,names=F),
                      LCI.C=quantile(N[(burn.in+1):iterations]+w-w.tau,probs=.025,names=F),
                      UCI.C=quantile(N[(burn.in+1):iterations]+w-w.tau,probs=.975,names=F),
                      stddev.C=sqrt(var((N[(burn.in+1):iterations]+w-w.tau))),
                      mean.D=mean.D,
                      median.D=median.D,
                      DIC)

  final_results <- list()
  final_results$est <- results$median.C
  final_results$mean <- results$mean.C
  final_results$semeanest <- sd(N[(burn.in+1):iterations])+w-w.tau
  final_results$ci <- c("lower"=results$LCI.C, "upper"=results$UCI.C)
  final_results$dic  <- DIC
  final_results$fits <- fitted.values
  final_results$diagnostics<-c("acceptance rate N"=results$acceptance.rate.N,"acceptance rate lambda"=results$acceptance.rate.lambda)

  if (print) {
    # print results and fitted values
    print(final_results)
  }

  if (!is.null(write)) {
    # output results and fitted values
    write.csv(t(results),write[1],quote=F)
    write.csv(fitted.values,write[2],quote=F)
  }

  if (plot != FALSE) {
    # prepare trace plot for C

    # first thin values of C if there are more than 10,000 iterations
    # must be a divisor of (iterations-burn.in)
    iterations.trace<-min(10000,iterations-burn.in)
    N.thin<-rep(0,iterations.trace)
    N.thin[1:iterations.trace]<-N[1:iterations.trace*((iterations-burn.in)/iterations.trace)]

    if (plot == "all") {
      plot.new()
      mat <- matrix(c(1,2,3,3), 2, byrow=T)
      layout(mat, c(1,1), c(1,1))

      # make trace plot
      plot(1:iterations.trace,N.thin,type="l",xlab="Iteration Number",ylab="Total Number of Species", main="Trace plot")

      # autocorrelation plot for C
      acf(N[(burn.in+1):iterations],type="correlation",main="Autocorr plot",ylab="ACF",xlab="Lag")

    } else {
      par(mfrow=c(2,1))
      # make trace plot
      plot(1:iterations.trace,N.thin,type="l",xlab="Iteration Number",ylab="Total Number of Species", main="Trace plot")
    }
    # # a histogram with # bars for each discrete value
    hist(N[(burn.in+1):iterations]+w-w.tau,
         breaks=seq(w,max(N[(burn.in+1):iterations])+w-w.tau+1,
                    length=(max(N[(burn.in+1):iterations])-w.tau+1)/bars+1)-0.5,
         main="Posterior distribution",xlab="Total Number of Species",
         col='purple',freq=F,ylab="Density",
         xlim=c(w,quantile(N[(burn.in+1):iterations]+w-w.tau,probs=.975,names=F)))
  }

  if (answers) {
    return(final_results)
  }

}
