\name{chao_bunge}
\alias{chao_bunge}
\title{ function for species richness estimation }
\description{This function implements the species richness estimation procedure outlined in Chao and Bunge (2002).}
\usage{chao_bunge(data, cutoff=10, print=TRUE, answers=FALSE)}
\arguments{
  \item{data}{
    The sample frequency count table for the population of interest. See \samp{breakaway} for details.
  }
  \item{cutoff}{
    The maximum frequency count to use for prediction. Defaults to 10.
  }
  \item{print}{
    Logical: whether the results should be printed to screen.
  }
 \item{answers}{
    Logical: whether the function should return an argument.
  }
}
\value{
  \item{name}{
    The ``name'' of the model.
  }
  \item{est}{
    The estimate of total (observed plus unobserved) diversity.
  }
  \item{seest}{
    The standard error in the diversity estimate.
  }
  \item{ci}{
    An asymmetric 95\% confidence interval for diversity.
  }
}
\references{
  Chao, A. and Bunge, J.  (2002). Estimating the number of species in a stochastic abundance model. \emph{Biometrics}, \bold{58}.
}
\author{Amy Willis}
\note{ This estimator is based on the negative binomial model and for that reason generally produces poor fits to microbial data. The result is artificially low standard errors. Caution is advised.
}
\seealso{ \code{\link{breakaway}}; \code{\link{apples}}; \code{\link{wlrm_untransformed}}; \code{\link{chao1}} }
\examples{chao_bunge(apples)}
\keyword{ diversity }
\keyword{ models }
