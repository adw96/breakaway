\name{chao1}
\alias{chao1}
\title{ species richness lower bound }
\description{This function implements the Chao1 diversity index.}
\usage{chao1(data, print=TRUE, answers=FALSE)}
\arguments{
  \item{data}{
    The sample frequency count table for the population of interest. See \samp{breakaway} for details.
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
  Chao, A. (1984). Nonparametric estimation of the number of classes in a population. \emph{Scandinavian Journal of Statistics}, \bold{4}.
}
\author{Amy Willis}
\note{ This estimator is based on the assumption of equal species detection probabilities, that is,  every group in the population is equally likely to be observed. This assumption never applies, thus this should never be used for species richness estimation. However, it provides a reliable \emph{lower bound} for species richness, and has historical significance. \samp{breakaway}, \samp{breakaway_nof1} and \samp{chao_bunge} provide viable alternatives for species richness estimation and are included in this package.
}
\seealso{ \code{\link{breakaway}}; \code{\link{apples}}; \code{\link{wlrm_untransformed}}; \code{\link{chao1_bc}}; \code{\link{chao_bunge}} }
\examples{chao1(apples)}
\keyword{ diversity }
\keyword{ models }
