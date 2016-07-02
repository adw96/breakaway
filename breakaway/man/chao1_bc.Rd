\name{chao1_bc}
\alias{chao1_bc}
\title{ species richness estimator under equal detection probabilities }
\description{This function implements the bias-corrected Chao1 species richness estimator.}
\usage{chao1_bc(data, print=TRUE, answers=FALSE)}
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
\note{ This estimator is based on the assumption of equal species detection probabilities, that is,  every group in the population is equally likely to be observed. This assumption never applies, thus this estimator should never be used for species richness estimation. \samp{breakaway}, \samp{breakaway_nof1} and \samp{chao_bunge} provide viable alternatives and are included in this package.}
\seealso{ \code{\link{breakaway}}; \code{\link{apples}}; \code{\link{wlrm_untransformed}}; \code{\link{chao1}}; \code{\link{chao_bunge}}}
\examples{chao1_bc(apples)}
\keyword{ diversity }
\keyword{ models }
