\name{breakaway-package}
\alias{breakaway-package}
\docType{package}
\title{ Species richness estimation and modelling in the high-diversity setting }
\description{ Species richness estimation is an important problem in biodiversity analysis. This package provides methods for total species richness estimation (observed plus unobserved) and a method for modelling total diversity with covariates. \samp{breakaway} estimates total (observed plus unobserved) species richness. Microbial diversity datasets are characterized by a large number of rare species and a small number of highly abundant species. The class of models implemented by \samp{breakaway} is flexible enough to model both these features. \samp{breakaway_nof1} implements a similar procedure however does not require a singleton count. \samp{betta} provides a method for modelling total diversity with covariates in a way that accounts for its estimated nature and thus accounts for unobserved taxa, and \samp{betta_random} permits random effects modelling. }
\details{
\tabular{ll}{
Package: \tab breakaway\cr
Type: \tab Package\cr
Version: \tab 3.0\cr
Date: \tab 2016-03-29\cr
License: \tab GPL-2 \cr
}
The function \code{\link{breakaway}} estimates the total (observed plus unobserved) number of classes (usually, distinct species) based on a sample of the frequency counts. Standard errors and model fits are also given. The algorithm is based on theory of characterization of distributions by ratios of their probabilities. Parameter estimation is done via nonlinear regression. The class of models available is usually broad enough to account for the high-diversity case, which is often observed in microbial diversity datasets. Since many classical estimation procedures either fail to provide an estimate or provide poor fits in the microbial setting, \code{\link{breakaway}} addresses this data structure. Additionally, since sequencing errors may result in an inflated singleton count, \code{\link{breakaway_nof1}} performs a similar procedure but does not require a singleton count. It can be used as an exploratory tool for investigating the plausibility of the given singleton count. \code{\link{betta}} runs a regression-type analysis of estimated total diversity, thus permitting accounting for unobserved taxa. It does not enforce use of \code{\link{breakaway}} for diversity estimation. A mixed-model approach accounts for the differing levels of confidence in the diversity estimates, and covariates constitute the fixed effects. Support of this work from Cornell University's Department of Statistical Sciences is gratefully acknowledged.
}
\author{
Amy Willis & John Bunge

Maintainer: Amy Willis <adw96@cornell.edu>
}
\references{
	Willis, A. and Bunge, J. (2015). Estimating diversity via frequency ratios. \emph{Biometrics.}

	  Willis, A. (2015). Species richness estimation with high diversity but spurious singletons. \emph{Under review.}

	Willis, A., Bunge, J., and Whitman, T. (2015). Inference for changes in biodiversity. \emph{arXiv preprint.}

	Rocchetti, I., Bunge, J. and Bohning, D. (2011). Population size estimation based upon ratios of recapture probabilities. \emph{Annals of Applied Statistics}, \bold{5}.

  Chao, A. and Bunge, J.  (2002). Estimating the number of species in a stochastic abundance model. \emph{Biometrics}, \bold{58}.

  Chao, A. (1984). Nonparametric estimation of the number of classes in a population. \emph{Scandinavian Journal of Statistics}, \bold{4}.
}
\examples{
breakaway(apples)
breakaway(apples,plot=FALSE,print=FALSE,answers=TRUE)
breakaway_nof1(apples[-1,])
betta(c(1552,1500,884),c(305,675,205))
}
