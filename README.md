# breakaway

The goal of breakaway is to do ALL OF ALPHA DIVERSITY the best. But for now, it's species richness estimation and modelling, with some other alpha diversity tools, for the high diversity (microbial) case.

## Installation

You can install breakaway from github with:

```R
# install.packages("devtools")
devtools::install_github("adw96/breakaway")
```

## Humans

Authors: Amy Willis, Kathryn Barger and John Bunge

Maintainer: Amy Willis (adwillis(at)uw(dot)edu)

## Example

Some quick examples (Amy: to return to**)

```R
breakaway(apples)
shannon(apples)

sd(replicate(50, resample_estimate(otus1[,5], shannon)))
```

## Notes from Amy Dec' '16:

This is the development version of the R package breakaway. It is being continually maintained so please use this version rather than the CRAN version.
I'm currently looking to expand the scope of breakaway from a species richness estimation procedure into a full one-stop shop for modelling alpha diversity. Key projects include (1) automated data structure detection and estimator selection, (2) coding CatchAll into R and making the source available, (3) incorporating uncertainties in OTU tables into alpha diversity estimates using quality scores, (4) updating base graphics procedures to ggplot, (5) implementing unit tests and turning this into a proper repository, (6) expanding the use cases of betta, (7) exploring methods of imputation for breakaway... the list goes on and you probably have ideas that I haven't even thought of! If you would like to be involved, please get in contact with me, or jump right in with your pull requests!
Please post your issues and questions here. I'm hoping to create a full Wiki with vignettes and bells and whistles soon, so any suggestions are appreciated.
Thank you for your patience with the current development version, and especially for logging your issues.
~Amy

## From CRAN: (needs to be updated)

Species richness estimation is an important problem in biodiversity analysis. This package provides methods for total species richness estimation (observed plus unobserved) and a method for modelling total diversity with covariates. `breakaway` estimates total (observed plus unobserved) species richness. Microbial diversity datasets are characterized by a large number of rare species and a small number of highly abundant species. The class of models implemented by `breakaway` is flexible enough to model both these features. `breakaway_nof1` implements a similar procedure however does not require a singleton count. `betta` provides a method for modelling total diversity with covariates in a way that accounts for its estimated nature and thus accounts for unobserved taxa, and `betta_random` permits random effects modelling.