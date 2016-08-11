# R package breakaway
Species richness estimation and modelling with high diversity, version 4.0

Last update: 2016-08-11

Authors: Amy Willis, Kathryn Barger and John Bunge

Maintainer: Amy Willis <adw96@cornell.edu>

(needs to be updated to reflect broader scope)

Species richness estimation is an important problem in biodiversity analysis. This package provides methods for total species richness estimation (observed plus unobserved) and a method for modelling total diversity with covariates. breakaway() estimates total (observed plus unobserved) species richness. Microbial diversity datasets are characterized by a large number of rare species and a small number of highly abundant species. The class of models implemented by breakaway() is flexible enough to model both these features. breakaway_nof1() implements a similar procedure however does not require a singleton count. betta() provides a method for modelling total diversity with covariates in a way that accounts for its estimated nature and thus accounts for unobserved taxa, and betta_random() permits random effects modelling.

License: GPL-2
