# testing and compiling breakaway
# directory <- "/Users/adwillis/software/breakaway"
directory <- "/Users/adw96/Documents/software/breakaway"
# directory <- "/Users/amy/Documents/software/breakaway"
setwd(directory)
library(devtools)
library(roxygen2)
library(knitr)
library(rstudioapi)
library(Rd2roxygen)
# install_github("r-lib/fs")
# devtools::install_github("r-lib/pkgdown")
# devtools::install_github("r-lib/testthat")
# devtools::install_github("HenrikBengtsson/R.rsp")
library(fs)
library(pkgdown)
library(testthat)
library(R.rsp)
library(covr)
library(magrittr)
library(ggplot2)
library(tibble)


check() 

document()
roxygenise()
build(vignettes = F)
library(breakaway)

# TODO fix plot(breakaway(GlobalPatterns))
# x <- breakaway(GlobalPatterns)
plot(x, symmetric = T) +  ylim(0, 1e5) ## why 6 warnings?
plot(x, symmetric = F) +  ylim(0, 1e5) ## why no intervals?
plot(x, symmetric = F) +  ylim(0, 1e5) +
  coord_cartesian(ylim = c(0,25000))
summary(x)$upper

test(directory)

covr::package_coverage()

# rmarkdown::render("README.Rmd")
# pkgdown::build_site()

