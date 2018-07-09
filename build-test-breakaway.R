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


roxygenise(directory)
build(directory)
library(breakaway)

test(directory)

covr::package_coverage()

# rmarkdown::render("README.Rmd")
# pkgdown::build_site()

