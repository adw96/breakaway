# testing and compiling breakaway
# directory <- "/Users/adwillis/software/breakaway"
# directory <- "/Users/adw96/Documents/software/breakaway"
directory <- "/Users/amy/Documents/software/breakaway"
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

library(phyloseq)
data("GlobalPatterns")
phyloseq::plot_richness(GlobalPatterns, x="SampleType", color="SampleType")

yy <- breakaway(GlobalPatterns)
plot(yy, data = GlobalPatterns,
     xaxis = "SampleType", color = "SampleType")

plot(x = chao_bunge(GlobalPatterns), physeq = GlobalPatterns,
     xaxis = "SampleType", color = "SampleType")

test(directory)

cov <- covr::package_coverage()
plot(cov)
zero_coverage(cov)

# rmarkdown::render("README.Rmd")
# pkgdown::build_site()


breakaway(apples) %>% class
breakaway(GlobalPatterns) %>% class
sample_richness(apples)
sample_richness(GlobalPatterns)

