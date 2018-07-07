# testing and compiling breakaway
# directory <- "/Users/adwillis/software/breakaway"
# directory <- "/Users/adw96/Documents/software/breakaway"
directory <- "/Users/amy/Documents/software/breakaway"
# setwd(directory)
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

# check function can build
check() # builds namespace
# pkgdown::build_site()
rmarkdown::render("README.Rmd")
roxygenise(directory)
document(directory)
# build(directory)
build(directory, vignettes=F)
# build(directory, vignettes=T)
library(breakaway)
test(directory)

data("GlobalPatterns")
breakaway(GlobalPatterns)

(GlobalPatterns %>% 
  otu_table)[,1] %>% c %>% breakaway -> x
x
x$reasonable
x$other
x$plot + ylim(0,5)

?build_vignettes

#############################################
##### For development and testing
document(directory)
roxygenise(directory)
build(directory, vignettes=F)
library(breakaway)
breakaway(apples)
breakaway_nof1(apples[-1,]) %>% names
breakaway_nof1(hawaii[-1,])$estimate
breakaway_nof1(hawaii[-1,])$plot
breakaway_nof1(apples[-1,])$plot

breakaway(apples)$plot
test(directory)
check() # builds namespace

breakaway(toy_otu_table$Sample31)
breakaway(toy_otu_table$Sample4)
breakaway(toy_otu_table$Sample28)
chao_bunge(toy_otu_table$Sample31)
breakaway(apples)


#############################################


#check()
install(pkg = directory, build_vignettes=F)
library(breakaway)



data(apples)
breakaway(apples)

library(phyloseq)
library(ggplot2)

data(GlobalPatterns)
p <- plot_alpha(GlobalPatterns)
p
p + scale_y_continuous(limits=c(0, 10000))

plot_alpha(GlobalPatterns, x= "SampleType")


#
check()

roxygenise(directory)

devtools::load_all()
install()
build_site()
check()




objective_bayes_negbin(apples)


detach("breakaway", unload=TRUE)
unloadNamespace("breakaway")
breakaway::betta


### testing
setwd("/Users/adw96/Documents/Software/breakaway/")
library(devtools)
require(roxygen2)
require(testthat)
require(knitr)
require(rstudioapi)
require(devtools)
roxygen2::roxygenise()
devtools::build(vignettes = F)
devtools::install()
alpha_better(apples)
breakaway(hawaii)
hill(hawaii, q = c(0,1,2))


devtools::build(vignettes = F)
devtools::install()
hill_pic(hawaii)
hill_pic(hawaii, method = "delta")
hill_pic(hawaii, method = "delta?")
hill_pic(hawaii, method = "resample")
hill_pic(apples, method = "delta")
hill_pic(apples)
class(xx)
length(xx)
