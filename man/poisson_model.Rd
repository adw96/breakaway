% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/richness_poisson.R
\name{poisson_model}
\alias{poisson_model}
\title{PoissonModel}
\usage{
poisson_model(input_data, cutoff = 10)
}
\arguments{
\item{input_data}{A frequency count table}

\item{cutoff}{The largest frequency to use for predicting f0. Defaults to 10.}
}
\value{
An object of class \code{alpha_estimate}.
}
\description{
A model to estimate the number of missing taxa under a Poisson Model
}
