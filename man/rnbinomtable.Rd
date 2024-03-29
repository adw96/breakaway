% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tools_simulate.R
\name{rnbinomtable}
\alias{rnbinomtable}
\title{Negative binomially distributed frequency count tables.}
\usage{
rnbinomtable(C, size, probability)
}
\arguments{
\item{C}{species richness}

\item{size}{size parameter for the negative binomial distribution}

\item{probability}{probability parameter for the negative binomial
distribution}
}
\value{
A simulated frequency count table.
}
\description{
Simulate a frequency count table based on a negative binomial model.
Zero-truncated, obviously.
}
\author{
Amy Willis
}
