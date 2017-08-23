#' Negative binomially distributed frequency count tables.
#' 
#' Simulate a frequency count table based on a negative binomial model. Zero-truncated, obviously.
#' 
#' 
#' 
#' @param C species richness
#' @param size size parameter for the negative binomial distribution
#' @param probability probability parameter for the negative binomial distribution
#' 
#' @author Amy Willis
#' 
#' @export betta_pic
rnbinomtable <- function(C, size, probability) {
  x <- rnbinom(n=C, size, probability)
  frequency_counts <- data.frame(table(x))
  colnames(frequency_counts) <- c("Index", "Frequency")
  rownames(frequency_counts) <- frequency_counts$Index
  frequency_counts$Index <- as.numeric(frequency_counts$Index)
  frequency_counts <- frequency_counts[frequency_counts$Index != 0, ]
  frequency_counts
}

#' beta version: Zero-truncated negative binomially distributed frequency count tables.
#' 
#' Simulate a frequency count table based on a negative binomial model. Zero-truncated, obviously.
#' 
#' 
#' 
#' @param C species richness
#' @param size size parameter for the negative binomial distribution
#' @param probability probability parameter for the negative binomial distribution
#' 
#' @author Amy Willis
#' 
rztnbinomtable <- function(C, size, probability) {
  # need to think about whether or not this is right
  x <- rnbinom(n=C, size, probability)
  frequency_counts <- data.frame(table(x))
  colnames(frequency_counts) <- c("Index", "Frequency")
  frequency_counts$Index <- as.numeric(frequency_counts$Index)
  frequency_counts
}
rztnbinomtable(50000, 500, 0.99)
