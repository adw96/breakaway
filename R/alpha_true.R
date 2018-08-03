is_proportions <- function(input) {
  ((sum(input) - 1)^2) < 10^-8
}

#' Calculate the true Shannon index based on proportions
#' 
#' @param input A vector of proportions.
#' 
#' @return The Shannon index of the population given by input.
#' @note This function is intended for population-level data. If you are
#' dealing with a microbial sample, use DivNet instead.
#' @export 
true_shannon <- function(input) {
  
  if (!is_proportions(input)) {
    stop("Shannon can only be calculated on proportions.\nUse sample_shannon() instead.\n")
  } 

  -sum(input*log(input, base=exp(1)))
}


#' Calculate the true Hill numbers
#' 
#' @param input A vector of proportions.
#' @param q The Hill number of interest. q = 0 corresponds to species richness, q = 2 corresponds to inverse Simpson, etc.
#' @return The Hill number of the population given by input.
#' @export 
true_hill <- function(input, q) {
  
  if (!is_proportions(input)) {
    stop("Hill numbers can only be calculated on proportions.\nUse sample_shannon() instead.\n")
  } 
  
  if (q == 1) {
    hh <- exp(true_shannon(input))
  } else {
    hh <- (sum(input^q))^(1/(1-q))
  }
  hh
}

#' Calculate the true Inverse Simpson index
#' 
#' 
#' @param input A vector of proportions.
#' @return The inverse-Simpson index of the population given by input.
#' @note This function is intended for population-level data. If you are
#' dealing with a microbial sample, use DivNet instead.
#' @export 
true_inverse_simpson <- function(input) {
  true_hill(input, 2)
}

#' Calculate the true Simpson index
#' 
#' @param input A vector of proportions.
#' @return The Simpson index of the population given by input.
#' @note This function is intended for population-level data. If you are
#' dealing with a microbial sample, use DivNet instead.
#' @export 
true_simpson <- function(input) {
  1/true_hill(input, 2)
}

#' Calculate the true Gini-Simpson index
#' 
#' @param input A vector of proportions.
#' @return The Gini-Simpson index of the population given by input.
#' @note This function is intended for population-level data. If you are
#' dealing with a microbial sample, use DivNet instead.
#' @export 
true_gini <- function(input) {
  1-true_simpson(input)
}

#' Calculate the true Shannon's equitability index
#' 
#' @param input A vector of proportions.
#' @return The Shannon E's of the population given by input.
#' @note This function is intended for population-level data. If you are
#' dealing with a microbial sample, use DivNet instead.
#' @export 
true_shannon_e <- function(input) {
  
  if (!is_proportions(input)) {
    stop("Shannon can only be calculated on proportions.\nUse sample_shannon() instead.\n")
  } 
  
  true_shannon(input) / log(sum (input > 0))

}