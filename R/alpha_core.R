#' @export
frequency_count_or_proportion_or_column <- function(input) {
  if (class(input) %in% c("data.frame", "matrix")) {
    if (ncol(input) == 2) {
      if (all(input[,2] %% 1 == 0)) {
        "frequency count"
      } else {
        "?"
      }
    } else {
      "?"
    }
  } else {
    if (all(input %% 1 == 0)) {
      "column"
    } else {
      "proportion"
    }
  }
}

#' convert to proportions
#' 
#' TODO write as generic method
#' 
#' @param input the input data
#' @param type What time
#' 
#' @export
to_proportions <- function(input, type) {
  
  if (type == "frequency count") {
    total_reads <- sum(input[, 2]*input[,1])
    relative_abundance <- input[, 1] / total_reads
    rep(relative_abundance, input[, 2])
  } else if (type == "proportion") {
    if (sum(input) != 1) warning("Proportions that don't sum to one; renormalising...")
    input <- input/sum(input)
    input[input>0]
  } else if (type == "column") {
    input <- input/sum(input)
    input[input>0]
  } else {
    stop("Error in to_proportions: unknown type")
  }
}













#' Plug-in Shannon index
#' 
#' Plug-in Shannon index, for population-level data.
#' 
#' 
#' @param input A frequency count table, vector of abundances, or vector of
#' proportions.
#' @return The Shannon index of the population given by input.
#' @note This function is intended for population-level data. If you are
#' dealing with a microbial sample, use DivNet instead.
#' @export shannon
shannon <- function(input) {
  if(any(is.na(input))) stop("Input contains missing values")
  type <- frequency_count_or_proportion_or_column(input)
  
  proportions <- to_proportions(input, type)
  
  -sum(proportions*log(proportions, base=exp(1)))
}















































#' Plug-in Hill numbers
#' 
#' Plug-in Hill numbers, for population-level data.
#' 
#' @param input A frequency count table, vector of abundances, or vector of
#' proportions.
#' @param q The Hill number of interest. q = 0 corresponds to species richness, q = 2 corresponds to inverse Simpson, etc.
#' @return The Hill number of the population given by input.
#' @note This function is intended for population-level data. If you are
#' dealing with a microbial sample, use DivNet instead.
#' @export hill
hill <- function(input, q) {
  
  if (length(q) > 1) {
    #  sapply(q, FUN = hill_under, input = input)
    sapply(q, FUN = hill, input = input)
  } else {
    if (q == 1) {
      exp(shannon(input))
    } else {
      type <- frequency_count_or_proportion_or_column(input)
      proportions <- to_proportions(input, type)
      (sum(proportions^q))^(1/(1-q))
    }
  }
}













































#' Plug-in Inverse Simpson index
#' 
#' Plug-in Inverse Simpson index, for population-level data.
#' 
#' @param input A frequency count table, vector of abundances, or vector of
#' proportions.
#' @return The inverse-Simpson index of the population given by input.
#' @note This function is intended for population-level data. If you are
#' dealing with a microbial sample, use DivNet instead.
#' @export inverse_simpson
inverse_simpson <- function(input) {
  hill(input, 2)
}














































#' Plug-in Simpson index
#' 
#' Plug-in Simpson index, for population-level data.
#' 
#' @param input A frequency count table, vector of abundances, or vector of
#' proportions.
#' @return The Simpson index of the population given by input.
#' @note This function is intended for population-level data. If you are
#' dealing with a microbial sample, use DivNet instead.
#' @export simpson
simpson <- function(input) {
  1/hill(input, 2)
}













































#' Plug-in Gini-Simpson index
#' 
#' Plug-in Gini-Simpson index, for population-level data.
#' 
#' @param input A frequency count table, vector of abundances, or vector of
#' proportions.
#' @return The Gini-Simpson index of the population given by input.
#' @note This function is intended for population-level data. If you are
#' dealing with a microbial sample, use DivNet instead.
#' @export gini
gini <- function(input) {
  1-simpson(input)
}

































