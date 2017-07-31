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

#' @export
to_proportions <- function(input, type) {
  
  if (type == "frequency count") {
    total_reads <- sum(input[, 2]*input[,1])
    relative_abundance <- input[, 1] / total_reads
    rep(relative_abundance, input[, 2])
  } else if (type == "proportion") {
    input[input>0]
  } else if (type == "column") {
    input <- input/sum(input)
    input[input>0]
  } else {
    stop("Error in to_proportions: unknown type")
  }
}

#' @export
shannon <- function(input) {
  
  type <- frequency_count_or_proportion_or_column(input)
  proportions <- to_proportions(input, type)
  
  -sum(proportions*log(proportions))
}

hill_under <- function(input, q) {
  type <- frequency_count_or_proportion_or_column(input)
  proportions <- to_proportions(input, type)
  
  (sum(proportions^q))^(1/(1-q))
  
}

#' @title  Plug-in Hill numbers
#' 
#' @description  TODO
#' 
#' @param input TODO
#' @param q TODO
#' @export
hill <- function(input, q) {
  
  if (length(q) > 1) {
    sapply(q, FUN = hill_under, input = input)
  } else {
    hill_under(input, q)
  }
}

#' @title  Plug-in Hill numbers
#' 
#' @description  TODO
#' 
#' @param input TODO
#' @export
inverse_simpson <- function(input) {
  hill(input, 2)
}

#' @title  Plug-in Simpson
#' 
#' @description  TODO
#' 
#' @param input TODO
#' @export
simpson <- function(input) {
  1/hill(input, 2)
}

#' @title  Plug-in Gini-Simpson
#' 
#' @description  TODO
#' 
#' @param input TODO
#' @export
gini <- function(input) {
  1-simpson(input)
}

#' @title  alpha diversity estimates
#' 
#' @description  TODO
#' 
#' @param input TODO
#' @param q TODO
#' @param ccc TODO
#' @export
alpha_better <-  function(input, q, ccc = NA) {
    
  type <- frequency_count_or_proportion_or_column(input)
  proportions <- to_proportions(input, type)
  
  if ((is.na(ccc) | is.null(ccc)) & type == "frequency count") {
    baway <- breakaway(input, print = F, answers = T, plot = F)
    ccc <- round(baway$est)
  }
  
  cc <- sum(hawaii[,2])
  unobs <- ccc-cc
  
  unobs_props <- rep(1/ccc, unobs)
  obs_props <- cc/ccc * proportions
  new_props <- c(unobs_props, obs_props)
  
  hill(new_props, q)
}