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



















#' Plug-in Hill numbers
#' 
#' TODO
#' 
#' 
#' @param input TODO
#' @param q TODO
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

















#' Plug-in Hill numbers
#' 
#' TODO
#' 
#' 
#' @param input TODO
#' @export inverse_simpson
inverse_simpson <- function(input) {
  hill(input, 2)
}



















#' Plug-in Simpson
#' 
#' TODO
#' 
#' 
#' @param input TODO
#' @export simpson
simpson <- function(input) {
  1/hill(input, 2)
}

















#' Plug-in Gini-Simpson
#' 
#' TODO
#' 
#' 
#' @param input TODO
#' @export gini
gini <- function(input) {
  1-simpson(input)
}



















#' alpha diversity estimates
#' 
#' TODO
#' 
#' 
#' @param input TODO
#' @param q TODO
#' @param ccc TODO
#' @export alpha_better
alpha_better <-  function(input, q = 0, ccc = NA, ccc_se = NA) {
  
  type <- frequency_count_or_proportion_or_column(input)
  proportions <- to_proportions(input, type)
  
  if ((is.na(ccc) | is.null(ccc)) & type == "frequency count") {
    baway <- breakaway(input, print = F, answers = T, plot = F)
    ccc <- round(baway$est)
    ccc_se <- round(baway$seest)
  }
  
  cc <- sum(input[,2])
  unobs <- ccc-cc
  
  unobs_props <- rep(1/ccc, unobs)
  obs_props <- cc/ccc * proportions
  new_props <- c(unobs_props, obs_props)
  
  if (length(q) > 1) {
    estimates <- sapply(X = q, FUN = hill, input = new_props)
    ses <- mapply(FUN = hill_se, 
                  q = q, 
                  Dest = estimates,
                  MoreArgs = list(Cest = ccc, Cse = ccc_se,
                                  props = obs_props, c = cc))
    
  } else {
    estimates <- hill(new_props, q)
    ses <- hill_se(Cest = ccc, Cse = ccc_se, 
                   q, props = obs_props, Dest = estimates, 
                   c = cc)
  }
  data.frame("q" = q, "Estimate" = estimates, "StdError" = ses)
  
}

















#' alpha diversity std errors
#' 
#' TODO
#' 
#' 
#' @param Cest TODO
#' @param Cse TODO
#' @param q TODO
#' @param props TODO
#' @param Dest TODO
#' @param c TODO
#' @export hill_se
hill_se <- function(Cest, Cse, q, props, Dest, c) {
  
  derivative <- (Dest^q) / (1-q) * 
    ((1-q) * Cest^-q + c * q * Cest^(-q-1)  - q * c^q * Cest^(-q-1) * sum(props^q))
  
  variance_est <- derivative ^ 2 * Cse^2 
  sqrt(variance_est)
}
