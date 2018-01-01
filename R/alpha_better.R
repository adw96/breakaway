#' Hill number estimates
#' 
#' In many ecological surveys, only a small fraction of the total ecosystem is
#' observed. In order to draw conclusions about the total ecosystem, it is
#' necessary to adjust for inexhaustive sampling in estimates of alpha
#' diversity. This function implements such an adjustment for the Hill numbers.
#' 
#' 
#' @param input A frequency count table or vector of abundances.
#' @param q The index of the Hill number of interest. Vectors are acceptable.
#' @param ccc An estimate of the total species richness of the community under
#' study. If none is given, \code{\link{breakaway}} is used
#' @param ccc_se The standard error in the estimate of ccc.
#' @param quantile The desired coverage of the confidence interval. Defaults to 0.95. Note all confidence intervals are approximate.
#' @return An estimate, standard error, and confidence set of the Hill number.
#' @examples 
#'
#' hill_better(apples, q = 2)
#' 
#' @export hill_better
hill_better <-  function(input, q = 0, ccc = NA, ccc_se = NA, quantile = 0.95) {
  
  type <- frequency_count_or_proportion_or_column(input)
  proportions <- to_proportions(input, type)
  
  if ((is.na(ccc) | is.null(ccc)) & type == "frequency count") {
    baway <- breakaway(input,  output = F, answers = T, plot = F)
    ccc <- round(baway$est)
    ccc_se <- round(baway$seest)
  }
  
  cc <- sum(input[,2])
  unobs <- ccc-cc
  
  unobserved_proportions <- rep(1/ccc, unobs)
  observed_proportions <- cc/ccc * proportions
  amended_proportions <- c(unobserved_proportions, observed_proportions)
  
  if (length(q) > 1) {
    estimates <- sapply(X = q, FUN = hill, input = amended_proportions)
    ses <- mapply(FUN = hill_se, 
                  q = q, 
                  Dest = estimates,
                  MoreArgs = list(Cest = ccc, Cse = ccc_se,
                                  proportions = observed_proportions, c = cc))
    
  } else {
    estimates <- hill(amended_proportions, q)
    ses <- hill_se(Cest = ccc, Cse = ccc_se, 
                   q, proportions = observed_proportions, Dest = estimates, 
                   c = cc)
  }
  
  data.frame("index"  = q, "estimate" = estimates, "standard_error" = ses, 
             "lower" = max(0, estimates-qnorm(1-(1-quantile)/2)*ses), 
             "upper" = estimates+qnorm(1-(1-quantile)/2)*ses)
}






hill_se <- function(Cest, Cse, q, proportions, Dest, c) {
  
  derivative <- (Dest^q) / (1-q) * 
    ((1-q) * Cest^-q + c * q * Cest^(-q-1)  - q * c^q * Cest^(-q-1) * sum(proportions^q))
  
  variance_est <- derivative ^ 2 * Cse^2 
  sqrt(variance_est)
}











#' Shannon diversity estimates
#' 
#' In many ecological surveys, only a small fraction of the total ecosystem is
#' observed. In order to draw conclusions about the total ecosystem, it is
#' necessary to adjust for inexhaustive sampling in estimates of alpha
#' diversity. This function implements such an adjustment for the Shannon
#' index.
#' 
#' 
#' @param input A frequency count table or vector of abundances.
#' @param ccc An estimate of the total species richness of the community under
#' study. If none is given, \code{\link{breakaway}} is used
#' @param ccc_se The standard error in the estimate of ccc.
#' @param quantile The desired coverage of the confidence interval. Defaults to 0.95. Note all confidence intervals are approximate.
#' @return An estimate, standard error, and confidence set of the index.
#' @examples
#' 
#' 
#' 
#' 
#' shannon_better(apples)
#' 
#' 
#' 
#' 
#' @export shannon_better
shannon_better <- function(input, ccc = NA, ccc_se = NA, quantile = 0.95) {
  
  if (length(dim(input)) != 2) {
    stop("Incorrent format: I need a frequency count table")
  } 
  input <- check_format(input)
  
  type <- frequency_count_or_proportion_or_column(input)
  proportions <- to_proportions(input, type)
  
  nn <- sum(input[, 1]*input[, 2])
  
  if ((is.na(ccc) | is.null(ccc))) {
    baway <- breakaway(input,  output = F, answers = T, plot = F)
    ccc <- round(baway$est)
    ccc_se <- round(baway$seest)
  }
  
  cc <- sum(input[,2])
  unobs <- round(ccc-cc)

  unobserved_proportions <- rep(1/ccc, unobs)
  observed_proportions <- cc/ccc * proportions
  amended_proportions <- c(unobserved_proportions, observed_proportions)
  
  estimates <- shannon(amended_proportions)
  ses <- sqrt(sum((1 + log(amended_proportions))^2 * amended_proportions * (1-amended_proportions)/nn) )
  
  data.frame("index"  = "shannon", 
             "estimate" = estimates, 
             "standard_error" = ses, 
             "lower" = max(0, estimates-qnorm(1-(1-quantile)/2)*ses), 
             "upper" = estimates+qnorm(1-(1-quantile)/2)*ses)
  
}


























#' An estimate of the Inverse-Simpson index
#' 
#' In many ecological surveys, only a small fraction of the total ecosystem is
#' observed. In order to draw conclusions about the total ecosystem, it is
#' necessary to adjust for inexhaustive sampling in estimates of alpha
#' diversity. This function implements such an adjustment for the
#' Inverse-Simpson index.
#' 
#' 
#' @param input A frequency count table or vector of abundances.
#' @return An estimate, standard error, and confidence set of the index.
#' @examples 
#' 
#' inverse_simpson_better(apples)
#' 
#' @export inverse_simpson_better
inverse_simpson_better <- function(input) {
  hill_better(input, 2)
}














#' An estimate of the Simpson index
#' 
#' In many ecological surveys, only a small fraction of the total ecosystem is
#' observed. In order to draw conclusions about the total ecosystem, it is
#' necessary to adjust for inexhaustive sampling in estimates of alpha
#' diversity. This function implements such an adjustment for the Simpson
#' index.
#' 
#' 
#' @param input A frequency count table or vector of abundances.
#' @return An estimate, standard error, and confidence set of the index.
#' @export
simpson_better <- function(input, ccc = NA, ccc_se = NA, quantile = 0.95) {
  
  if (length(dim(input)) != 2) {
    stop("Incorrent format: I need a frequency count table")
  } 
  input <- check_format(input)
  
  type <- frequency_count_or_proportion_or_column(input)
  proportions <- to_proportions(input, type)
  
  nn <- sum(input[, 1]*input[, 2])
  
  if ((is.na(ccc) | is.null(ccc))) {
    baway <- breakaway(input,  output = F, answers = T, plot = F)
    ccc <- round(baway$est)
    ccc_se <- round(baway$seest)
  }
  
  cc <- sum(input[,2])
  unobs <- round(ccc-cc)
  
  unobserved_proportions <- rep(1/ccc, unobs)
  observed_proportions <- cc/ccc * proportions
  amended_proportions <- c(unobserved_proportions, observed_proportions)
  
  estimates <- simpson(amended_proportions)
  ses <- sqrt(sum((1 + log(amended_proportions))^2 * amended_proportions * (1-amended_proportions)/nn) )
  
  data.frame("index"  = "shannon", 
             "estimate" = estimates, 
             "standard_error" = ses, 
             "lower" = max(0, estimates-qnorm(1-(1-quantile)/2)*ses), 
             "upper" = estimates+qnorm(1-(1-quantile)/2)*ses)
}
















#' An estimate of the Gini-Simpson index
#' 
#' In many ecological surveys, only a small fraction of the total ecosystem is
#' observed. In order to draw conclusions about the total ecosystem, it is
#' necessary to adjust for inexhaustive sampling in estimates of alpha
#' diversity. This function implements such an adjustment for the Gini-Simpson
#' index.
#' 
#' 
#' @param input A frequency count table or vector of abundances.
#' @return An estimate, standard error, and confidence set of the index.
gini_better <- function(input) {
  stop("You need to email Amy telling her that you want this implemented!")
}


