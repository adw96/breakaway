#' Wrapper for phyloseq objects
#' 
#' @param data Frequency count table
#' 
#' @return Cutoff value
#' @export
cutoff_wrap <- function(data) {
  iss <- data$index
  fis <- data$frequency
  length_fis <- length(fis)
  
  breaks <- which(iss[-1] - iss[-length_fis] > 1)
  cutoff <- ifelse(is.na(breaks[1]), length_fis, breaks[1])
  return(cutoff)
}