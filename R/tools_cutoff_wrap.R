#' Find a cut-off for estimates relying on contiguous counts
#' 
#' @param my_data Frequency count table
#' @param requested The user-requested cutoff
#' 
#' @return Cutoff value
cutoff_wrap <- function(my_data, requested = NA) {
  
  iss <- my_data$index
  fis <- my_data$frequency
  length_fis <- length(fis)
  
  breaks <- which(iss[-1] - iss[-length_fis] > 1)
  cutoff <- ifelse(is.na(breaks[1]), length_fis, breaks[1])
  
  if (!is.na(requested)) {
    # if the requested cutoff is lte cutoff, use cutoff, o/w, cutoff
    if (requested <= cutoff) {
      cutoff <- requested # ok
    } else {
      warning("ignoring requested cutoff; it's too low")
      cutoff <- cutoff
    }
  }
  
  return(cutoff)
}