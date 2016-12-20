#' @export
shannon <- function(data) {
  if (class(data) %in% c("data.frame", "matrix") ) {
    if (ncol(data) == 2) {
      data <- data[,2]
    }
  }
  if (max(data) > 1) {
    # normalize if not proportions
    data <- data/sum(data)
  }
  data <- data[data>0]
  -sum(data*log2(data))
}

#' @export
shannon_e  <- function(data) {
  if (class(data) %in% c("data.frame", "matrix") ) {
    if (ncol(data) == 2) {
      data <- data[,2]
    }
  }
  if (max(data) > 1) {
    # normalize if not proportions
    data <- data/sum(data)
  }
  data <- data[data>0]
  -sum(data*log2(data))/log2(length(data)) 
  ## is this the best way to do it? I'm not sure
  ## really need to scale up with the unseen data
  ## i.e. number of observed species is the worst possible estimate for total species
  ## please be very careful;  this function is in development
}

#' @export
simpson <- function(data) {
  if (class(data) %in% c("data.frame", "matrix") ) {
    if (ncol(data) == 2) {
      data <- data[,2]
    }
  }
  if (max(data) > 1) {
    # normalize if not proportions
    data <- data/sum(data)
  }
  data <- data[data>0]
  sum(data^2)
}
