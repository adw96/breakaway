#' @export
hill <- function(data, q) {
  if (class(data) %in% c("data.frame", "matrix") ) {
    if (ncol(data) == 2) {
      # convert to long form
      data <- rep(data[, 1], times = data[, 2])
    }
  }
  data <- data[data>0]
  if (max(data) > 1) {
    # normalize if not proportions
    data <- data/sum(data)
  }
  if (q == 1) {
    shannon(data)
  } else {
    (sum(data^q))^(1/(1-q))
  }
  
}

#' @export
shannon <- function(data) {
  if (class(data) %in% c("data.frame", "matrix") ) {
    if (ncol(data) == 2) {
      data <- rep(data[, 1], times = data[, 2])
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
      data <- rep(data[, 1], times = data[, 2])
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
      data <- rep(data[, 1], times = data[, 2])
    }
  }
  if (max(data) > 1) {
    # normalize if not proportions
    data <- data/sum(data)
  }
  data <- data[data>0]
  1/sum(data^2)
}
