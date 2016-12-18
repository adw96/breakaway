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
