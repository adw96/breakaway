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
shannon <- function(input) {
  
  type <- frequency_count_or_proportion_or_column(input)
  
  if (type == "frequency count") {
    total_reads <- sum(input[, 2]*input[,1])
    relative_abundance <- input[, 1] / total_reads
    -sum(input[, 2] * relative_abundance*log(relative_abundance))
  } else if (type == "proportion") {
    input_nonzero <- input[input>0]
    -sum(input_nonzero*log(input_nonzero))
  } else if (type == "column") {
    input <- input/sum(input)
    input_nonzero <- input[input>0]
    -sum(input_nonzero*log(input_nonzero))
  } else {
    error("Invalid input")
  }
}

#' @export
hill <- function(input, q) {
  type <- frequency_count_or_proportion_or_column(input)
  
  if (type == "frequency count") {
    total_reads <- sum(input[, 2]*input[,1])
    relative_abundance <- input[, 1] / total_reads
    sum(input[, 2] * relative_abundance^q)^(1/(1-q))
  } else if (type == "proportion") {
    input_nonzero <- input[input>0]
    (sum(input_nonzero^q))^(1/(1-q))
  } else if (type == "column") {
    input <- input/sum(input)
    input_nonzero <- input[input>0]
    (sum(input_nonzero^q))^(1/(1-q))
  } else {
    error("Invalid input")
  }
}

#' @export
inverse_simpson <- function(input) {
  hill(input, 2)
}

#' @export
simpson <- function(input) {
  1/hill(input, 2)
}

#' @export
gini <- function(input) {
  1-simpson(input)
}
