#' convert between different inputs for alpha-diversity estimates
#' 
#' Inputs slated for development include phyloseq and otu_table
#' 
#' @param input_data Supported types include filenames, data frames, matrices, vectors...
#' 
#' @return Frequency count able 
#' 
#' @export
convert <- function(input_data) {
  
  if ("phyloseq" %in% class(input_data)) {
    
    stop("phyloseq not yet supported as an input type to function `convert`")
    
  } else if ("otu_table" %in% class(input_data)) {
    
    stop("otu_table not yet supported as an input type to function `convert`")
    
  } else if (class(input_data) == "character") {
    
    filename <- input_data
    ext <- substr(filename, nchar(filename) - 2, nchar(filename))
    
    if (ext %in% c("csv", "txt")) {
      
      output_data <- read.table(file=filename, header=0, sep=",")
      # TODO figure out if there is a header

    } else {
      
      stop("Input is a string, but not a valid file path.")
      
    }
    
  } else if ("data.frame" %in% class(input_data) | "matrix" %in% class(input_data)) {
    
    # determine if frequency count table or vector
    if (dim(input_data)[2] != 2) {
      stop("input_data is a data.frame or matrix but not a frequency count table.\n")
    } 
    
    
    output_data <- input_data
    
  } else if (class(input_data) == "numeric") {
    
    # must be vector of counts
    if (isTRUE(all.equal(sum(input_data), 1))) {
      stop("species richness estimates cannot accept relative abundances")
    }
    
    if (any(input_data %% 1 != 0)) {
      stop("input_data is a vector but not a vector of integers")
    }
    
    output_data <- as.data.frame(table(input_data))
    
  } else {
    
    stop("Unsupported input type to function `convert`")
    
  }
  
  checked_output_data <- check_format(output_data)
  checked_output_data
  
}

#' Run some basic checks on a possible frequency count table
#' 
#' @param output_data A matrix to test
#' 
#' @return A checked frequency count table
check_format <- function(output_data) {
  
  if (length(output_data) <= 1) {
    stop("Input data to `check_format` is of length 1 or 0.")
  }
  
  if (is.factor(output_data[, 1]) ) {
    
    ## if table is used to create the frequency tables, the frequency index column is usually a factor, so fix this here
    
    fs <- as.numeric(as.character(output_data[, 1]))
    output_data <- cbind(fs, output_data[, 2])
    output_data <- output_data[output_data[, 1] != 0, ]
    
  }
  
  # remove rows with zero
  output_data <- output_data[output_data[, 2] != 0, ]
  
  colnames(output_data) <- c("index", "frequency")
  
  output_data
}