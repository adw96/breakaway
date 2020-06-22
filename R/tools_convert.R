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
  
  if (all(class(input_data) == "character")) {
    
    stop("breakaway no longer supports file paths as inputs")
    
  } else if ("data.frame" %in% class(input_data) | 
             "matrix" %in% class(input_data)) {

    # determine if frequency count table or vector
    if (dim(input_data)[2] != 2) {
      stop("input_data is a data.frame or matrix but not a frequency count table.\n")
    }
    
    if ("tbl_df" %in% class(input_data)) {
      input_data <- input_data %>% as.data.frame
    }

    output_data <- input_data
    
  } else if (class(input_data) %in% c("numeric", "integer")) {
    
    # must be vector of counts
    if (isTRUE(all.equal(sum(input_data), 1)) &
        length(unique(input_data)) > 2) {
      stop("species richness estimates cannot accept relative abundances")
    }
    
    if (any(input_data %% 1 != 0)) {
      stop("input_data is a vector but not a vector of integers")
    }
    
    input_data_remove_zeros <- input_data[input_data != 0]
    output_data <- as.data.frame(table(input_data_remove_zeros))
    
  } else {
    stop(paste("Unsupported input type to function `convert`.",
               "You passed in an object of class", class(input_data)))
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
    warning("Input data to `check_format` is of length 1 or 0.")
    return(NULL)
  }
  
  if(length(intersect(class(output_data), c("matrix", "data.frame"))) == 0) stop("Input should be a matrix or a data frame")
  
  if(length(dim(output_data)) != 2) stop("Input should have 2 columns")
  
  if(any(output_data[,2] %% 1 != 0)) stop("Second input column not integer-valued; should be counts")
  
  if(!all(rank(output_data[,1]) == 1:length(output_data[,1]))) warning("Frequency count format, right?")
  
  
  ## if table is used to create the frequency tables, the frequency index column is usually a factor, so fix this here
  if (output_data[,1] %>% class == "factor") {
    output_data[,1] %<>% as.character %>% as.integer
  }
  
  # remove rows with zero
  output_data <- output_data[output_data[, 2] != 0, ]
  output_data %<>% data.frame
  
  colnames(output_data) <- c("index", "frequency")
  
  if (!all(sort(output_data[, 1]) == output_data[, 1])) {
    stop("frequency counts not in order in `convert`")
  }
  
  output_data
}