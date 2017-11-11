#' Build frequency count tables from an OTU table
#' 
#' Build frequency count tables from an OTU table
#' 
#' 
#' @param otu_table An OTU table as a data frame or a matrix. Columns are the
#' samples and rows give the taxa.
#' @return A list of frequency count tables corresponding to the columns.
#' @export build_frequency_count_tables
build_frequency_count_tables <- function(otu_table) {
  if (dim(otu_table)[1] < dim(otu_table)[2]) warning('More columns then rows. You probably need to transpose your data.')
  
  frequency_table_list_with_zeros <- lapply(apply(otu_table, 2, table), as.data.frame)
  
  frequency_table_list_without_zeros <- lapply(frequency_table_list_with_zeros, function(x) x[x[,1]!=0,])
  
  ## correct table making indices categorical
  frequencytablelist <- lapply(frequency_table_list_without_zeros, convert_freq_indices)
  frequencytablelist
}











































#' Draw frequency count subtables from an OTU table
#' 
#' Draw frequency count subtables from an OTU table
#' 
#' @param labels A vector of counts of the taxa; i.e. a vector giving the number of times each taxon was observed.
#' 
#' @export make_frequency_count_table
make_frequency_count_table <- function(labels) {
  x <- as.data.frame(table(labels[labels != 0]))
  x[, 1] <- as.numeric(as.character(x[, 1]))
  x[, 2] <- as.numeric(as.character(x[, 2]))
  x
}

convert_freq_indices <- function(frequency_table) {
  checked_frequency_table <- frequency_table
  if (!is.numeric(frequency_table[, 1])) {
    fs <- as.numeric(as.character(frequency_table[, 1]))
    checked_frequency_table <- cbind("Index"=fs, "Frequency" = as.numeric(as.character(frequency_table[, 2])))
  }
  checked_frequency_table
}

#' @export
proportions_instead <- function(otu_table) {
  if (!(is.null(dim(otu_table)) | is.vector(otu_table))) {
    #  is in frequency count form; convert to long form
    long_form <- rep(otu_table[, 1], times = otu_table[, 2])
    long_form/sum(long_form)
  } else {
    ## must be a list of counts
    otu_table/sum(otu_table)
  }
}

check_format <- function(my_data) {
  
  ## read in data
  if( !(is.matrix(my_data) || is.data.frame(my_data))) {
    print(head(my_data))
    stop(paste("Input data of class ", class(my_data), ". ",
               "Please input your data as a data frame or matrix.", sep = ""))
  }
  
  if(length(my_data) <= 1) {
    stop("Input data is of length 1 or 0. Huh?")
  }
  

  
  ## if tables() is used to create the frequency tables, the frequency index column is usually a factor, so fix this here
  if ( is.factor(my_data[,1]) ) {
    fs <- as.numeric(as.character(my_data[,1]))
    my_data <- cbind(fs,my_data[,2])
    my_data <- my_data[my_data[,1]!=0,]
  }
  
  my_data[!(my_data[,2]==0 | is.na(my_data[,2])),]
  
}