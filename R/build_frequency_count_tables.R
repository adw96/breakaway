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
#' 
#' @export make_frequency_count_table
make_frequency_count_table <- function(labels) {
  x <- as.data.frame(table(table(labels)))
  x[, 1] <- as.numeric(as.character(x[, 1]))
  x
}

convert_freq_indices <- function(frequency_table) {
  checked_frequency_table <- frequency_table
  if (is.factor(frequency_table[, 1])) {
    fs <- as.numeric(as.character(frequency_table[, 1]))
    checked_frequency_table <- cbind("Index"=fs, "Frequency" = frequency_table[, 2])
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
