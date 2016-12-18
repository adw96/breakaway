build_frequency_count_tables <- function(otu_table) {
  if (dim(otu_table)[1] < dim(otu_table)[2]) warning('More columns then rows. You probably need to transpose your data.')
  
  frequency_table_list_with_zeros <- lapply(apply(otu_table, 2, table), as.data.frame)
  
  frequency_table_list_without_zeros <- lapply(frequency_table_list_with_zeros, function(x) x[x[,1]!=0,])
  
  ## correct table making indices categorical
  frequencytablelist <- lapply(frequency_table_list_without_zeros, .convert_freq_indices)
  frequencytablelist
}

.convert_freq_indices <- function(frequency_table) {
  checked_frequency_table <- frequency_table
  if (is.factor(frequency_table[, 1])) {
    fs <- as.numeric(as.character(frequency_table[, 1]))
    checked_frequency_table <- cbind("Index"=fs, "Frequency" = frequency_table[, 2])
  }
  checked_frequency_table
}
