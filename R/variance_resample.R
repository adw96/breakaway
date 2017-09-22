# a function to resample otu vector and recalculate the value
# of an estimate
#' @export resample_estimate
resample_estimate <- function(my_data, my_function, my_sample_size=NA, ...) {
  # assume  my_data is a list of frequencies from an otu table
  if (max(my_data) < 1) {
    warning("Are you sure these are OTU counts?")
  }
  if (is.na(my_sample_size)[1]) my_sample_size <- sum(my_data)
  my_data2 <- my_data[my_data > 0]
  if (length(my_data2) > 0) {
    category_names <- paste("sample", 1:length(my_data2), sep = "")
    my_data_long <- rep(category_names, times = my_data2)
    my_resample <- sample(my_data_long, size = sample(my_sample_size), replace = TRUE)
    my_function(table(my_resample), ...) 
  } else {
    my_function(my_data)
  }
}

#' @export subsample_otu
subsample_otu <- function(my_data, my_sample_size=NA, ...) {
  # assume  my_data is an OTU table
  if (is.na(my_sample_size)[1]) {
    my_sample_size <- sum(my_data[, 1] * my_data[,2])
    my_sample_size <- round(runif(n = 1, min = 0.6*my_sample_size, max = my_sample_size))
  }
  
  num_taxa <- sum(my_data[, 1] * my_data[,2])
  num_species <- sum(my_data[,2])
  abundances <- rep(my_data[,1], times = my_data[,2])
  relative_abundances <- abundances / sum(abundances)
  
  # new ID labels
  new_sample <- sample(x = 1:num_species, size = my_sample_size, replace = T, prob = relative_abundances)
  make_frequency_count_table(new_sample)
  
}
