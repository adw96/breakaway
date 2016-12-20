# a function to resample otu vector and recalculate the value
# of an estimate
#' @export
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