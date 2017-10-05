#' Resampling-based estimates of alpha diversity
#' 
#' This function resamples OTUs to a specified read count, then calculates a
#' given alpha diversity metric based on the resample. Useful for determining
#' variance of alpha diversity estimates.
#' 
#' 
#' @param my_data The columns of your otu table, giving the number of counts
#' for each taxa
#' @param my_function The alpha diversity metric that you're interested in
#' (e.g. Shannon)
#' @param my_sample_size The number of reads that you want to resample. Common
#' estimates of alpha diversity (such as the plug-in estimate of Shannon
#' diversity) are very sensitive to sample size. To gauge the effect of
#' differing sample sizes on the variance in estimating alpha diversity, it is
#' (currently) advisable to include multiple different sample sizes in your
#' resampling
#' @param ... Other arguments to pass to \code{my_function}
#' @return A (independent-bootstrap resampled) estimate of the alpha diversity
#' metric given by my_function
#' @author Amy Willis
#' @examples
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' # Bootstrap observations to calculate a re-sample estimate of the Shannon index
#' resample_estimate(toy_otu_table[,1], shannon)
#' 
#' # To see the sampling distribution of the resample estimates, repeat this 200 times
#' shannon_estimates <- replicate(200, resample_estimate(toy_otu_table[,1], shannon))
#' hist(shannon_estimates)
#' 
#' # To investigate the effect of sample size on the resample estimate of the Shannon index,
#' # re-sample over multiple sample sizes; in this case, the number of reads
#' ns <- apply(toy_otu_table, 2, sum)
#' shannon_estimates <- replicate(200, resample_estimate(toy_otu_table[,1], shannon, my_sample_size = ns))
#' 
#' hist(shannon_estimates)
#' sd(shannon_estimates) # a not too terrible estimate of the standard error in estimating Shannon
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
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




















#' Resampling the taxa of your OTU table
#' 
#' This function resamples OTUs to a specified read count
#' 
#' 
#' @param my_data The columns of your OTU table, giving the number of counts
#' for each taxa
#' @param my_sample_size The number of reads of the bootstrap.
#' @param ... Other arguments, which are ignored.
#' @return A independent-bootstrap resampled OTU frequencies, in frequency
#' count table form.
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
