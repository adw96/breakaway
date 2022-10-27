#' Build frequency count tables from an OTU table
#' 
#' Build frequency count tables from an OTU table
#' 
#' 
#' @param the_table An OTU table as a data frame or a matrix. Columns are the
#' samples and rows give the taxa.
#' @return A list of frequency count tables corresponding to the columns.
#' @export build_frequency_count_tables
build_frequency_count_tables <- function(the_table) {
  
  if (intersect(c("phyloseq", "otu_table"),
                class(the_table)) %>% length > 0) {
    
    if (the_table %>% otu_table %>% taxa_are_rows) {
      otus <- the_table %>% otu_table %>% get_taxa %>% t
    } else {
      otus <- the_table %>%  otu_table
    }
    
    tbls <- otus %>%
      apply(1, make_frequency_count_table)
    
  } else {
    
    if (dim(the_table)[1] < dim(the_table)[2]) {
      warning('More columns then rows. You probably need to transpose your data.')
    }
    
    tbls <- apply(the_table, 2, convert)
    names(tbls) <- colnames(the_table)
  }
  
  tbls
}




#' Draw frequency count subtables from an OTU table
#' 
#' Draw frequency count subtables from an OTU table
#' 
#' @param labels A vector of counts of the taxa; i.e. a vector giving the number of times each taxon was observed.
#' @return A frequency count table.
#'
#' @export make_frequency_count_table
make_frequency_count_table <- function(labels) {
  x <- as.data.frame(table(labels[labels != 0]))
  x[, 1] <- as.numeric(as.character(x[, 1]))
  x[, 2] <- as.numeric(as.character(x[, 2]))
  x
}

#' OTU table to relative abundances
#' 
#' @param the_table An OTU table
#' 
#' @return A proportion table or vector.
#' 
#' @export
proportions_instead <- function(the_table) {
  if (!(is.null(dim(the_table)) | is.vector(the_table))) {
    #  is in frequency count form; convert to long form
    long_form <- rep(the_table[, 1], times = the_table[, 2])
    long_form/sum(long_form)
  } else {
    ## must be a list of counts
    the_table/sum(the_table)
  }
}

