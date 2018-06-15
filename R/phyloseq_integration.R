################################################################################
#' Plot alpha diversity with less incorrect error bars
#'
#' For no good reason, we almost never plot alpha diversity with its appropriate
#' standard errors. This function is intended to fix that.
#'
#'  NOTE: This code is modify from the phyloseq functions plot_richness
#'  and estimate_richness. Please cite phyloseq as well as breakaway
#'  if you use this function!
#' 
#' @param physeq (Required). \code{\link{phyloseq-class}}, or alternatively, 
#'  an \code{\link{otu_table-class}}. The data about which you want to estimate.
#'
#' @param x (Optional). A variable to map to the horizontal axis. The vertical
#'  axis will be mapped to the alpha diversity index/estimate
#'  and have units of total taxa, and/or index value (dimensionless).
#'  This parameter (\code{x}) can be either a character string indicating a
#'  variable in \code{sample_data} 
#'  (among the set returned by \code{sample_variables(physeq)} );
#'  or a custom supplied vector with length equal to the number of samples
#'  in the dataset (nsamples(physeq)).
#'
#'  The default value is \code{"samples"}, which will map each sample's name
#'  to a separate horizontal position in the plot.
#'
#' @param color (Optional). Default \code{NULL}. 
#'  The sample variable to map to different colors.
#'  Like \code{x}, this can be a single character string of the variable name in 
#'  \code{sample_data} 
#'  (among the set returned by \code{sample_variables(physeq)} );
#'  or a custom supplied vector with length equal to the number of samples
#'  in the dataset (nsamples(physeq)).
#'  The color scheme is chosen automatically by \code{link{ggplot}},
#'  but it can be modified afterward with an additional layer using
#'  \code{\link[ggplot2]{scale_color_manual}}.
#'
#' @param shape (Optional). Default \code{NULL}. The sample variable to map
#'  to different shapes. Like \code{x} and \code{color},
#'  this can be a single character string 
#'  of the variable name in 
#'  \code{sample_data} 
#'  (among the set returned by \code{sample_variables(physeq)} );
#'  or a custom supplied vector with length equal to the number of samples
#'  in the dataset (nsamples(physeq)).
#'  The shape scale is chosen automatically by \code{link{ggplot}},
#'  but it can be modified afterward with an additional layer using
#'  \code{\link[ggplot2]{scale_shape_manual}}.
#'
#' @param title (Optional). Default \code{NULL}. Character string.
#'  The main title for the graphic.
#'  
#' @param nrow (Optional). Default is \code{1},
#'  meaning that all plot panels will be placed in a single row,
#'  side-by-side. 
#'  This argument is passed to \code{\link[ggplot2]{facet_wrap}}.
#'  If \code{NULL}, the number of rows and columns will be 
#'  chosen automatically (wrapped) based on the number of panels
#'  and the size of the graphics device.
#'
#' @param measures (Optional). Default is \code{NULL}, meaning that
#'  all available alpha-diversity measures will be included in plot panels.
#'  Alternatively, you can specify one or more measures
#'  as a character vector of measure names.
#'  Values must be among those supported:
#'  \code{c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher")}.
#'
#'
#' @return A \code{\link{ggplot}} plot object summarizing
#'  the richness estimates, and their standard error.
#' 
#' @seealso 
#'  \code{\link{estimate_richness}}
#' 
#' @import reshape2
#' @import phyloseq
#' 
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 geom_errorbar
#' @importFrom ggplot2 facet_wrap
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_text
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 geom_path
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 ggtitle
#' @importFrom ggplot2 ylab
#' 
#' @export
plot_alpha <- function(physeq, x = "samples", color = NULL, shape = NULL,
                       nrow = 1, 
                       measures = NULL,
                       estimators = NULL,
                       upper_limits = NULL) { 
  
  # Calculate the relevant alpha-diversity measures
  estimate_df <- estimate_alpha_better(physeq, split = TRUE, measures = measures, estimators = estimators)
  
  
  # Measures may have been renamed in `estimate_df`. Replace it w with the name from estimate_df
  # Set up non-estimate columns
  measures_ind <- grep("[A-Za-z].*estimate", colnames(estimate_df))
  colnames(estimate_df)[measures_ind] <- unlist(strsplit(colnames(estimate_df)[measures_ind], split = ".estimate"))
  measures <- colnames(estimate_df)[measures_ind]
  
  
  ses <- colnames(estimate_df)[grep("[A-Za-z].*standard_error", colnames(estimate_df))]
  lowers <- colnames(estimate_df)[grep("[A-Za-z].*lower", colnames(estimate_df))]
  uppers <- colnames(estimate_df)[grep("[A-Za-z].*upper", colnames(estimate_df))]
  
  # Make the plotting data.frame.
  # This coerces to data.frame, required for reliable output from reshape2::melt()
  if( !is.null(sample_data(physeq, errorIfNULL=FALSE)) ){
    # Include the sample data, if it is there.
    DF <- data.frame(estimate_df, sample_data(physeq))
  } else {
    # If no sample data, leave it out.
    DF <- data.frame(estimate_df)
  }
  if( !"samples" %in% colnames(DF) ){
    # If there is no "samples" variable in DF, add it
    DF$samples <- sample_names(physeq)
  }
  
  # Make backwards compatible with old formats
  if( !is.null(x) ){
    if( x %in% c("sample", "samples", "sample_names", "sample.names") ){
      x <- "samples"
    }
  } else {
    # If x was NULL for some reason, set it to "samples"
    x <- "samples"
  }
  
  # melt to display different alpha-measures separately
  mdf <- reshape2::melt(DF, measure.vars=measures)
  
  ## Merge s.e. into one "se" column
  # Define conversion vector, `selabs`
  selabs <- ses
  names(selabs) <- strsplit(ses, split = ".standard_error")
  lowerslabs <- lowers
  names(lowerslabs) <- strsplit(lowers, split = ".lower")
  upperslabs <- uppers
  names(upperslabs) <- strsplit(uppers, split = ".upper")
  
  # use selabs conversion vector to process `mdf`
  mdf$wse <- sapply(as.character(mdf$variable), function(i, selabs){selabs[i]}, selabs)
  mdf$wlower <- sapply(as.character(mdf$variable), function(i, lowerslabs){lowerslabs[i]}, lowerslabs)
  mdf$wupper <- sapply(as.character(mdf$variable), function(i, upperslabs){upperslabs[i]}, upperslabs)
  for( i in 1:nrow(mdf) ){
    if( !is.na(mdf[i, "wse"]) ){
      mdf[i, "se"] <- mdf[i, (mdf[i, "wse"])]
      mdf[i, "lower"] <- mdf[i, (mdf[i, "wlower"])]
      mdf[i, "upper"] <- mdf[i, (mdf[i, "wupper"])]
    }
  }
  # prune the redundant columns
  mdf <- mdf[, -which(colnames(mdf) %in% c(selabs, "wse", lowerslabs, "wlower", upperslabs, "wupper"))]
  
  # Allow truncation of axes
  # 
  # Ideally, want to do this in the plot, facet_wrap makes this really challenging.
  # If anyone knows how to do this in the facet_wrap, please let me know!
  if (any(!is.null(upper_limits))) {
    for (measure in measures) {
      index <- which(measure == measures)
      if (!(is.na(upper_limits[index]) | is.null(upper_limits[index]))) {
        mdf$value[which(mdf$variable == measure & mdf$value >  upper_limits[index])] <- NA
        mdf$upper[which(mdf$variable == measure & mdf$upper >  upper_limits[index])] <- upper_limits[index]
      }
    }
  }

  
  # Define variable mapping
  richness_map <- aes_string(x=x, y="value", colour=color, shape=shape)
  
  # Make the ggplot.
  p <- ggplot(mdf, richness_map) + geom_point(na.rm=TRUE)  
  
  # Add error bars
  p <- p + geom_errorbar(aes(ymin=lower, ymax=upper), width=0.1)
  
  # Rotate horizontal axis labels, and adjust
  p <- p + theme(axis.text.x=element_text(angle=-90, vjust=0.5, hjust=0))
  
  # Add y-label
  p <- p + ylab('Alpha Diversity Measure')
  
  # Facet wrap using user-options
  p <- p + facet_wrap(~variable, nrow=nrow, scales="free_y")
  
  ## Long-term goal: allow user-defined limits
  # if (!is.null(upper_limits)) {
  #   q <- ggplot_build(p)
  #   for (i in which(!(is.na(upper_limits) | is.null(upper_limits)))) {
  #     q$layout$panel_scales$y[[i]]$range$range[2] <- upper_limits[i]
  #   }
  #   
  #   grid.newpage()
  #   grid.draw(ggplot_gtable(q))
  # }
  # 
  
  return(p)
}



################################################################################
#' Estimate alpha diversity
#'
#' Estimates a number of standard alpha diversity in a statistically rigourous way, 
#' and returns the results as a \code{data.frame}.
#' Strictly speaking, this function is not only estimating richness,
#' despite its name. 
#' It can operate on the cumulative population of all
#' samples in the dataset, or by repeating the richness estimates for each
#' sample individually.
#' NOTE: You must use untrimmed datasets
#' for meaningful results, as these estimates (and even the ``observed'' richness)
#' are highly dependent on the number of singletons. You can always trim the data
#' later on if needed, just not before using this function.
#' 
#' @param physeq (Required). \code{\link{phyloseq-class}}, or alternatively, 
#'  an \code{\link{otu_table-class}}. The data about which you want to estimate
#'  the richness.
#'
#' @param split (Optional). Logical. Should a separate set of richness estimates
#'  be performed for each sample? Or alternatively, pool all samples and 
#'  estimate richness of the entire set.
#'  
#' @param measures (Optional). Default is \code{NULL}, meaning that
#'  all available alpha-diversity measures will be included.
#'  Alternatively, you can specify one or more measures
#'  as a character vector of measure names.
#'  Values must be among those supported:
#'  \code{c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher")}.
#'
#' @return A \code{data.frame} of the richness estimates, and their standard error.
#' 
#' @import phyloseq
#' 
#' @seealso 
#'  Check out the custom plotting function, \code{\link{plot_richness}},
#'  for easily showing the results of different estimates, 
#'  with method-specific error-bars.
#'  Also check out the internal functions borrowed from the \code{vegan} package:
#'  
#' @export estimate_alpha_better
estimate_alpha_better <- function(physeq, split=TRUE, measures=NULL, estimators = NULL){
  
  if( !any(otu_table(physeq)==1) ){
    # Check for singletons, and then warning if they are missing.
    # These metrics only really meaningful if singletons are included.
    warning(
      "Your data doesn't have any singletons. \n",
      "This is generally makes alpha diversity estimation extremely difficult.\n",
      "It is recommended that you find the un-trimmed data and retry.\n", 
      "Otherwise, proceed with extreme caution."
    )
  }
  
  # If we are not splitting sample-wise, sum the species. Else, enforce orientation.
  if( !split ){
    OTU <- taxa_sums(physeq)		
  } else if( split ){
    OTU <- as(otu_table(physeq), "matrix")
    if( taxa_are_rows(physeq) ){ OTU <- t(OTU) }
  }
  
  # Define renaming vector:
  renamevec <- c("richness", "shannon", "simpson", "invsimpson")
  names(renamevec) <- renamevec
  
  # If measures was not explicitly provided (is NULL), set to all supported methods
  if( is.null(measures) ){
    measures <- renamevec
  }
  
  # Stop with error if no measures are supported
  if( !any(measures %in% renamevec) ){
    stop("None of the `measures` you provided are supported.\n Try default `NULL` instead.")
  }
  
  # If estimators was not explicitly provided, set to the following defaults
  if( is.null(estimators) ){
    estimators <- list("richness" = richness_clean, 
                       "shannon" = shannon_better,
                       "simpson" = simpson_better, 
                       "invsimpson" = inverse_simpson_better)
  }
  
  # Initialize to NULL
  outlist <- vector("list")
  
  freq_counts <- apply(OTU, 1, make_frequency_count_table)
  freq_props <- apply(OTU, 2, to_proportions, type = "column")
  
  if ("richness" %in% measures) {
    richness_ests <- lapply(freq_counts, FUN = estimators$richness) 
    outlist <- c(outlist, list(richness = richness_ests))
  }
  if( "shannon" %in% measures ){
    shannon_ests <- lapply(freq_counts, FUN = estimators$shannon)
    outlist <- c(outlist, list(shannon = shannon_ests))
  }
  if( "simpson" %in% measures ){
    simpson_ests <-  lapply(freq_counts, FUN = estimators$simpson)
    outlist <- c(outlist, list(simpson = simpson_ests))
  }
  if( "invsimpson" %in% measures ){
    invsimpson_ests <-  lapply(freq_counts, FUN = estimators$invsimpson)
    outlist <- c(outlist, list(invsimpson = invsimpson_ests))
  }
  
  out_list <- do.call("cbind", outlist)
  
  out <- t(apply(out_list, 1, unlist))
  
  redundant <- sapply(FUN = grepl, X = colnames(out), pattern = "index", fixed = TRUE)
  out <- out[, !redundant, drop=FALSE]
  
  # Make sure that you return a data.frame for reliable performance.
  out <- as.data.frame(out)
  return(out)
}

#' @export test_alpha
test_alpha <- function(physeq, x, measure = "richness") {
  estimate_df <- estimate_alpha_better(physeq, measures = measure)
  
  measures_ind <- grep("[A-Za-z].*estimate", colnames(estimate_df))
  colnames(estimate_df)[measures_ind] <- unlist(strsplit(colnames(estimate_df)[measures_ind], split = ".estimate"))
  measures <- colnames(estimate_df)[measures_ind]
  
  ses <- colnames(estimate_df)[grep("[A-Za-z].*standard_error", colnames(estimate_df))]
  
  # Make the plotting data.frame.
  # This coerces to data.frame, required for reliable output from reshape2::melt()
  if( !is.null(sample_data(physeq, errorIfNULL=FALSE)) ){
    # Include the sample data, if it is there.
    DF <- data.frame(estimate_df, sample_data(physeq))
  } else {
    stop("No sample data in your phyloseq object!")
  }
  if( !"samples" %in% colnames(DF) ){
    # If there is no "samples" variable in DF, add it
    DF$samples <- sample_names(physeq)
  }
  
  # melt to display different alpha-measures separately
  mdf <- reshape2::melt(DF, measure.vars=measures)
  
  ## Merge s.e. into one "se" column
  # Define conversion vector, `selabs`
  selabs <- ses
  names(selabs) <- strsplit(ses, split = ".standard_error")
  
  # use selabs conversion vector to process `mdf`
  mdf$wse <- sapply(as.character(mdf$variable), function(i, selabs){selabs[i]}, selabs)
  for( i in 1:nrow(mdf) ){
    if( !is.na(mdf[i, "wse"]) ){
      mdf[i, "se"] <- mdf[i, (mdf[i, "wse"])]
    }
  }
  
  X <- model.matrix(lm(value ~ mdf[[x]], data = mdf))
  colnames(X) <- c("Intercept", levels(mdf[[x]])[-1])
  betta(mdf$value, mdf$se, X)$table
  
}