#' plot_alpha_better <- function(physeq, estimators, x="samples", 
#'                               color=NULL, shape=NULL, title=NULL,
#'                               scales="free_y", nrow=1, sortby=NULL) {
#'   
#'   # estimate alpha
#'   erDF = estimate_alpha_better(physeq, estimators=estimators)
#'   
#'   # estimators may have been renamed in `erDF`. Replace it with the name from erDF
#'   estimators = colnames(erDF)
#'   # Define "measure" variables and s.e. labels, for melting.
#'   ses = colnames(erDF)[grep("^se\\.", colnames(erDF))]
#'   # Remove any S.E. from `estimators`
#'   estimators = estimators[!estimators %in% ses]
#'   # Make the plotting data.frame.
#'   # This coerces to data.frame, required for reliable output from reshape2::melt()
#'   if( !is.null(sample_data(physeq, errorIfNULL=FALSE)) ){
#'     # Include the sample data, if it is there.
#'     DF <- data.frame(erDF, sample_data(physeq))
#'   } else {
#'     # If no sample data, leave it out.
#'     DF <- data.frame(erDF)
#'   }
#'   if( !"samples" %in% colnames(DF) ){
#'     # If there is no "samples" variable in DF, add it
#'     DF$samples <- sample_names(physeq)
#'   }
#'   # sample_names used to be default, and should also work.
#'   # #backwardcompatibility
#'   if( !is.null(x) ){
#'     if( x %in% c("sample", "samples", "sample_names", "sample.names") ){
#'       x <- "samples"
#'     }
#'   } else {
#'     # If x was NULL for some reason, set it to "samples"
#'     x <- "samples"
#'   }
#'   # melt to display different alpha-estimators separately
#'   mdf = reshape2::melt(DF, measure.vars=estimators)
#'   # Initialize the se column. Helpful even if not used.
#'   mdf$se <- NA_integer_
#'   if( length(ses) > 0 ){
#'     ## Merge s.e. into one "se" column
#'     # Define conversion vector, `selabs`
#'     selabs = ses
#'     # Trim the "se." from the names
#'     names(selabs) <- substr(selabs, 4, 100)
#'     # Make first letter of selabs' names uppercase
#'     substr(names(selabs), 1, 1) <- toupper(substr(names(selabs), 1, 1))
#'     # use selabs conversion vector to process `mdf`
#'     mdf$wse <- sapply(as.character(mdf$variable), function(i, selabs){selabs[i]}, selabs)
#'     for( i in 1:nrow(mdf) ){
#'       if( !is.na(mdf[i, "wse"]) ){
#'         mdf[i, "se"] <- mdf[i, (mdf[i, "wse"])]
#'       }
#'     }
#'     # prune the redundant columns
#'     mdf <- mdf[, -which(colnames(mdf) %in% c(selabs, "wse"))]
#'   }
#'   ## Interpret estimators
#'   # If not provided (default), keep all 
#'   if( !is.null(estimators) ){
#'     if( any(estimators %in% as.character(mdf$variable)) ){
#'       # If any estimators were in mdf, then subset to just those.
#'       mdf <- mdf[as.character(mdf$variable) %in% estimators, ]
#'     } else {
#'       # Else, print warning about bad option choice for estimators, keeping all.
#'       warning("Argument to `estimators` not supported. All alpha-diversity estimators (should be) included in plot.")
#'     }
#'   }
#'   # Address `sortby` argument
#'   if(!is.null(sortby)){
#'     if(!all(sortby %in% levels(mdf$variable))){
#'       warning("`sortby` argument not among `estimators`. Ignored.")
#'     }
#'     # if(!is.discrete(mdf[, x])){
#'     #   warning("`sortby` argument provided, but `x` not a discrete variable. `sortby` is ignored.")
#'     # }
#'     if(all(sortby %in% levels(mdf$variable)) & is.discrete(mdf[, x])){
#'       # Replace x-factor with same factor that has levels re-ordered according to `sortby`
#'       wh.sortby = which(mdf$variable %in% sortby)
#'       mdf[, x] <- factor(mdf[, x],
#'                          levels = names(sort(tapply(X = mdf[wh.sortby, "value"],
#'                                                     INDEX = mdf[wh.sortby, x],
#'                                                     mean,
#'                                                     na.rm=TRUE, simplify = TRUE))))
#'     }
#'   }
#'   # Define variable mapping
#'   richness_map = aes_string(x=x, y="value", colour=color, shape=shape)
#'   # Make the ggplot.
#'   p = ggplot(mdf, richness_map) + geom_point(na.rm=TRUE)  
#'   # Add error bars if mdf$se is not all NA
#'   if( any(!is.na(mdf[, "se"])) ){
#'     p = p + geom_errorbar(aes(ymax=value + se, ymin=value - se), width=0.1) 
#'   }
#'   # Rotate horizontal axis labels, and adjust
#'   p = p + theme(axis.text.x=element_text(angle=-90, vjust=0.5, hjust=0))
#'   # Add y-label 
#'   p = p + ylab('Alpha Diversity Measure') 
#'   # Facet wrap using user-options
#'   p = p + facet_wrap(~variable, nrow=nrow, scales=scales)
#'   # Optionally add a title to the plot
#'   if( !is.null(title) ){
#'     p <- p + ggtitle(title)
#'   }
#'   return(p)
#'   
#' }

format_breakaway <- function(counts) {
  frequency_table <- build_frequency_count_tables(counts)
  baway <- try(breakaway(frequency_table))
  if (class(baway) != "try-error") {
    c("estimate" = baway$est, 
      "stderror" = baway$seest, 
      "lcb" = baway$ci[1], 
      "ucb" = baway$ci[2])
  } else {
    c("estimate" = NA,
      "stderror" = NA, 
      "lcb" = sum(frequency_table[,2]), 
      "ucb" = Inf)
  }
}

#' @export
estimate_alpha_better <- function(physeq, estimators) {
  
  if( !any(otu_table(physeq)==1) ){
    warning(
      "The data you have provided does not have\n",
      "any singletons. Quality control is important,\n",
      "... TODO..."
    )
  }
  
  # OTU is a matrix with taxa in the columns
  OTU <- as(otu_table(physeq), "matrix")
  if( taxa_are_rows(physeq) ){ OTU <- t(OTU)  }
  
  # Initialize to NULL
  outlist = vector("list")
  if( "breakaway" %in% estimators ){
    
    baway_matrix <- apply(OTU, 1, format_breakaway)
    
    outlist <- c(outlist, list(breakaway = baway_matrix))
  } 
  
  
  if ("objectiveBayes" %in% estimators) {
    warning("OB in progress")
  }
  
  if ("shannon_better" %in% estimators) {
    warning("shannon_better in progress")
  }
  
  
  out = do.call("cbind", outlist)
  # # Rename columns per renamevec
  # namechange = intersect(colnames(out), names(renamevec))
  # colnames(out)[colnames(out) %in% namechange] <- renamevec[namechange]
  # Final prune to just those columns related to "estimators". Use grep.
  colkeep = sapply(paste0("(se\\.){0,}", estimators), grep, colnames(out), ignore.case=TRUE)
  out = out[, sort(unique(unlist(colkeep))), drop=FALSE]
  # Make sure that you return a data.frame for reliable performance.
  out <- as.data.frame(out)
  return(out)
}

# for (file in list.files("/Users/amy/Documents/software/breakaway/R", full.names = T)) source(file)
# plot_alpha_better(GlobalPatterns, estimators = c("breakaway", "shannon_better"))

