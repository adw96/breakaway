#' Plot function for alpha_estimates class
#' 
#' @param x Object of class \code{alpha_estimates}.
#' @param method Optional. If there are multiple richness methods included in \code{x}, this can be set to the the desired method to be plotted. Defaults to the method of the first alpha diversity estimate.
#' @param fullCI Optional boolean. Indicates whether you want the plot to include the full confidence intervals. Defaults to \code{FALSE}.
#' @param ... See details
#' 
#' @details ... does not currently have any implemented options. Optional arguments currently include "fullCI", a  Optional
#' 
#' @examples
#' \dontrun{
#' library(phyloseq)
#' data(GlobalPatterns)
#' alphas <- breakaway(GlobalPatterns)
#' plot(alphas)
#' }
#' @export
plot.alpha_estimates <- function(x, method = NULL, fullCI = FALSE, ...) {
  
  if (is.null(method)) {
    all_methods <- x %>% lapply(function(x) x$name) %>% unlist %>% unique
    method <- all_methods[1]
  }

  
  pts_tmp <- x %>% lapply(function(x) x$estimate) %>% unlist %>% na.omit
  if (length(pts_tmp) == 0) {
    stop("There are no estimates in this alpha_estimates object!")
  }
  
  names_tmp <- factor(names(pts_tmp), levels = names(pts_tmp))
  if (length(names_tmp) == 0) {
    names_tmp <- paste("Sample", 1:length(pts_tmp))
  }
  df <- data.frame(pts = pts_tmp,
                   lci = x %>% lapply(function(x) x$interval[1]) %>% unlist,
                   uci = x %>% lapply(function(x) x$interval[2]) %>% unlist,
                   names = names_tmp)
  yname1 <- method
  yname2 <- x[[1]]$estimand
  
  if (!fullCI) {
    fiven <- stats::fivenum(df$uci, na.rm = TRUE)
    iqr <- diff(fiven[c(2, 4)])
    if (!is.na(iqr)) {
      out <- df$uci < (fiven[2L] - 1.5 * iqr) | df$uci > (fiven[4L] + 1.5 * iqr)
      ylower <- min(0, 0.95*min(df$uci[!out]))
      yupper <- 1.05*max(df$uci[!out])
      
      ggplot2::ggplot(df, ggplot2::aes(x = names, xend = names)) +
        ggplot2::geom_point(ggplot2::aes(x = names, y = pts)) +
        ggplot2::coord_cartesian(ylim = c(ylower,yupper)) + 
        ggplot2::geom_segment(ggplot2::aes(y = lci, yend = uci)) +
        ggplot2::ylab(paste(yname1, "estimate of", yname2)) +
        ggplot2::xlab("") +
        ggplot2::theme_bw() + 
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    } else {
      ggplot2::ggplot(df, ggplot2::aes(x = names, xend = names)) +
        ggplot2::geom_point(ggplot2::aes(x = names, y = pts)) +
        ggplot2::geom_segment(ggplot2::aes(y = lci, yend = uci)) +
        ggplot2::ylab(paste(yname1, "estimate of", yname2)) +
        ggplot2::xlab("") +
        ggplot2::theme_bw() + 
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    } # Close last else - happens with no CI or fullCI = TRUE
  } # Close if (!fullCI)
}
