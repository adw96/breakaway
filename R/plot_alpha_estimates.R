#' Plot function for alpha_estimates class
#' 
#' @param x Object of class \code{alpha_estimates}
#' @param ... See details
#' 
#' @details ... Can currently include "fullCI", a boolean to indicate whether you want the plot to include the full confidence intervals. Defaults to FALSE.
#' 
#' @examples
#' \dontrun{
#' library(phyloseq)
#' data(GlobalPatterns)
#' alphas <- breakaway(GlobalPatterns)
#' plot(alphas)
#' }
#' @export
plot.alpha_estimates <- function(x, ...) {
  
  args <- match.call(expand.dots = TRUE)
  if (is.null(args$fullCI)) {
    args$fullCI <- FALSE
  }
  
  pts_tmp <- x %>% lapply(function(x) x$estimate) %>% unlist %>% na.omit
  df <- data.frame(pts = pts_tmp,
                   lci = x %>% lapply(function(x) x$interval[1]) %>% unlist,
                   uci = x %>% lapply(function(x) x$interval[2]) %>% unlist,
                   names = factor(names(pts_tmp), levels = names(pts_tmp)))
  yname1 <- x[[1]]$name
  yname2 <- x[[1]]$estimand
  xx <- names(pts_tmp)
  
  if (args$fullCI) {
    ggplot2::ggplot(df, ggplot2::aes(x = names, xend = names)) +
      ggplot2::geom_point(ggplot2::aes(x = names, y = pts)) +
      ggplot2::geom_segment(ggplot2::aes(y = lci, yend = uci)) +
      ggplot2::ylab(paste(yname1, "estimate of", yname2)) +
      ggplot2::xlab("") +
      ggplot2::theme_bw() + 
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  } else {
    fiven <- stats::fivenum(df$uci, na.rm = TRUE)
    iqr <- diff(fiven[c(2, 4)])
    out <- if (!is.na(iqr)) {
      df$uci < (fiven[2L] - 1.5 * iqr) | df$uci > (fiven[4L] + 1.5 * iqr)
    }
    
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
  }

}
