#' Plot function for alpha_estimates class
#' 
#' @param x Object of class \code{alpha_estimates}.
#' @param data (Optional). Default \code{NULL}. Required object of class \code{phyloseq} if including a \code{sample_data} variable for \code{color} or \code{shape}.
#' @param measure (Optional). If there are multiple richness measures included in \code{x}, this can be set to the the desired measure to be plotted. Defaults to the measure of the first alpha diversity estimate.
#' @param color (Optional). Default \code{NULL}. The sample variable to map to different colors. Can be a single character string of the variable name in \code{sample_data} or a custom supplied vector with length equal to the number of samples.
#' @param shape (Optional). Default \code{NULL}. The sample variable to map to different shapes. Can be a single character string of the variable name in \code{sample_data} or a custom supplied vector with length equal to the number of samples.
#' @param title (Optional). Default NULL. Character string. The main title for the graphic.
#' @param trim_plot (Optional). Default \code{FALSE}. Boolean indicator for whether you want the plot to include the full confidence intervals. 
#' @param ... See details
#' 
#' @details ... does not currently have any implemented options. Optional arguments currently include "trim_plot", a  Optional
#' 
#' @import magrittr
#' 
#' @examples
#' \dontrun{
#' library(phyloseq)
#' data(GlobalPatterns)
#' alphas <- breakaway(GlobalPatterns)
#' plot(alphas)
#' }
#' @export
plot.alpha_estimates <- function(x, data = NULL, measure = NULL, color = NULL, shape = NULL, title = NULL, trim_plot = FALSE, ...) {
  
  if (is.null(measure)) {
    all_measures <- x %>% lapply(function(x) x$name) %>% unlist %>% unique
    measure <- all_measures[1]
  }
  
  df <- summary(x) %>% 
    dplyr::mutate("Sample" = data %>% sample_names)
  
  if (all(is.na(df$estimate))) {
    stop("There are no estimates in this alpha_estimates object!")
  }
  
  if (!is.null(color)) {
    if (color %in% phyloseq::sample_variables(data)) {
      df$color <- phyloseq::get_variable(data, color)
    } else if (length(color) == nrow(df)) {
      df$color <- color
    } else {
      stop("color must either match a variable or be a custom vector of correct length!")
    }
  } # End if (!is.null(color))
  
  if (!is.null(shape)) {
    if (shape %in% phyloseq::sample_variables(data)) {
      df$shape <- phyloseq::get_variable(data, shape)
    } else if (length(shape) == nrow(df)) {
      df$shape <- shape
    } else {
      stop("shape must either match a variable or be a custom vector of correct length!")
    }
  } # End if (!is.null(color))
  
  yname1 <- measure
  yname2 <- x[[1]]$estimand
  
  my_gg <- ggplot2::ggplot(df, ggplot2::aes(x = names, xend = names, 
                                            colour = color, shape = shape)) +
    ggplot2::geom_point(ggplot2::aes(x = names, y = estimate)) + 
    ggplot2::geom_segment(ggplot2::aes(y = lower, yend = upper)) +
    ggplot2::ylab(paste(yname1, "estimate of", yname2)) +
    ggplot2::xlab("") +
    ggplot2::labs(title = title) +
    ggplot2::theme_bw() + 
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  
  if (!trim_plot) {
    fiven <- stats::fivenum(df$uci, na.rm = TRUE)
    iqr <- diff(fiven[c(2, 4)])
    if (!is.na(iqr)) {
      out <- df$uci < (fiven[2L] - 1.5 * iqr) | df$uci > (fiven[4L] + 1.5 * iqr)
      ylower <- min(0, 0.95*min(df$uci[!out]))
      yupper <- 1.05*max(df$uci[!out])
      
      my_gg <- my_gg +
        ggplot2::coord_cartesian(ylim = c(ylower,yupper)) 
    } 
  } 
  
  my_gg
}
