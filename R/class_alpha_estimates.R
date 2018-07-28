#' alpha_estimates
#' 
#' Build objects of class alpha_estimates from their components. \code{alpha_estimates()} is a constructor method
#'
#' @param ... Objects of class alpha_estimate, or a list of objects of class alpha_estimate
#'
#' @return An object of class alpha_estimates

#' @import tibble
#' @import ggplot2
#'
#' @export
alpha_estimates <- function(...) {
  
  alpha_ests_object <- list(...)
  
  if (alpha_ests_object %>% length == 1) {
    if (class(alpha_ests_object[[1]]) == "list") {
      alpha_ests_object <- alpha_ests_object[[1]]
    }
  } 
  
  if (!possible(alpha_ests_object)) {
    stop(paste("Attempt to make alpha_estimates object", 
               "out of a collection of objects",
               "that aren't all of class alpha_estimate"))
  } 
  
  class(alpha_ests_object) <- append("alpha_estimates", class(alpha_ests_object))
  
  return(alpha_ests_object)
}

possible <- function(candidate) {
  candidate %>% 
    lapply(class) %>%
    lapply({function(x) {"alpha_estimate" %in% x}}) %>%
    unlist %>% all
}

#' @export
print.alpha_estimates <- function(x, ...) {
  n <- length(x)
  
  cat(paste("A collection of", n, "alpha diversity estimates:\n\n"))
  
  lapply(x, print)
  return(NULL)
}

#' @export
summary.alpha_estimates <- function(object, ...) {
  # data frame of points to plot
  
  to_vector <- function(object, piece) {
    object %>% 
      lapply(function(x) {x[[piece]]}) %>%
      unlist 
  }

  get_interval <- function(x) {
    if (is.null(x[["interval"]])) {
      c(NA, NA)
    } else {
      x[["interval"]]
    }
  }
  intervals_df <- lapply(object, get_interval) %>% rbind.data.frame
  
  tibble::tibble("estimate" = to_vector(object, "estimate"),
                 "error" = to_vector(object, "error"),
                 "lower" = intervals_df[1, ] %>% c %>% unlist,
                 "upper" = intervals_df[2, ] %>% c %>% unlist)
}

#' @export
plot.alpha_estimates <- function(x, 
                                 xaxis = NULL, 
                                 symmetric = TRUE,
                                 ...) {
  
  if (is.null(xaxis)) {
    xaxis <- 1:length(x)
  }
  
  df <- summary(x) 
  
  if (symmetric) {
    # TODO (Amy) make lower bound observed richness
    df$lower_y <- pmax(0, df$estimate - 1.96*df$error)
    df$upper_y <- df$estimate + 1.96*df$error
    
  } else {
    df$lower_y <- df$lower 
    df$upper_y <- df$upper 
    
  }

  upper_ylim <- max(df$upper_y*1.2, na.rm = T)

  ggplot2::ggplot(df, aes_string(x = "xaxis")) +
    ylab("Alpha-diversity estimate") +
    geom_point(aes_string(y = "estimate")) +
    geom_linerange(aes_string(ymin = "lower_y", ymax = "upper_y")) +
    theme_bw() +
    xlim(0.5, length(x)+0.5) +
    coord_cartesian(ylim = c(0,upper_ylim))
}


