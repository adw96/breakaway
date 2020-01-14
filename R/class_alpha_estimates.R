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
  
  # Unfortunately I don't know how to 
  for (i in 1:n) {
    cat(paste("$", names(x)[i]), "\n", sep = "")
    print(x[[i]])
    cat(paste("\n"))
  }
  # lapply(x, print)
  return(NULL)
}

#' @export
summary.alpha_estimates <- function(object, ...) {
  
  dots <- list(...)
  
  to_vector <- function(object, piece) {
    unlisted_object <- object %>% 
      lapply(function(x) {x[[piece]]}) %>%
      unlist 
    
    if (unlisted_object %>% is.null %>% all) {
      rep(NA, object %>% length)
    } else {
      unlisted_object
    }
  }
  
  get_interval <- function(x) {
    if (is.null(x[["interval"]])) {
      c(NA, NA)
    } else {
      x[["interval"]]
    }
  }
  intervals_df <- lapply(object, get_interval) %>% rbind.data.frame
  
  tb <- tibble::tibble("estimate" = to_vector(object, "estimate"),
                       "error" = to_vector(object, "error"),
                       "lower" = intervals_df[1, ] %>% c %>% unlist,
                       "upper" = intervals_df[2, ] %>% c %>% unlist)
  
  if (!is.null(names(object))) {
    tb %<>%
      tibble::add_column("sample_names" = names(object))
  }
  
  tb %<>%
    tibble::add_column("name" = lapply(object, function(x) x$name) %>% unlist)
  
  tb %<>%
    tibble::add_column("model" = lapply(object, function(x) x$model) %>% unlist)
  
  
  tb 
  
}
