################################################################################
#' Build alpha-class objects from their components.
#'
#' \code{alpha()} is a constructor method
#'
#' @usage alpha(...)
#'
#' @param ... One or more components
#'
#' @return An object of class alpha. 
#'
#' @seealso \code{\link{breakaway}}
alpha <- function(estimate = NULL,
                  error = NULL,
                  interval = NULL,
                  type = NULL,
                  model = NULL,
                  warnings = NULL
                  frequentist = NULL, # frequentist or Bayesian
                  parametric = NULL,
                  other = NULL
                  ){
  
  alpha_object <- list()
  alpha_object@estimate <- estimate
  alpha_object@error <- error
  
  
  class(alpha_object) <- append("alpha", class(alpha_object))
  
  return(alpha_object)
}

print.alpha <- function(my_alpha) {
  alpha.df <- data.frame("estimate" = my_alpha$estimate, "error" = my_alpha$estimate)
  print(alpha.df)
}


