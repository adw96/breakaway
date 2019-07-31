#' alpha_estimate
#' 
#' Build objects of class alpha_estimate from their components. \code{alpha_estimate()} is a constructor method
#'
#' @param estimate The estimate
#' @param error The standard error in the estimate
#' @param estimand What is the estimate trying to estimate? (richness, Shannon...)
#' @param name The name of the method
#' @param interval An interval estimate
#' @param interval_type Type of interval estimate
#' @param type TODO(Amy): Deprecate?
#' @param model  What model is fit
#' @param warnings  Any warnings?
#' @param frequentist Logical. Frequentist or Bayesian?
#' @param parametric Logical. Parametric or not?
#' @param plot A ggplot associated with the estimate
#' @param reasonable Is the estimate likely to be reasonable?
#' @param other Any other relevant objects
#' @param ... Any other objects
#'
#' @return An object of class alpha_estimate 
#'
#' @export
alpha_estimate <- function(estimate = NULL,
                           error = NULL,
                           estimand = NULL,
                           name = NULL,
                           interval = NULL,
                           interval_type = NULL,
                           type = NULL,
                           model = NULL,
                           warnings = NULL,
                           frequentist = NULL, 
                           parametric = NULL,
                           plot = NULL,
                           reasonable = NULL,
                           other = NULL,
                           ...) {
  
  # if (is.null(ci) & !is.null(error)) {
  # # TODO need f0
  #   d <- exp(1.96*sqrt(log(1 + error^2 / f0)))
  #   
  # }
  
  alpha_object <- list(estimate = estimate,
                       error = error,
                       estimand = tolower(estimand),
                       name = name,
                       interval = interval, #ifelse(is.na(estimate), c(NA, NA), interval),
                       interval_type = interval_type,
                       type = type,
                       model = model,
                       warnings = warnings,
                       frequentist = frequentist, 
                       parametric = parametric,
                       plot = plot,
                       reasonable = reasonable,
                       other = other,
                       est = estimate,
                       seest = error,
                       ci = interval,
                       ...)
  
  
  class(alpha_object) <- append("alpha_estimate", class(alpha_object))
  
  return(alpha_object)
}

#' @export
print.alpha_estimate <- function(x, ...) {
  
  if (is.null(x$estimand) | is.null(x$name) | is.null(x$estimate)) {
    cat("Attempt to print estimate with unknown estimand, name, or estimate.\n")
    cat("Stripping alpha_estimate from class and then printing gives: \n")
    class(x) <- x
    print(x)
    cat("\n")
    
  } else {
    cat(paste("Estimate of ", x$estimand,
              " from method ", x$name, ":\n", sep=""))
    cat(paste("  Estimate is ", round(x$estimate, ifelse(x$estimand == "richness", 0, 2)), "\n", sep=""))
    if (!is.null(x$error)) {
      cat(paste(" with standard error ", round(x$error, 2), "\n", sep=""))
    }
    if (!is.null(x$interval)) {
      cat(paste("  Confidence interval: (", 
                round(x$interval[1], ifelse(x$estimand == "richness", 0, 2)), ", ", 
                round(x$interval[2], ifelse(x$estimand == "richness", 0, 2)), ")\n", sep=""))
    }
    if (!is.null(x$other$cutoff)) {
      cat(paste("  Cutoff: ", x$other$cutoff)) 
    }
    cat("\n")
  }
}

#' @export
summary.alpha_estimate <- function(object, ...) {
  # output just like a list
  
  # don't plot
  y <- object
  class(y) <- setdiff(class(y), "alpha_estimate")
  y$plot <- NULL
  print(y)
}


#' @export
plot.alpha_estimate <- function(x, ...) {
  x$plot
}


