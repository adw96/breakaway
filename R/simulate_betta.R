
#' Simulate from a fitted betta model
#'
#'
#' @param fitted_betta A fitted betta object
#' @param nsim Number of times to simulate
#' @return A list of length nsim, each element of which is a vector of
#' simulated Y-values under the fitted betta model
simulate_betta <- function(fitted_betta,
                           nsim){
  obs_vars <- fitted_betta$function.args$ses^2 + fitted_betta$ssq_u
  fitted_values <- fitted_betta$function.args$X%*%fitted_betta$table[,1,drop = F]

  colnames(fitted_values) <- NULL
  sims <- lapply(1:nsim,
  function(x)  fitted_values + rnorm(length(obs_vars),0,obs_vars))
  return(sims)
}
