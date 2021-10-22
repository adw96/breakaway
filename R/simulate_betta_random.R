
#' Simulate from a fitted betta_random model
#'
#'
#' @param fitted_betta A fitted betta_random object
#' @param nsim Number of times to simulate
#' @return A list of length nsim, each element of which is a vector of
#' simulated Y-values under the fitted betta model
simulate_betta_random <- function(fitted_betta,
                           nsim){
  obs_vars <- unlist(fitted_betta$function.args$ses)^2 + fitted_betta$ssq_u
  ssq_group <- fitted_betta$ssq_group
  groups <- fitted_betta$function.args$groups
  groups <- unlist(groups)
  unique_groups <- names(ssq_group)
  ngroups <- length(unique_groups)
  prematrix <- lapply(unique_groups,function(x) as.numeric(groups == x))
  group_matrix <- do.call(rbind,prematrix)
  ssq_group <- matrix(ssq_group,nrow = 1)

  fitted_values <- fitted_betta$function.args$X%*%fitted_betta$table[,1,drop = F]
  colnames(fitted_values) <- NULL
  sims <- lapply(1:nsim,
                 function(x)  fitted_values + rnorm(length(obs_vars),0,obs_vars) +
                   as.numeric(ssq_group %*% group_matrix)
                   )
  return(sims)
}
