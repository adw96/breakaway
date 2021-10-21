



test_submodel <- function(fitted_betta,
                          submodel_formula){
  # use submodel_formula to derive L
  submodel_X <- stats::model.matrix(submodel_formula,
                           fitted_betta$function.args$data)

  L <- matrix(nrow = 0, ncol = ncol(fitted_betta$function.args$X))

  for(k in 1:ncol(L)){
    if( !(colnames(fitted_betta$function.args$X)[k] %in% colnames(submodel_X))
    ){
      new_row <- matrix(0, nrow = 1, ncol = ncol(L))
      new_row[1,k] <- 1
      L <- rbind(L,new_row)
    }
  }

  return(F_test(fitted_betta, L))

}
