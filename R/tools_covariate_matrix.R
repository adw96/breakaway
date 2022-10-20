#' Make design matrix
#' 
#' @param phyloseq_object A phyloseq object
#' @param variables variable names
#' 
#' @return A matrix object giving the design matrix from the desired variables.
#' 
#' @importFrom phyloseq sample_data
#' @importFrom phyloseq get_variable
#' 
#' @export
make_design_matrix <- function(phyloseq_object, variables) {
  predictors <- phyloseq_object %>% sample_data %>% get_variable(variables)
  model.matrix( ~predictors, data = predictors %>% as.data.frame)
}