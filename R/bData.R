################################################################################
#' (Data) Frequency count table of soil microbes in an apples orchard.
#'
#' @format A data frame with 88 rows and 2 variables:
#' \describe{
#'   \item{index}{an index variable}
#'   \item{frequency}{number of taxa that were observed with this frequency}
#'   ...
#' }
#' 
#' @source \url{http://onlinelibrary.wiley.com/doi/10.1111/biom.12332/abstract}
#' 
#' @references 
#' Willis, A. and Bunge, J. (2015). Estimating diversity via
#' frequency ratios. \emph{Biometrics}, \bold{71}(4), 1042--1049.
#' \url{http://onlinelibrary.wiley.com/doi/10.1111/biom.12332/abstract}
#'
#' Walsh, F. et al. (2014). (2014). Restricted streptomycin use in apple 
#' orchards did not adversely alter the soil bacteria communities. 
#' \emph{Frontiers in Microbiology} \bold{4}, 383.
#'
################################################################################
"apples"
################################################################################
################################################################################
#' (Data) Frequency count table of soil microbes in Hawaii.
#'
#' @format A data frame with 198 rows and 2 variables:
#' \describe{
#'   \item{index}{an index variable}
#'   \item{frequency}{number of taxa that were observed with this frequency}
#'   ...
#' }
#' 
#' @source \url{http://onlinelibrary.wiley.com/doi/10.1111/biom.12332/abstract}
#' 
#' @references 
#' Willis, A. and Bunge, J. (2015). Estimating diversity via
#' frequency ratios. \emph{Biometrics}, \bold{71}(4), 1042--1049.
#' \url{http://onlinelibrary.wiley.com/doi/10.1111/biom.12332/abstract}
#'
################################################################################
"hawaii"
################################################################################
################################################################################
#' (Data) Data frame of covariate information about toy_otu_table.
#'
#' @format A data frame with 143 rows and 4 variables:
#' \describe{
#'   \item{Years}{Year of sampling}
#'   \item{bloom2}{Did the sample correspond to a bloom event?}
#'   \item{Period}{What season was sampled?}
#'   \item{Site}{Where was the sample taken from?}
#'   ...
#' }
#' 
#'
################################################################################
"toy_metadata"
################################################################################
################################################################################
#' (Data) A toy OTU table. Covariate info available in `toy_metadata`.
#' 
#' @format A data frame with 448 rows and 143 columns:
#' \describe{
#'   \item{Rows:}{The abundance of each taxon...}
#'   \item{Columns:}{...in each sample}
#'   ...
#' }
#' 
################################################################################
"toy_otu_table"
################################################################################
################################################################################
#' (Data) The taxonomy of the OTUs in `toy_otu_table`.
#' 
################################################################################
"toy_taxonomy"
################################################################################