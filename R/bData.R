#' A canonical QIIME2 dataset
#'
#' Data used for internal testing
#' retrieved from https://raw.githubusercontent.com/paulinetrinh/data/master/otu_table_atacama.txt
#'
#' @format character vector
"atacama"

#' Output from test_submodel() for use in hypothesis testing vignette
#'
#' @format A list containing the bootstrap p-value, the observed F-statistic,
#' and a vector of bootstrapped F-statistics
"submodel_test"

#' DivNet model fitted to \code{soil_phylo} dataset
#'
#' DivNet model fit via call
#' #dv <- DivNet::divnet(soil_phylum, X = NULL)
#' Further details in diversity hypothesis testing vignette
#'
#' @format An object of class \code{diversityEstimates}
#' @references Whitman, T., Pepe-Ranney, C., Enders, A., Koechli, C., Campbell, A.,  Buckley, D. H., Lehmann, J. (2016). \emph{Dynamics of microbial community composition and soil organic carbon mineralization in soil following addition of pyrogenic and fresh organic matter}. The ISME journal, 10(12):2918. <doi: 10.1038/ismej.2016.68>.
"dv"

#' DivNet model fitted to \code{soil_phylo} dataset including observations on day 1
#'
#' DivNet model fit via call
#' #dv_days_1_2 <- DivNet::divnet(soil_phylum_days_1_2, X = NULL)
#' Further details in diversity hypothesis testing vignette
#'
#' @format An object of class \code{diversityEstimates}
#' @references Whitman, T., Pepe-Ranney, C., Enders, A., Koechli, C., Campbell, A.,  Buckley, D. H., Lehmann, J. (2016). \emph{Dynamics of microbial community composition and soil organic carbon mineralization in soil following addition of pyrogenic and fresh organic matter}. The ISME journal, 10(12):2918. <doi: 10.1038/ismej.2016.68>.
"dv_days_1_2"

#' DivNet model fitted to \code{soil_phylo} data from corncob
#'
#' DivNet model fit via call
#' #dv <- DivNet::divnet(soil_phylum, X = NULL)
#'
#' @format An object of class \code{diversityEstimates}
#' @references Whitman, T., Pepe-Ranney, C., Enders, A., Koechli, C., Campbell, A.,  Buckley, D. H., Lehmann, J. (2016). \emph{Dynamics of microbial community composition and soil organic carbon mineralization in soil following addition of pyrogenic and fresh organic matter}. The ISME journal, 10(12):2918. <doi: 10.1038/ismej.2016.68>.
"dv"

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
#' @source \url{https://onlinelibrary.wiley.com/doi/abs/10.1111/biom.12332}
#'
#' @references
#' Willis, A. and Bunge, J. (2015). Estimating diversity via
#' frequency ratios. \emph{Biometrics}, \bold{71}(4), 1042--1049.
#' \url{https://onlinelibrary.wiley.com/doi/abs/10.1111/biom.12332}
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
#' @source \url{https://onlinelibrary.wiley.com/doi/abs/10.1111/biom.12332}
#'
#' @references
#' Willis, A. and Bunge, J. (2015). Estimating diversity via
#' frequency ratios. \emph{Biometrics}, \bold{71}(4), 1042--1049.
#' \url{https://onlinelibrary.wiley.com/doi/abs/10.1111/biom.12332}
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
#' (Data) A toy OTU table.
#'
#' Covariate info available in `toy_metadata`. A data frame with 448 rows and 143 columns.
#' Rows give the abundance of each taxon; columns give the samples
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
