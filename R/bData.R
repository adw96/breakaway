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
#' 
#' @references 
#' Willis, A. and Bunge, J. (2015). Estimating diversity via
#' frequency ratios. \emph{Biometrics}, \bold{71}(4), 1042--1049.
#' \doi{doi:10.1111/biom.12332}
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
#' @references 
#' Willis, A. and Bunge, J. (2015). Estimating diversity via
#' frequency ratios. \emph{Biometrics}, \bold{71}(4), 1042--1049.
#' \doi{doi:10.1111/biom.12332}
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
################################################################################
#' (Data) Data frame of covariate information about pasolli_et_al.
#'
#' @format A data frame with 4930 rows and 9 variables:
#' \describe{
#'   \item{SGB.ID}{sample ID}
#'   \item{AverageAbundance}{average abundance in the sample}
#'   \item{std}{standard deviation}
#'   \item{q1}{first quartile}
#'   \item{MedianAbundance}{median abundance}
#'   \item{q3}{third quartile}
#'   \item{min}{minimum}
#'   \item{max}{maximum}
#'   \item{#.Samples}{number of samples}
#'   ...
#' }
#' 
#'
################################################################################
"pasolli_et_al"
################################################################################
################################################################################
#' (Data) Data frame of soil data from whitman_et_al.
#'
#' A phyloseq object with an OTU table and sample data from a soil microbiome study.
#'
#' @format A phyloseq-class experiment-level object with an OTU table and sample data.
#' \describe{
#' \item{otu_table}{OTU table with 7,770 taxa and 119 samples}
#' \item{tax_table}{taxonomy table}
#' \item{sam_data}{sample data with the following covariates:
#' \itemize{
#' \item \code{Plants}, values \code{0} and \code{1}. Index for different plants
#' \item \code{Day}, values \code{0} (initial sampling point), \code{1} (12 days after treatment additions), and \code{2} (82 days after treatment additions). Index for different days of measurement
#' \item \code{Amdmt}, values \code{0} (no additions), \code{1} (biochar additions), and \code{2} (fresh biomass additions). Index for different soil additives.
#' \item \code{DayAmdmt}, values \code{00}, \code{01}, \code{02}, \code{10}, \code{11}, \code{12}, \code{20}, \code{21}, and \code{22}. A single index for the combination of \code{Day} and \code{Amdmt} with \code{Day} as the first digit and \code{Amdmt} as the second digit.
#' \item \code{ID}, values \code{A}, \code{B}, \code{C}, \code{D},  and \code{F}. Index for different soil plots.
#' }}
#' ...
#' }
#' @references Whitman, T., Pepe-Ranney, C., Enders, A., Koechli, C., Campbell, A.,  Buckley, D. H., Lehmann, J. (2016). \emph{Dynamics of microbial community composi-tion and soil organic carbon mineralization in soil following addition of pyrogenic andfresh organic matter}. The ISME journal, 10(12):2918. <doi: 10.1038/ismej.2016.68>.
################################################################################
"soil_phylo"
################################################################################