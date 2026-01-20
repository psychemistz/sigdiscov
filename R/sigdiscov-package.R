#' sigdiscov: Spatial Correlation Analysis for Spatial Transcriptomics
#'
#' Compute pairwise Moran's I statistics and signed delta I signatures
#' for spatial transcriptomics data.
#'
#' @section Main Functions:
#' \describe{
#'   \item{\code{\link{pairwise_moran}}}{Compute pairwise Moran's I between all genes}
#'   \item{\code{\link{compute_signed_delta_I}}}{Compute signed delta I signatures}
#'   \item{\code{\link{get_moran_curve}}}{Get I(r) curve for a single gene pair}
#'   \item{\code{\link{plot_moran_curve}}}{Visualize I(r) curve}
#' }
#'
#' @section Data Loading:
#' \describe{
#'   \item{\code{\link{load_visium_data}}}{Load 10x Visium data from Space Ranger output}
#'   \item{\code{\link{filter_in_tissue}}}{Filter to in-tissue spots}
#'   \item{\code{\link{vst_transform}}}{VST normalization using sctransform}
#'   \item{\code{\link{run_pipeline}}}{Complete analysis pipeline}
#' }
#'
#' @useDynLib sigdiscov, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom Matrix sparseMatrix readMM
#' @importFrom methods is as
#' @importFrom stats sd
#' @importFrom utils read.csv read.table
#' @importFrom graphics plot points lines abline legend
"_PACKAGE"
