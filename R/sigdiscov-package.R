#' sigdiscov: Spatial Signature Discovery for Spatial Transcriptomics
#'
#' Unified spatial signature discovery for spatial transcriptomics data.
#' Supports Visium (10x Genomics) and single-cell resolution platforms
#' (CosMx, Xenium, MERFISH).
#'
#' @section Core Functions:
#' \itemize{
#'   \item \code{\link{compute_moran_from_lag}}: Moran's I from spatial lag
#'   \item \code{\link{compute_ind_from_lag}}: I_ND (cosine similarity)
#'   \item \code{\link{compute_delta_i}}: Signed Delta I
#'   \item \code{\link{batch_permutation_test}}: Permutation testing
#'   \item \code{\link{standardize}}, \code{\link{standardize_matrix}}: Normalization
#' }
#'
#' @section Visium Functions:
#' \itemize{
#'   \item \code{\link{load_data_visium}}: Load Space Ranger output
#'   \item \code{\link{as_data_visium}}: Convert Seurat object
#'   \item \code{\link{create_weights_visium}}: Binary weight matrix
#'   \item \code{\link{split_by_expression}}: Expression-based sender/receiver
#'   \item \code{\link{compute_signature_visium}}: Main Visium analysis
#' }
#'
#' @section Single-Cell Functions:
#' \itemize{
#'   \item \code{\link{load_data_cosmx}}: Load CosMx data
#'   \item \code{\link{load_data_anndata}}: Load H5AD files
#'   \item \code{\link{as_data_sc}}: Convert Seurat object
#'   \item \code{\link{create_weights_sc}}: Gaussian weight matrix
#'   \item \code{\link{define_sender_receiver_sc}}: Cell type + expression filter
#'   \item \code{\link{compute_signature_sc}}: Single pair analysis
#'   \item \code{\link{compute_batch_sc}}: All pairs batch analysis
#'   \item \code{\link{save_hdf5_sc}}, \code{\link{load_hdf5_sc}}: HDF5 I/O
#' }
#'
#' @section High-Level API:
#' \itemize{
#'   \item \code{\link{compute_signature}}: Auto-detect platform
#' }
#'
#' @useDynLib sigdiscov, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom Matrix rowSums colSums sparseMatrix t
#' @importFrom methods new as
#' @importFrom stats quantile sd p.adjust
#' @importFrom utils head read.csv
#' @importFrom graphics abline legend lines points
#' @keywords internal
"_PACKAGE"
