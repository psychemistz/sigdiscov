#' Permutation Test for Spatial Correlation (Single Gene)
#'
#' Performs a permutation test to assess the significance of spatial
#' correlation between a factor and gene expression.
#'
#' @param z_f Numeric vector. Standardized factor expression (length n).
#' @param lag_g Numeric vector. Pre-computed spatial lag (W * z_g).
#' @param metric Character. Either "moran" or "ind".
#' @param n_perm Integer. Number of permutations (default: 999).
#' @param seed Integer. Random seed for reproducibility (optional).
#'
#' @return A list with components:
#' \describe{
#'   \item{I_obs}{Observed statistic value}
#'   \item{p_value}{Two-sided p-value (bounded between 0 and 1)}
#'   \item{z_score}{Z-score relative to null distribution}
#'   \item{null_mean}{Mean of null distribution}
#'   \item{null_sd}{Standard deviation of null distribution}
#'   \item{n_perm}{Number of permutations performed}
#' }
#'
#' @details
#' The null hypothesis is that there is no spatial association between
#' the factor and gene expression. Under the null, factor expression
#' is randomly distributed in space.
#'
#' The permutation strategy:
#' \enumerate{
#'   \item Compute the observed statistic
#'   \item Keep spatial lag fixed (expensive computation done once)
#'   \item Permute factor expression B times
#'   \item p-value = (1 + count(|I_perm| >= |I_obs|)) / (B + 1)
#' }
#'
#' The continuity correction ensures p-values are never exactly 0.
#'
#' @examples
#' set.seed(42)
#' n <- 100
#' z_f <- rnorm(n)
#' lag_g <- rnorm(n)
#'
#' # Test with I_ND metric
#' result <- permutation_test_core(z_f, lag_g, "ind", n_perm = 999)
#' result$p_value
#'
#' # P-value is always valid
#' stopifnot(result$p_value >= 0 && result$p_value <= 1)
#'
#' @seealso \code{\link{batch_permutation_test}} for testing multiple genes
#'
#' @export
permutation_test_core <- function(z_f, lag_g, metric = c("moran", "ind"),
                                   n_perm = 999L, seed = NULL) {
    metric <- match.arg(metric)

    if (!is.null(seed)) {
        set.seed(seed)
    }

    permutation_test_cpp(z_f, lag_g, metric, as.integer(n_perm))
}

#' Batch Permutation Test for Multiple Genes
#'
#' Efficiently performs permutation tests for one factor against all genes.
#' Pre-generates permutation indices and reuses them across genes for
#' consistency and efficiency.
#'
#' @param z_f Numeric vector. Standardized factor expression (length n_obs).
#' @param Z_g Numeric matrix. Gene expression matrix (n_obs x n_genes),
#'   genes in columns.
#' @param W Sparse matrix. Weight matrix for computing spatial lags.
#' @param metric Character. Either "moran" or "ind".
#' @param n_perm Integer. Number of permutations (default: 999).
#' @param seed Integer. Random seed for reproducibility (optional).
#'
#' @return A data.frame with columns:
#' \describe{
#'   \item{gene_idx}{Gene index (1-based)}
#'   \item{I_obs}{Observed statistic}
#'   \item{p_value}{Two-sided p-value (bounded between 0 and 1)}
#'   \item{z_score}{Z-score relative to null}
#' }
#'
#' @details
#' This function is optimized for genome-wide analysis:
#' \enumerate{
#'   \item Pre-generates all permutation indices once
#'   \item Reuses permutations across all genes (ensures consistency)
#'   \item Computes spatial lag for each gene on-the-fly
#'   \item Checks for user interrupts every 100 genes
#' }
#'
#' For a typical analysis with 20,000 genes and 999 permutations,
#' this processes approximately 20 million permutation statistics.
#'
#' @examples
#' \dontrun{
#' library(Matrix)
#' set.seed(42)
#'
#' # Simulate data
#' n <- 100
#' n_genes <- 50
#' z_f <- rnorm(n)
#' Z_g <- matrix(rnorm(n * n_genes), n, n_genes)
#' W <- rsparsematrix(n, n, density = 0.1)
#' W <- W / rowSums(abs(W))  # row normalize
#'
#' # Run batch test
#' result <- batch_permutation_test(z_f, Z_g, W, "ind", n_perm = 999)
#'
#' # All p-values are valid
#' stopifnot(all(result$p_value >= 0 & result$p_value <= 1, na.rm = TRUE))
#'
#' # Adjust for multiple testing
#' result$p_adj <- adjust_pvalues(result$p_value, "BH")
#' }
#'
#' @seealso \code{\link{permutation_test_core}} for single gene test,
#'   \code{\link{adjust_pvalues}} for multiple testing correction
#'
#' @export
batch_permutation_test <- function(z_f, Z_g, W, metric = c("moran", "ind"),
                                    n_perm = 999L, seed = NULL) {
    metric <- match.arg(metric)

    if (!is.null(seed)) {
        set.seed(seed)
    }

    # Ensure Z_g is a matrix with genes in columns
    if (!is.matrix(Z_g)) {
        Z_g <- as.matrix(Z_g)
    }

    # Call C++ implementation
    batch_permutation_test_cpp(z_f, Z_g, W, metric, as.integer(n_perm))
}
