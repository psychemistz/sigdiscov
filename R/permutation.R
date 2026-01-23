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
#' @param weight_sum Numeric. Sum of weights for normalization. If NULL or <= 0,
#'   uses n_obs for normalization (default: NULL).
#' @param seed Integer. Random seed for reproducibility (optional).
#'
#' @return A data.frame with columns:
#' \describe{
#'   \item{gene_idx}{Gene index (1-based)}
#'   \item{I_obs}{Observed statistic (Moran's I or I_ND)}
#'   \item{p_value}{Two-sided p-value (bounded between 0 and 1)}
#'   \item{z_score}{Z-score relative to null}
#'   \item{null_sd}{Standard deviation of the null distribution}
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
                                    n_perm = 999L, weight_sum = NULL,
                                    seed = NULL) {
    metric <- match.arg(metric)

    if (!is.null(seed)) {
        set.seed(seed)
    }

    # Ensure Z_g is a matrix with genes in columns
    if (!is.matrix(Z_g)) {
        Z_g <- as.matrix(Z_g)
    }

    # Use weight_sum for normalization if provided
    ws <- if (is.null(weight_sum) || weight_sum <= 0) -1.0 else weight_sum

    # Call C++ implementation
    batch_permutation_test_cpp(z_f, Z_g, W, metric, as.integer(n_perm), ws)
}


# =============================================================================
# ALL-PAIRS PERMUTATION TESTING (Pairwise Moran's I)
# =============================================================================

#' Permutation Test for Pairwise Moran's I Statistical Significance
#'
#' Computes p-values and z-scores for the full pairwise Moran's I matrix using
#' permutation testing. The null hypothesis is that there is no spatial
#' autocorrelation (values are randomly distributed across spots).
#'
#' @param data Gene expression matrix (genes x spots). Will be z-normalized
#'   internally if not already normalized.
#' @param spot_coords Spot coordinates matrix (n x 2) with columns for row/col
#'   or x/y coordinates.
#' @param n_perm Number of permutations. Default: 1000. Higher values give
#'   more precise p-values but take longer.
#' @param seed Random seed for reproducibility. Default: NULL (random).
#' @param max_radius Maximum grid radius for neighbor search. Default: 5.
#' @param platform Platform type: "visium" or "old". Default: "visium".
#' @param same_spot Include same-spot weights (diagonal). Default: FALSE.
#' @param return_null Return null distribution statistics. Default: FALSE.
#' @param verbose Print progress messages. Default: TRUE.
#'
#' @return A list containing:
#'   \item{moran}{Observed Moran's I matrix (genes x genes)}
#'   \item{z_scores}{Z-score matrix from permutation null}
#'   \item{p_values}{Two-sided p-value matrix}
#'   \item{null_sd}{Standard deviation of null distribution}
#'   \item{n_perm}{Number of permutations performed}
#'   \item{null_mean}{Mean of null distribution (if return_null = TRUE)}
#'
#' @details
#' The permutation test works as follows:
#' \enumerate{
#'   \item Compute observed Moran's I for all gene pairs: M = Z W Z^T / S0
#'   \item For each permutation, shuffle spot labels and recompute Moran's I
#'   \item Calculate z-score: (observed - mean(permuted)) / sd(permuted)
#'   \item Calculate p-value: proportion of |permuted| >= |observed|
#' }
#'
#' The p-value formula includes a continuity correction:
#' p = (1 + count_extreme) / (n_perm + 1)
#'
#' @examples
#' \dontrun{
#' # Load Visium data
#' visium <- load_visium_data("spaceranger_output")
#' data_vst <- vst_transform(visium$counts)
#' spot_coords <- get_spot_coords(visium)
#'
#' # Run permutation test
#' result <- moran_permutation_test(
#'   data_vst,
#'   spot_coords,
#'   n_perm = 1000,
#'   seed = 42
#' )
#'
#' # Find significant pairs (FDR < 0.05)
#' p_adj <- adjust_moran_pvalues(result$p_values, method = "BH")
#' sig_pairs <- which(p_adj < 0.05 & lower.tri(p_adj), arr.ind = TRUE)
#' }
#'
#' @seealso \code{\link{batch_permutation_test}}, \code{\link{adjust_moran_pvalues}},
#'   \code{\link{extract_significant_pairs}}
#'
#' @export
moran_permutation_test <- function(
    data,
    spot_coords,
    n_perm = 1000,
    seed = NULL,
    max_radius = 5,
    platform = "visium",
    same_spot = FALSE,
    return_null = FALSE,
    verbose = TRUE
) {
    # Input validation
    if (!is.matrix(data)) {
        data <- as.matrix(data)
    }
    if (!is.matrix(spot_coords)) {
        spot_coords <- as.matrix(spot_coords)
    }

    stopifnot(
        "Number of spots must match" = ncol(data) == nrow(spot_coords),
        "n_perm must be at least 100" = n_perm >= 100,
        "spot_coords must have 2 columns" = ncol(spot_coords) == 2
    )

    # Set seed if provided
    if (!is.null(seed)) {
        set.seed(seed)
    }

    n_genes <- nrow(data)
    n_spots <- ncol(data)

    if (verbose) {
        message(sprintf(
            "Permutation test: %d genes, %d spots, %d permutations",
            n_genes, n_spots, n_perm
        ))
    }

    # Create weight matrix using existing package function
    if (verbose) message("Creating weight matrix...")
    W_result <- create_weight_matrix(spot_coords, max_radius, platform, same_spot)
    W <- W_result$W
    S0 <- W_result$weight_sum

    if (S0 == 0) {
        stop("Weight matrix has no non-zero weights. Try increasing max_radius.")
    }

    # Z-normalize data
    if (verbose) message("Z-normalizing expression data...")
    data_z <- t(apply(data, 1, function(x) {
        sd_x <- sd(x)
        if (sd_x < .Machine$double.eps) {
            return(rep(0, length(x)))
        }
        (x - mean(x)) / sd_x
    }))

    # Compute observed Moran's I: M = Z * W * Z^T / S0
    if (verbose) message("Computing observed Moran's I...")
    observed <- tcrossprod(data_z %*% W, data_z) / S0

    # Run permutation test
    if (verbose) message("Running permutations...")
    perm_result <- run_permutation_loop_internal(
        data_z, W, S0, n_perm,
        abs(observed), verbose
    )

    # Compute null distribution statistics
    mean_null <- perm_result$sum / n_perm
    var_null <- perm_result$sum_sq / n_perm - mean_null^2
    sd_null <- sqrt(pmax(var_null, .Machine$double.eps))

    # Compute z-scores
    z_scores <- (observed - mean_null) / sd_null

    # Compute p-values (two-sided with continuity correction)
    p_values <- (1 + perm_result$count_extreme) / (n_perm + 1)

    # Add row/column names
    gene_names <- rownames(data)
    if (is.null(gene_names)) {
        gene_names <- paste0("Gene", seq_len(n_genes))
    }

    rownames(observed) <- colnames(observed) <- gene_names
    rownames(z_scores) <- colnames(z_scores) <- gene_names
    rownames(p_values) <- colnames(p_values) <- gene_names
    rownames(sd_null) <- colnames(sd_null) <- gene_names

    # Build result
    result <- list(
        moran = observed,
        z_scores = z_scores,
        p_values = p_values,
        null_sd = sd_null,
        n_perm = n_perm
    )

    if (return_null) {
        rownames(mean_null) <- colnames(mean_null) <- gene_names
        result$null_mean <- mean_null
    }

    result
}


#' Run Permutation Loop (Internal)
#'
#' Internal function to run the permutation loop using pure R implementation.
#' Can be replaced with C++ implementation for better performance.
#'
#' @param data_z Z-normalized expression matrix (genes x spots)
#' @param W Weight matrix (spots x spots)
#' @param S0 Sum of weights
#' @param n_perm Number of permutations
#' @param abs_observed Absolute value of observed Moran's I matrix
#' @param verbose Print progress messages
#'
#' @return List with sum, sum_sq, count_extreme matrices
#'
#' @keywords internal
run_permutation_loop_internal <- function(data_z, W, S0, n_perm,
                                           abs_observed, verbose) {
    n_genes <- nrow(data_z)
    n_spots <- ncol(data_z)

    # Initialize accumulators
    perm_sum <- matrix(0, n_genes, n_genes)
    perm_sum_sq <- matrix(0, n_genes, n_genes)
    count_extreme <- matrix(0, n_genes, n_genes)

    # Progress reporting interval
    report_interval <- max(1, n_perm %/% 10)

    for (p in seq_len(n_perm)) {
        if (verbose && (p %% report_interval == 0)) {
            message(sprintf("  Permutation %d / %d", p, n_perm))
        }

        # Permute columns (spots)
        perm_idx <- sample.int(n_spots)
        data_perm <- data_z[, perm_idx, drop = FALSE]

        # Compute permuted Moran's I
        moran_perm <- tcrossprod(data_perm %*% W, data_perm) / S0

        # Update accumulators
        perm_sum <- perm_sum + moran_perm
        perm_sum_sq <- perm_sum_sq + moran_perm * moran_perm
        count_extreme <- count_extreme + (abs(moran_perm) >= abs_observed)
    }

    list(
        sum = perm_sum,
        sum_sq = perm_sum_sq,
        count_extreme = count_extreme
    )
}


#' Multiple Testing Correction for Moran's I P-values
#'
#' Applies multiple testing correction to p-values from permutation test.
#' Only adjusts the unique pairs (lower triangle) to avoid over-correction.
#'
#' @param p_values P-value matrix from \code{moran_permutation_test}
#' @param method Correction method passed to \code{p.adjust}:
#'   "BH" (Benjamini-Hochberg, default), "bonferroni", "holm", "BY", etc.
#' @param include_diagonal Include diagonal (univariate p-values) in
#'   correction. Default: FALSE.
#'
#' @return Matrix of adjusted p-values with same dimensions as input.
#'
#' @examples
#' \dontrun{
#' result <- moran_permutation_test(data, coords, n_perm = 1000)
#' p_adj <- adjust_moran_pvalues(result$p_values, method = "BH")
#'
#' # Find significant pairs at FDR < 0.05
#' sig <- which(p_adj < 0.05 & lower.tri(p_adj), arr.ind = TRUE)
#' }
#'
#' @export
adjust_moran_pvalues <- function(p_values, method = "BH",
                                  include_diagonal = FALSE) {
    stopifnot(is.matrix(p_values))
    stopifnot(nrow(p_values) == ncol(p_values))

    n <- nrow(p_values)

    # Extract lower triangle (unique pairs)
    if (include_diagonal) {
        mask <- lower.tri(p_values, diag = TRUE)
    } else {
        mask <- lower.tri(p_values, diag = FALSE)
    }

    p_vec <- p_values[mask]

    # Apply correction
    p_adj_vec <- p.adjust(p_vec, method = method)

    # Reconstruct symmetric matrix
    result <- matrix(NA, n, n)
    result[mask] <- p_adj_vec
    result <- t(result)
    result[mask] <- p_adj_vec

    # Copy names
    rownames(result) <- rownames(p_values)
    colnames(result) <- colnames(p_values)

    result
}


#' Extract Significant Gene Pairs from Permutation Test
#'
#' Convenience function to extract significant gene pairs from permutation
#' test results as a data frame.
#'
#' @param perm_result Result from \code{moran_permutation_test}
#' @param threshold P-value or FDR threshold. Default: 0.05.
#' @param method Multiple testing correction method. Default: "BH".
#' @param min_moran Minimum absolute Moran's I to include. Default: 0.
#'
#' @return Data frame with columns: gene1, gene2, moran_i, z_score, p_value,
#'   p_adjusted.
#'
#' @examples
#' \dontrun{
#' result <- moran_permutation_test(data, coords, n_perm = 1000, seed = 42)
#' sig_pairs <- extract_significant_pairs(result, threshold = 0.05)
#' head(sig_pairs)
#' }
#'
#' @export
extract_significant_pairs <- function(
    perm_result,
    threshold = 0.05,
    method = "BH",
    min_moran = 0
) {
    # Adjust p-values
    p_adj <- adjust_moran_pvalues(perm_result$p_values, method = method)

    # Get lower triangle indices
    mask <- lower.tri(perm_result$moran, diag = FALSE)

    # Find significant pairs
    sig_mask <- mask & (p_adj <= threshold) & (abs(perm_result$moran) >= min_moran)
    sig_idx <- which(sig_mask, arr.ind = TRUE)

    if (nrow(sig_idx) == 0) {
        message("No significant pairs found at threshold ", threshold)
        return(data.frame(
            gene1 = character(0),
            gene2 = character(0),
            moran_i = numeric(0),
            z_score = numeric(0),
            p_value = numeric(0),
            p_adjusted = numeric(0)
        ))
    }

    gene_names <- rownames(perm_result$moran)

    data.frame(
        gene1 = gene_names[sig_idx[, 1]],
        gene2 = gene_names[sig_idx[, 2]],
        moran_i = perm_result$moran[sig_idx],
        z_score = perm_result$z_scores[sig_idx],
        p_value = perm_result$p_values[sig_idx],
        p_adjusted = p_adj[sig_idx],
        stringsAsFactors = FALSE
    )
}
