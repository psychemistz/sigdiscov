#' Permutation Test for Pairwise Moran's I (Single-Cell Data)
#'
#' Computes p-values for pairwise Moran's I matrix using permutation testing.
#' Designed for single-cell spatial transcriptomics (CosMx, Xenium, MERFISH).
#'
#' For large datasets, automatically subsamples cells to make computation
#' tractable while preserving spatial structure.
#'
#' @param data Gene expression matrix (genes x cells). Will be z-normalized.
#' @param coords Cell coordinates (n x 2 matrix)
#' @param radius Radius for Gaussian weights (in coordinate units)
#' @param sigma Gaussian sigma (default: radius/3)
#' @param n_perm Number of permutations (default: 999)
#' @param max_cells Maximum cells for permutation test (default: 10000).
#'   If n > max_cells, cells are randomly subsampled.
#' @param seed Random seed for reproducibility (default: NULL)
#' @param verbose Print progress messages (default: TRUE)
#'
#' @return A list containing:
#'   \item{moran}{Observed Moran's I matrix (genes x genes)}
#'   \item{z_scores}{Z-score matrix}
#'   \item{p_values}{Two-sided p-value matrix}
#'   \item{null_sd}{Standard deviation of null distribution}
#'   \item{n_perm}{Number of permutations}
#'   \item{n_cells_used}{Number of cells used (after subsampling if applied)}
#'   \item{subsampled}{Logical, whether subsampling was applied}
#'
#' @details
#' The permutation test procedure:
#' \enumerate{
#'   \item Build sparse weight matrix W from cell coordinates
#'   \item Compute observed Moran's I: M = Z * W * Z^T / S0
#'   \item For each permutation, shuffle cell labels and recompute M
#'   \item P-value = (1 + count(|M_perm| >= |M_obs|)) / (n_perm + 1)
#' }
#'
#' For datasets larger than \code{max_cells}, random subsampling is performed
#' to maintain computational tractability. The spatial structure is preserved
#' because subsampling is uniform in space.
#'
#' @examples
#' \dontrun{
#' # Load CosMx data
#' data <- load_data_cosmx("expression.csv", "metadata.csv")
#'
#' # Run permutation test (will subsample if > 10k cells)
#' result <- moran_permutation_test_sc(
#'   data$expr,
#'   data$coords,
#'   radius = 50,
#'   n_perm = 999,
#'   seed = 42
#' )
#'
#' # Find significant pairs
#' sig_pairs <- extract_significant_pairs(result, threshold = 0.05)
#' }
#'
#' @seealso \code{\link{pairwise_moran_sc_large}}, \code{\link{extract_significant_pairs}}
#'
#' @export
moran_permutation_test_sc <- function(
    data,
    coords,
    radius,
    sigma = NULL,
    n_perm = 999,
    max_cells = 10000,
    seed = NULL,
    verbose = TRUE
) {
    # Input validation
    if (!is.matrix(data)) data <- as.matrix(data)
    if (!is.matrix(coords)) coords <- as.matrix(coords)

    n_genes <- nrow(data)
    n_cells <- ncol(data)

    stopifnot(
        "Number of cells must match" = n_cells == nrow(coords),
        "n_perm must be positive" = n_perm > 0,
        "coords must have 2 columns" = ncol(coords) == 2
    )

    if (is.null(sigma)) sigma <- radius / 3
    if (!is.null(seed)) set.seed(seed)

    # Subsample if necessary
    subsampled <- FALSE
    if (n_cells > max_cells) {
        if (verbose) {
            message(sprintf("Subsampling %d -> %d cells for permutation test",
                            n_cells, max_cells))
        }
        idx <- sample(n_cells, max_cells)
        data <- data[, idx, drop = FALSE]
        coords <- coords[idx, , drop = FALSE]
        n_cells <- max_cells
        subsampled <- TRUE
    }

    if (verbose) {
        message(sprintf("Permutation test: %d genes x %d cells, %d permutations",
                        n_genes, n_cells, n_perm))
    }

    # Create sparse weight matrix
    if (verbose) message("Creating sparse weight matrix...")
    W <- create_weights_sc(coords, coords, radius = radius, sigma = sigma)
    S0 <- sum(W)

    if (S0 < 1e-10) {
        stop("Weight matrix has no non-zero weights. Try increasing radius.")
    }

    if (verbose) {
        avg_neighbors <- length(W@x) / n_cells
        message(sprintf("W: %d non-zeros, %.1f avg neighbors/cell",
                        length(W@x), avg_neighbors))
    }

    # Z-normalize data (population SD for SpaCET compatibility)
    if (verbose) message("Z-normalizing expression data...")
    data_z <- t(apply(data, 1, function(x) {
        mu <- mean(x)
        sd_pop <- sqrt(mean((x - mu)^2))  # population SD
        if (sd_pop < 1e-10) return(rep(0, length(x)))
        (x - mu) / sd_pop
    }))

    # For sparse W, convert to dense for permutation test
    # (necessary for efficient matrix multiplication)
    if (verbose) message("Computing observed Moran's I...")

    # Use dense W for permutation test (sparse -> dense is OK for subsampled data)
    W_dense <- as.matrix(W)

    # Compute observed: M = Z * W * Z^T / S0
    observed <- tcrossprod(data_z %*% W_dense, data_z) / S0

    # Run permutation test
    if (verbose) message("Running permutation test...")

    perm_result <- run_pairwise_permutation_cpp(
        data_z, W_dense, S0, n_perm,
        abs(observed), verbose
    )

    # Compute statistics
    mean_null <- perm_result$sum / n_perm
    var_null <- perm_result$sum_sq / n_perm - mean_null^2
    sd_null <- sqrt(pmax(var_null, .Machine$double.eps))

    z_scores <- (observed - mean_null) / sd_null
    p_values <- (1 + perm_result$count_extreme) / (n_perm + 1)

    # Add gene names
    gene_names <- rownames(data)
    if (is.null(gene_names)) gene_names <- paste0("Gene", seq_len(n_genes))

    rownames(observed) <- colnames(observed) <- gene_names
    rownames(z_scores) <- colnames(z_scores) <- gene_names
    rownames(p_values) <- colnames(p_values) <- gene_names
    rownames(sd_null) <- colnames(sd_null) <- gene_names

    if (verbose) message("Done.")

    list(
        moran = observed,
        z_scores = z_scores,
        p_values = p_values,
        null_sd = sd_null,
        n_perm = n_perm,
        n_cells_used = n_cells,
        subsampled = subsampled
    )
}


#' Run Pairwise Permutation Loop (Internal, uses C++)
#'
#' @param data_z Z-normalized matrix (genes x cells)
#' @param W Dense weight matrix
#' @param S0 Weight sum
#' @param n_perm Number of permutations
#' @param abs_observed Absolute observed values (unused, kept for compatibility)
#' @param verbose Print progress
#' @return List with sum, sum_sq, count_extreme
#' @keywords internal
run_pairwise_permutation_cpp <- function(data_z, W, S0, n_perm,
                                          abs_observed, verbose) {
    # Use C++ implementation for all-pairs permutation test
    result <- allpairs_permutation_test_cpp(
        data_z, W, S0, n_perm,
        seed = 0L, verbose = verbose
    )
    list(
        sum = result$sum,
        sum_sq = result$sum_sq,
        count_extreme = result$count_extreme
    )
}


#' Permutation Test for Factor-Gene Analysis (Single-Cell)
#'
#' Tests spatial correlation between a factor gene and all other genes.
#' More efficient than pairwise test when analyzing one factor.
#'
#' @param factor_expr Factor gene expression (length n_cells)
#' @param gene_expr Gene expression matrix (genes x cells)
#' @param coords Cell coordinates (n x 2 matrix)
#' @param radius Radius for Gaussian weights
#' @param sigma Gaussian sigma (default: radius/3)
#' @param metric "moran" or "ind" (I_ND)
#' @param n_perm Number of permutations (default: 999)
#' @param seed Random seed (default: NULL)
#' @param verbose Print progress (default: TRUE)
#'
#' @return Data frame with columns: gene_idx, I_obs, p_value, z_score, p_adj
#'
#' @details
#' This function is optimized for testing one factor against many genes:
#' \enumerate{
#'   \item Computes spatial lag for each gene once
#'   \item Permutes factor expression (not gene expression)
#'   \item Uses same permutation indices across all genes
#' }
#'
#' This is more efficient than pairwise test when you only care about
#' correlations with a specific factor gene.
#'
#' @export
factor_permutation_test_sc <- function(
    factor_expr,
    gene_expr,
    coords,
    radius,
    sigma = NULL,
    metric = c("moran", "ind"),
    n_perm = 999,
    seed = NULL,
    verbose = TRUE
) {
    metric <- match.arg(metric)

    if (!is.matrix(gene_expr)) gene_expr <- as.matrix(gene_expr)
    if (!is.matrix(coords)) coords <- as.matrix(coords)

    n_genes <- nrow(gene_expr)
    n_cells <- length(factor_expr)

    stopifnot(
        "Cells must match" = n_cells == ncol(gene_expr),
        "Cells must match coords" = n_cells == nrow(coords)
    )

    if (is.null(sigma)) sigma <- radius / 3
    if (!is.null(seed)) set.seed(seed)

    if (verbose) {
        message(sprintf("Factor permutation test: %d genes x %d cells",
                        n_genes, n_cells))
    }

    # Create sparse weight matrix
    if (verbose) message("Creating weight matrix...")
    W <- create_weights_sc(coords, coords, radius = radius, sigma = sigma)

    # Z-normalize
    z_factor <- (factor_expr - mean(factor_expr)) /
                sqrt(mean((factor_expr - mean(factor_expr))^2))

    gene_expr_z <- t(apply(gene_expr, 1, function(x) {
        mu <- mean(x)
        sd_pop <- sqrt(mean((x - mu)^2))
        if (sd_pop < 1e-10) return(rep(0, length(x)))
        (x - mu) / sd_pop
    }))

    # Run batch permutation test
    if (verbose) message("Running permutation test...")
    result <- batch_permutation_test(
        z_f = z_factor,
        Z_g = t(gene_expr_z),  # n_cells x n_genes
        W = W,
        metric = metric,
        n_perm = n_perm,
        weight_sum = sum(W)
    )

    # Add adjusted p-values
    result$p_adj <- p.adjust(result$p_value, method = "BH")

    # Add gene names
    gene_names <- rownames(gene_expr)
    if (!is.null(gene_names)) {
        result$gene <- gene_names[result$gene_idx]
    }

    if (verbose) {
        n_sig <- sum(result$p_adj < 0.05, na.rm = TRUE)
        message(sprintf("Significant genes (FDR < 0.05): %d", n_sig))
    }

    result
}


#' Streaming Permutation Test for Pairwise Moran's I (Memory Efficient)
#'
#' Performs permutation test without storing dense weight matrix in memory.
#' Uses KD-tree neighbor lists to compute spatial lag on-the-fly for each
#' permutation. Suitable for large datasets (>10k cells) where the dense
#' approach would run out of memory.
#'
#' @param data Gene expression matrix (genes x cells). Will be z-normalized.
#' @param coords Cell coordinates (n x 2 matrix)
#' @param radius Radius for Gaussian weights (in coordinate units)
#' @param sigma Gaussian sigma (default: radius/3)
#' @param n_perm Number of permutations (default: 999)
#' @param seed Random seed for reproducibility (default: NULL)
#' @param verbose Print progress messages (default: TRUE)
#'
#' @return A list containing:
#'   \item{moran}{Observed Moran's I matrix (genes x genes)}
#'   \item{z_scores}{Z-score matrix}
#'   \item{p_values}{Two-sided p-value matrix}
#'   \item{null_sd}{Standard deviation of null distribution}
#'   \item{n_perm}{Number of permutations}
#'   \item{n_cells}{Number of cells}
#'   \item{weight_sum}{Sum of spatial weights (S0)}
#'   \item{n_edges}{Number of neighbor pairs}
#'
#' @details
#' Memory comparison for 98k cells:
#' \itemize{
#'   \item Dense approach: ~77 GB (stores full weight matrix)
#'   \item Streaming approach: ~1.7 GB (stores only neighbor lists)
#' }
#'
#' The streaming approach stores sparse neighbor lists and recomputes
#' spatial lag for each permutation. This trades computation time for
#' memory efficiency.
#'
#' Computational cost per permutation:
#' \itemize{
#'   \item O(n_cells × avg_neighbors × n_genes) for spatial lag
#'   \item O(n_genes² × n_cells) for Moran's I matrix
#' }
#'
#' @examples
#' \dontrun{
#' # Load CosMx data (98k cells)
#' data <- load_data_cosmx("expression.csv", "metadata.csv")
#'
#' # Run streaming permutation test (no subsampling needed)
#' result <- moran_permutation_test_sc_streaming(
#'   data$expr,
#'   data$coords,
#'   radius = 0.1,  # 100 um
#'   n_perm = 999,
#'   seed = 42
#' )
#'
#' # Find significant pairs
#' sig_pairs <- extract_significant_pairs(result, threshold = 0.05)
#' }
#'
#' @seealso \code{\link{moran_permutation_test_sc}} for dense approach with subsampling
#'
#' @export
moran_permutation_test_sc_streaming <- function(
    data,
    coords,
    radius,
    sigma = NULL,
    n_perm = 999,
    seed = NULL,
    verbose = TRUE
) {
    # Input validation
    if (!is.matrix(data)) data <- as.matrix(data)
    if (!is.matrix(coords)) coords <- as.matrix(coords)

    n_genes <- nrow(data)
    n_cells <- ncol(data)

    stopifnot(
        "Number of cells must match" = n_cells == nrow(coords),
        "n_perm must be positive" = n_perm > 0,
        "coords must have 2 columns" = ncol(coords) == 2
    )

    if (is.null(sigma)) sigma <- radius / 3
    seed_val <- if (is.null(seed)) 0L else as.integer(seed)

    # Call C++ streaming permutation test
    result <- allpairs_permutation_streaming_cpp(
        data, coords, radius, sigma,
        as.integer(n_perm), seed_val, verbose
    )

    # Compute statistics from accumulated values
    observed <- result$observed
    mean_null <- result$sum / n_perm
    var_null <- result$sum_sq / n_perm - mean_null^2
    sd_null <- sqrt(pmax(var_null, .Machine$double.eps))

    z_scores <- (observed - mean_null) / sd_null
    p_values <- (1 + result$count_extreme) / (n_perm + 1)

    # Add gene names
    gene_names <- rownames(data)
    if (is.null(gene_names)) gene_names <- paste0("Gene", seq_len(n_genes))

    rownames(observed) <- colnames(observed) <- gene_names
    rownames(z_scores) <- colnames(z_scores) <- gene_names
    rownames(p_values) <- colnames(p_values) <- gene_names
    rownames(sd_null) <- colnames(sd_null) <- gene_names

    list(
        moran = observed,
        z_scores = z_scores,
        p_values = p_values,
        null_sd = sd_null,
        n_perm = n_perm,
        n_cells = n_cells,
        weight_sum = result$weight_sum,
        n_edges = result$n_edges
    )
}
