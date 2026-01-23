#' Normalize Expression Data
#'
#' Normalize raw count data using various methods. LogNorm is recommended for
#' large single-cell datasets (CosMx, Xenium, MERFISH) where VST is too slow.
#'
#' @param counts A matrix or sparse matrix of raw counts (genes x cells/spots).
#' @param method Character. Normalization method: "lognorm" (default), "logcpm",
#'   or "vst".
#' @param scale_factor Numeric. Scale factor for lognorm/logcpm (default: 10000
#'   for lognorm, 1e6 for logcpm).
#' @param verbose Logical. Print progress messages (default: TRUE).
#'
#' @return Normalized expression matrix (same dimensions as input).
#'
#' @details
#' Three normalization methods are supported:
#'
#' **lognorm (recommended for large datasets)**:
#' \deqn{y_{ij} = \log(1 + \frac{x_{ij}}{\sum_j x_{ij}} \times 10000)}
#' Fast and provides ~80% of VST performance for spatial correlation analysis.
#' Default in Seurat for single-cell data.
#'
#' **logcpm**:
#' \deqn{y_{ij} = \log(1 + \frac{x_{ij}}{\sum_j x_{ij}} \times 10^6)}
#' Similar to lognorm but with scale factor of 1 million (CPM).
#'
#' **vst**:
#' Variance stabilizing transformation using sctransform. Most accurate but
#' computationally expensive for large datasets. Requires sctransform package.
#'
#' Performance comparison (CytoSig correlation):
#' \itemize{
#'   \item VST: 0.051 (baseline)
#'   \item lognorm: 0.040 (78% of VST)
#'   \item logcpm: 0.020 (39% of VST)
#' }
#'
#' @examples
#' \dontrun{
#' # Fast normalization for CosMx (recommended)
#' norm_data <- normalize_data(counts, method = "lognorm")
#'
#' # VST for smaller Visium datasets
#' norm_data <- normalize_data(counts, method = "vst")
#'
#' # Custom scale factor
#' norm_data <- normalize_data(counts, method = "lognorm", scale_factor = 1e4)
#' }
#'
#' @seealso \code{\link{load_data_visium}}, \code{\link{load_data_cosmx}}
#'
#' @export
normalize_data <- function(counts,
                           method = c("lognorm", "logcpm", "vst"),
                           scale_factor = NULL,
                           verbose = TRUE) {

    method <- match.arg(method)

    # Convert to matrix if needed
    if (methods::is(counts, "sparseMatrix")) {
        counts <- as.matrix(counts)
    }

    if (verbose) {
        message("Normalizing ", nrow(counts), " genes x ", ncol(counts),
                " cells using ", method)
    }

    # Set default scale factors
    if (is.null(scale_factor)) {
        scale_factor <- switch(method,
            lognorm = 10000,
            logcpm = 1e6,
            vst = NA
        )
    }

    result <- switch(method,
        lognorm = normalize_lognorm(counts, scale_factor, verbose),
        logcpm = normalize_lognorm(counts, scale_factor, verbose),
        vst = normalize_vst(counts, verbose)
    )

    return(result)
}

#' LogNorm Normalization
#'
#' Library size normalization followed by log transformation.
#'
#' @param counts Raw count matrix (genes x cells).
#' @param scale_factor Scale factor (default: 10000).
#' @param verbose Print progress.
#'
#' @return Normalized matrix.
#'
#' @keywords internal
normalize_lognorm <- function(counts, scale_factor = 10000, verbose = TRUE) {
    if (verbose) message("  Computing library sizes...")
    lib_sizes <- colSums(counts)

    if (any(lib_sizes == 0)) {
        warning("Some cells have zero total counts. These will produce NaN values.")
    }

    if (verbose) message("  Normalizing to scale factor ", scale_factor, "...")
    # Normalize: counts / lib_size * scale_factor
    norm_counts <- sweep(counts, 2, lib_sizes, "/") * scale_factor

    if (verbose) message("  Log transforming...")
    result <- log1p(norm_counts)

    if (verbose) message("  Done.")
    return(result)
}

#' VST Normalization
#'
#' Variance stabilizing transformation using sctransform.
#'
#' @param counts Raw count matrix (genes x cells).
#' @param verbose Print progress.
#'
#' @return VST-normalized matrix.
#'
#' @keywords internal
normalize_vst <- function(counts, verbose = TRUE) {
    # Check if sctransform is available
    if (!requireNamespace("sctransform", quietly = TRUE)) {
        stop("VST normalization requires the 'sctransform' package.\n",
             "Install with: install.packages('sctransform')\n",
             "Or use method = 'lognorm' as a faster alternative.")
    }

    if (verbose) message("  Running sctransform VST (this may take a while)...")

    # Run VST
    vst_result <- sctransform::vst(
        umi = counts,
        return_corrected_umi = FALSE,
        return_cell_attr = FALSE,
        verbosity = if (verbose) 1 else 0
    )

    if (verbose) message("  Done.")
    return(vst_result$y)
}

#' Standardize Expression Matrix (Z-score)
#'
#' Gene-wise z-score normalization. Each gene is standardized to have
#' mean = 0 and sd = 1 across cells/spots.
#'
#' @param expr Expression matrix (genes x cells/spots). Can be raw or normalized.
#' @param verbose Logical. Print progress messages (default: FALSE).
#'
#' @return Standardized matrix with same dimensions.
#'
#' @details
#' This function is automatically called internally before computing spatial
#' metrics (Moran's I, I_ND). You typically don't need to call it directly.
#'
#' Formula: \eqn{z_{ij} = (x_{ij} - \bar{x}_i) / \sigma_i}
#'
#' Genes with zero variance (constant expression) are set to zero.
#'
#' @examples
#' \dontrun{
#' # Standardize normalized data
#' expr_z <- standardize_matrix(norm_data)
#' }
#'
#' @export
standardize_matrix <- function(expr, verbose = FALSE) {
    if (verbose) message("Z-score standardizing ", nrow(expr), " genes...")

    # Use C++ implementation if available, otherwise R
    if (exists("standardize_matrix_cpp")) {
        result <- standardize_matrix_cpp(as.matrix(expr))
    } else {
        # R fallback
        result <- t(scale(t(expr)))
        # Handle NaN from zero-variance genes
        result[is.nan(result)] <- 0
    }

    rownames(result) <- rownames(expr)
    colnames(result) <- colnames(expr)

    return(result)
}

#' Standardize a Vector
#'
#' Z-score standardization of a single vector.
#'
#' @param x Numeric vector.
#'
#' @return Standardized vector (mean = 0, sd = 1).
#'
#' @keywords internal
standardize <- function(x) {
    n <- length(x)
    m <- mean(x)
    s <- sqrt(sum((x - m)^2) / n)  # Population SD
    if (s < 1e-10) return(rep(0, n))
    return((x - m) / s)
}
