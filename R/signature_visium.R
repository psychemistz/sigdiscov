#' Compute Spatial Signature for Visium Data
#'
#' Main analysis function for computing spatial correlation signatures in
#' Visium data. Computes Moran's I and/or I_ND across multiple distance radii
#' for all genes against a specified factor gene.
#'
#' @param data A VisiumData list (from \code{\link{load_data_visium}} or
#'   \code{\link{as_data_visium}}).
#' @param factor_gene Character. Name of the factor gene (e.g., "IL1B").
#' @param radii Numeric vector. Distance radii in coordinate units
#'   (default: seq(100, 500, 100)).
#' @param metric Character. Metric to compute: "moran", "ind", or "both"
#'   (default: "both").
#' @param mode Character. Analysis mode: "bivariate" (all spots) or
#'   "directional" (sender/receiver split) (default: "directional").
#' @param sender_percentile Numeric. Percentile threshold for sender selection
#'   (default: 75, meaning top 25 percent of expressors become senders).
#' @param compute_delta Logical. Compute signed delta I (default: TRUE).
#' @param n_perm Integer. Number of permutations for significance testing.
#'   Set to 0 to skip permutation test (default: 0).
#' @param seed Integer. Random seed for reproducibility (optional).
#' @param verbose Logical. Print progress messages (default: TRUE).
#'
#' @return A VisiumSignature data.frame with columns: gene, moran_r1/r2/...,
#'   ind_r1/r2/..., moran_r1_val, moran_max, ind_r1_val, ind_max, delta_I,
#'   delta_I_sign, p_value, p_adj, z_score (depending on metric and options).
#'
#' @details
#' This function is the main entry point for Visium spatial signature analysis.
#'
#' Analysis Modes: (1) Bivariate - All spots participate as both senders and
#' receivers, using a symmetric weight matrix. Appropriate for general
#' spatial autocorrelation analysis. (2) Directional - Spots are split into
#' senders (high factor expression) and receivers (low factor expression),
#' using an asymmetric weight matrix. Appropriate for paracrine signaling analysis.
#'
#' Metrics: (1) Moran's I - Classic spatial autocorrelation measure, unbounded,
#' typically ranges from about -1 to +1. (2) I_ND - Normalized directional
#' Moran's I (cosine similarity), bounded exactly between -1 and +1.
#'
#' @examples
#' \dontrun{
#' # Load data
#' data <- load_data_visium("path/to/spaceranger")
#'
#' # Bivariate analysis
#' sig_bi <- compute_signature_visium(
#'     data,
#'     factor_gene = "IL1B",
#'     metric = "moran",
#'     mode = "bivariate"
#' )
#'
#' # Directional analysis with permutation test
#' sig_dir <- compute_signature_visium(
#'     data,
#'     factor_gene = "IL1B",
#'     metric = "ind",
#'     mode = "directional",
#'     sender_percentile = 75,
#'     n_perm = 999
#' )
#'
#' # Top correlated genes
#' head(sig_dir[order(-sig_dir$ind_r1_val), ], 20)
#'
#' # Significant genes
#' sig_dir[sig_dir$p_adj < 0.05, ]
#' }
#'
#' @seealso
#'   \code{\link{load_data_visium}}, \code{\link{as_data_visium}},
#'   \code{\link{split_by_expression}}, \code{\link{create_weights_visium}}
#'
#' @export
compute_signature_visium <- function(data,
                                      factor_gene,
                                      radii = seq(100, 500, 100),
                                      metric = c("both", "moran", "ind"),
                                      mode = c("directional", "bivariate"),
                                      sender_percentile = 75,
                                      compute_delta = TRUE,
                                      n_perm = 0,
                                      seed = NULL,
                                      verbose = TRUE) {

    metric <- match.arg(metric)
    mode <- match.arg(mode)

    # Validate data
    if (!inherits(data, "VisiumData") && !is.list(data)) {
        stop("data must be a VisiumData object or list with expr, coords, gene_names")
    }

    # Validate factor gene
    if (!factor_gene %in% data$gene_names) {
        stop("Factor gene '", factor_gene, "' not found in expression matrix.\n",
             "Available genes (first 10): ",
             paste(head(data$gene_names, 10), collapse = ", "))
    }

    # Set seed if provided
    if (!is.null(seed)) {
        set.seed(seed)
    }

    # Standardize expression matrix (gene-wise)
    expr_norm <- standardize_matrix(data$expr)

    # Define sender/receiver based on mode
    if (mode == "bivariate") {
        sender_idx <- seq_len(data$n_spots)
        receiver_idx <- seq_len(data$n_spots)
        sender_coords <- data$coords
        receiver_coords <- data$coords

        if (verbose) {
            message("Bivariate mode: all ", data$n_spots, " spots")
        }
    } else {
        # Expression-based split
        factor_expr_raw <- data$expr[factor_gene, ]
        split <- split_by_expression(
            factor_expr_raw,
            percentile = sender_percentile,
            min_spots = 10  # Lower threshold for warnings
        )

        sender_idx <- split$sender_idx
        receiver_idx <- split$receiver_idx
        sender_coords <- data$coords[sender_idx, , drop = FALSE]
        receiver_coords <- data$coords[receiver_idx, , drop = FALSE]

        if (verbose) {
            message("Directional mode:")
            message("  Senders: ", split$n_senders, " spots (expr >= ",
                    round(split$threshold, 3), ")")
            message("  Receivers: ", split$n_receivers, " spots")
        }
    }

    if (verbose) {
        message("Computing Visium signature for ", factor_gene)
        message("  Radii: ", paste(radii, collapse = ", "))
        message("  Metric: ", metric)
    }

    # Extract and standardize factor expression for senders
    factor_expr <- standardize(data$expr[factor_gene, sender_idx])

    # Extract receiver expression (already normalized)
    receiver_expr_norm <- expr_norm[, receiver_idx, drop = FALSE]

    # Initialize result matrices
    n_genes <- data$n_genes
    n_radii <- length(radii)
    gene_names <- data$gene_names

    compute_moran <- metric %in% c("both", "moran")
    compute_ind <- metric %in% c("both", "ind")

    I_moran <- if (compute_moran) matrix(NA_real_, n_genes, n_radii) else NULL
    I_ind <- if (compute_ind) matrix(NA_real_, n_genes, n_radii) else NULL

    # Compute for each radius
    for (r in seq_along(radii)) {
        radius <- radii[r]

        if (verbose) message("  Processing radius ", radius, "...")

        # Build weight matrix
        if (mode == "bivariate") {
            W <- create_weights_visium(data$coords, radius)
        } else {
            W <- create_directional_weights_visium(
                sender_coords, receiver_coords, radius
            )
        }

        # Compute spatial lag for all genes at once: W * Z_g (transpose for genes in columns)
        lag_G <- compute_spatial_lag_batch(W, t(receiver_expr_norm))

        # Compute metrics
        if (compute_moran) {
            I_moran[, r] <- compute_metric_batch(
                factor_expr, lag_G, "moran"
            )
        }
        if (compute_ind) {
            I_ind[, r] <- compute_metric_batch(
                factor_expr, lag_G, "ind"
            )
        }
    }

    # Build result data frame
    result <- data.frame(gene = gene_names, stringsAsFactors = FALSE)

    if (compute_moran) {
        colnames(I_moran) <- paste0("moran_r", seq_along(radii))
        result <- cbind(result, I_moran)
        result$moran_r1_val <- I_moran[, 1]
        result$moran_max <- apply(I_moran, 1, function(x) {
            if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
        })
    }

    if (compute_ind) {
        colnames(I_ind) <- paste0("ind_r", seq_along(radii))
        result <- cbind(result, I_ind)
        result$ind_r1_val <- I_ind[, 1]
        result$ind_max <- apply(I_ind, 1, function(x) {
            if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
        })
    }

    # Compute Delta I
    if (compute_delta && n_radii > 1) {
        I_for_delta <- if (compute_ind) I_ind else I_moran

        delta_results <- apply(I_for_delta, 1, function(curve) {
            compute_delta_i(curve, radii)
        })

        result$delta_I <- sapply(delta_results, `[[`, "delta_I")
        result$delta_I_sign <- sapply(delta_results, `[[`, "sign")
    }

    # Permutation test (at first radius)
    if (n_perm > 0) {
        if (verbose) message("  Running permutation test (", n_perm, " permutations)...")

        # Build weight matrix for first radius
        if (mode == "bivariate") {
            W1 <- create_weights_visium(data$coords, radii[1])
        } else {
            W1 <- create_directional_weights_visium(
                sender_coords, receiver_coords, radii[1]
            )
        }

        perm_metric <- if (compute_ind) "ind" else "moran"

        perm_result <- batch_permutation_test(
            z_f = factor_expr,
            Z_g = t(receiver_expr_norm),
            W = W1,
            metric = perm_metric,
            n_perm = n_perm
        )

        result$p_value <- perm_result$p_value
        result$p_adj <- adjust_pvalues(perm_result$p_value, "BH")
        result$z_score <- perm_result$z_score
    }

    # Return with class and attributes
    structure(
        result,
        class = c("VisiumSignature", "data.frame"),
        factor_gene = factor_gene,
        mode = mode,
        metric = metric,
        radii = radii,
        sender_percentile = sender_percentile,
        n_spots = data$n_spots,
        n_perm = n_perm
    )
}

#' Pairwise Moran's I for Spatial Transcriptomics
#'
#' Compute pairwise Moran's I statistics between all genes in spatial
#' transcriptomics data using Gaussian distance decay weights.
#'
#' @param data A matrix or sparse matrix (dgCMatrix) of gene expression values.
#'   Rows are genes, columns are spots.
#' @param spot_coords A data frame or matrix with spot coordinates. Should have
#'   columns named 'row' and 'col', or be a 2-column matrix.
#' @param max_radius Maximum grid radius for computing distances. Default: 5.
#' @param platform Platform type: "visium" (default) or "old".
#' @param same_spot Whether to consider the same spot in computation. Default: FALSE.
#' @param mode Computation mode: "paired" (all pairs), "first" (pairs with first gene),
#'   or "single" (diagonal only). Default: "paired".
#' @param verbose Print progress messages. Default: TRUE.
#'
#' @return A matrix of pairwise Moran's I values. For mode="paired", returns a
#'   symmetric matrix where element [i,j] is the Moran's I between genes i and j.
#'
#' @details
#' This function computes the classic pairwise Moran's I matrix where each
#' entry (i,j) represents the spatial correlation between gene i and gene j.
#'
#' Uses Gaussian distance decay weights with sigma=100 for Visium data,
#' accounting for the hexagonal grid geometry.
#'
#' @examples
#' \dontrun{
#' # Load VST-transformed data
#' data <- read.table("vst.tsv", header=TRUE, row.names=1)
#' colnames(data) <- gsub("X", "", colnames(data))
#' data <- as.matrix(data)
#'
#' # Parse spot coordinates
#' spot_coords <- parse_spot_names(colnames(data))
#'
#' # Compute pairwise Moran's I
#' result <- pairwise_moran(data, spot_coords, max_radius = 3)
#' result[1:5, 1:5]
#' }
#'
#' @export
pairwise_moran <- function(data,
                           spot_coords,
                           max_radius = 5,
                           platform = c("visium", "old"),
                           same_spot = FALSE,
                           mode = c("paired", "first", "single"),
                           verbose = TRUE) {

    platform <- match.arg(platform)
    mode <- match.arg(mode)

    # Convert platform to integer
    platform_int <- if (platform == "visium") 0L else 1L

    # Convert mode to booleans
    paired_genes <- mode != "single"
    all_genes <- mode == "paired"

    # Parse spot coordinates
    if (is.data.frame(spot_coords)) {
        if ("row" %in% names(spot_coords) && "col" %in% names(spot_coords)) {
            spot_row <- as.integer(spot_coords$row)
            spot_col <- as.integer(spot_coords$col)
        } else {
            spot_row <- as.integer(spot_coords[[1]])
            spot_col <- as.integer(spot_coords[[2]])
        }
    } else if (is.matrix(spot_coords)) {
        spot_row <- as.integer(spot_coords[, 1])
        spot_col <- as.integer(spot_coords[, 2])
    } else {
        stop("spot_coords must be a data frame or matrix")
    }

    # Validate dimensions
    if (ncol(data) != length(spot_row)) {
        stop("Number of columns in data must match number of spots in spot_coords")
    }

    # Convert to matrix if needed and call C++ function
    if (methods::is(data, "sparseMatrix")) {
        # Sparse matrix path
        data <- methods::as(data, "dgCMatrix")
        result <- cpp_compute_moran_full_sparse(
            data,
            spot_row, spot_col,
            as.integer(max_radius),
            platform_int,
            same_spot,
            paired_genes,
            all_genes,
            verbose
        )
    } else {
        # Dense matrix path
        data <- as.matrix(data)
        result <- cpp_compute_moran_full(
            data,
            spot_row, spot_col,
            as.integer(max_radius),
            platform_int,
            same_spot,
            paired_genes,
            all_genes,
            verbose
        )
    }

    # Set row/column names if available
    gene_names <- rownames(data)
    if (!is.null(gene_names)) {
        if (mode == "paired") {
            rownames(result) <- gene_names
            colnames(result) <- gene_names
        } else {
            names(result) <- gene_names
        }
    }

    return(result)
}

#' Compute Moran's I for a Single Gene Pair
#'
#' Compute Moran's I between two specific genes.
#'
#' @param x Numeric vector of expression values for gene 1.
#' @param y Numeric vector of expression values for gene 2.
#' @param W Weight matrix (spots x spots).
#' @param normalize Whether to z-normalize the input vectors. Default: TRUE.
#'
#' @return A single numeric value (Moran's I).
#'
#' @export
moran_I <- function(x, y, W, normalize = TRUE) {
    if (length(x) != length(y)) {
        stop("x and y must have the same length")
    }

    if (nrow(W) != length(x) || ncol(W) != length(x)) {
        stop("W must be a square matrix with dimensions matching x and y")
    }

    # Z-normalize if requested
    if (normalize) {
        x <- (x - mean(x)) / stats::sd(x)
        y <- (y - mean(y)) / stats::sd(y)
    }

    # Compute Moran's I: (x' W y) / sum(W)
    weight_sum <- sum(W)
    if (weight_sum == 0) {
        warning("Weight sum is zero")
        return(0)
    }

    result <- as.numeric(t(x) %*% W %*% y) / weight_sum
    return(result)
}

#' Create Spatial Weight Matrix
#'
#' Create a weight matrix based on spatial distances between spots.
#'
#' @param spot_coords A data frame or matrix with spot coordinates.
#' @param max_radius Maximum grid radius. Default: 5.
#' @param platform Platform type: "visium" or "old". Default: "visium".
#' @param same_spot Whether to include same-spot weights. Default: FALSE.
#'
#' @return A list with components:
#'   \item{W}{The weight matrix (spots x spots)}
#'   \item{weight_sum}{Sum of all weights}
#'
#' @export
create_weight_matrix <- function(spot_coords,
                                 max_radius = 5,
                                 platform = c("visium", "old"),
                                 same_spot = FALSE) {

    platform <- match.arg(platform)
    platform_int <- if (platform == "visium") 0L else 1L

    # Parse spot coordinates
    if (is.data.frame(spot_coords)) {
        if ("row" %in% names(spot_coords) && "col" %in% names(spot_coords)) {
            spot_row <- as.integer(spot_coords$row)
            spot_col <- as.integer(spot_coords$col)
        } else {
            spot_row <- as.integer(spot_coords[[1]])
            spot_col <- as.integer(spot_coords[[2]])
        }
    } else if (is.matrix(spot_coords)) {
        spot_row <- as.integer(spot_coords[, 1])
        spot_col <- as.integer(spot_coords[, 2])
    } else {
        stop("spot_coords must be a data frame or matrix")
    }

    # Create distance lookup
    distance <- cpp_create_distance(as.integer(max_radius), platform_int)
    if (!same_spot) {
        distance[1, 1] <- 0
    }

    # Create weight matrix
    result <- cpp_create_weight_matrix(
        spot_row, spot_col, distance,
        as.integer(max_radius), same_spot
    )

    return(result)
}

#' Print Method for VisiumSignature
#'
#' @param x A VisiumSignature object.
#' @param ... Additional arguments (ignored).
#'
#' @return Invisibly returns the input object.
#'
#' @export
print.VisiumSignature <- function(x, ...) {
    cat("VisiumSignature object\n")
    cat("  Factor gene:", attr(x, "factor_gene"), "\n")
    cat("  Mode:", attr(x, "mode"), "\n")
    cat("  Metric:", attr(x, "metric"), "\n")
    cat("  Radii:", paste(attr(x, "radii"), collapse = ", "), "\n")
    cat("  Genes tested:", nrow(x), "\n")

    if ("p_value" %in% names(x)) {
        n_sig <- sum(x$p_adj < 0.05, na.rm = TRUE)
        cat("  Significant (p_adj < 0.05):", n_sig, "\n")
    }

    cat("\nTop 5 genes by first radius metric:\n")
    metric_col <- if ("ind_r1_val" %in% names(x)) "ind_r1_val" else "moran_r1_val"
    top_idx <- order(x[[metric_col]], decreasing = TRUE)[1:min(5, nrow(x))]
    print(x[top_idx, c("gene", metric_col)], row.names = FALSE)

    invisible(x)
}
