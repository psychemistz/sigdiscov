#' Compute Spatial Signature for Single-Cell Data
#'
#' Main analysis function for computing spatial correlation signatures in
#' single-cell spatial transcriptomics data. Matches Python genomewide_interaction
#' implementation.
#'
#' @param data SCData list (from load_data_cosmx, load_data_anndata, or as_data_sc)
#' @param factor_gene Character. Name of the factor gene (e.g., "IL1B")
#' @param sender_type Character. Sender cell type name
#' @param receiver_type Character. Receiver cell type name
#' @param radii Numeric vector. Distance radii in micrometers
#'   (default: c(10, 20, 30, 50, 100, 200, 300, 500))
#' @param metric Character. Metric to compute: "ind" or "moran" (default: "ind")
#' @param quantile_prob Numeric. Quantile for sender expression filter
#'   (default: 0.25, meaning keep top 75 percent of sender-type cells)
#' @param compute_delta Logical. Compute delta I (default: TRUE)
#' @param n_perm Integer. Number of permutations (default: 0)
#' @param min_cells Integer. Minimum cells (default: 10)
#' @param min_connections Integer. Minimum neighbors (default: 10)
#' @param seed Integer. Random seed (optional)
#' @param verbose Logical. Print progress (default: TRUE)
#'
#' @return ScSignature data.frame with columns: gene, I_r1/r2/..., I_r1 (first radius),
#'   I_max, delta_I, p_value, p_adj, z_score (depending on options)
#'
#' @details
#' This function implements the single-cell spatial analysis workflow:
#'
#' (1) Sender Definition: Cells of sender_type with factor expression above
#' the quantile threshold (default: top 75 percent).
#'
#' (2) Receiver Definition: All cells of receiver_type.
#'
#' (3) Weight Matrix: Gaussian weights with sigma = radius/3, row-normalized.
#'
#' (4) Metric: I_ND (cosine similarity) or Moran's I.
#'
#' @examples
#' \dontrun{
#' # Load data
#' data <- load_data_cosmx("expression.csv", "metadata.csv")
#'
#' # Compute signature
#' sig <- compute_signature_sc(
#'     data,
#'     factor_gene = "IL1B",
#'     sender_type = "Macrophage",
#'     receiver_type = "Fibroblast",
#'     n_perm = 999
#' )
#'
#' # Top responding genes
#' head(sig[order(-sig$I_r1), ], 20)
#' }
#'
#' @seealso
#'   \code{\link{load_data_cosmx}}, \code{\link{define_sender_receiver_sc}},
#'   \code{\link{create_weights_sc}}
#'
#' @export
compute_signature_sc <- function(data,
                                  factor_gene,
                                  sender_type,
                                  receiver_type,
                                  radii = c(10, 20, 30, 50, 100, 200, 300, 500),
                                  metric = c("ind", "moran"),
                                  quantile_prob = 0.25,
                                  compute_delta = TRUE,
                                  n_perm = 0,
                                  min_cells = 10,
                                  min_connections = 10,
                                  seed = NULL,
                                  verbose = TRUE) {

    metric <- match.arg(metric)

    # Validate inputs
    if (!inherits(data, "SCData") && !is.list(data)) {
        stop("data must be an SCData object or list with expr, coords, cell_types")
    }
    if (!factor_gene %in% data$gene_names) {
        stop("Factor gene '", factor_gene, "' not found in expression matrix.\n",
             "Available genes (first 10): ",
             paste(head(data$gene_names, 10), collapse = ", "))
    }
    if (!sender_type %in% data$cell_types) {
        stop("Sender type '", sender_type, "' not found in cell types.\n",
             "Available types: ", paste(unique(data$cell_types), collapse = ", "))
    }
    if (!receiver_type %in% data$cell_types) {
        stop("Receiver type '", receiver_type, "' not found in cell types.\n",
             "Available types: ", paste(unique(data$cell_types), collapse = ", "))
    }

    # Set seed if provided
    if (!is.null(seed)) {
        set.seed(seed)
    }

    # Define sender/receiver (using RAW expression for threshold)
    sr <- define_sender_receiver_sc(
        cell_types = data$cell_types,
        sender_type = sender_type,
        receiver_type = receiver_type,
        factor_expr = data$expr[factor_gene, ],  # RAW
        quantile_prob = quantile_prob,
        min_cells = min_cells
    )

    if (verbose) {
        message("Computing single-cell signature:")
        message("  Factor: ", factor_gene)
        message("  Sender: ", sender_type, " (", sr$n_senders, " high-expr of ",
                sr$n_sender_type, " total)")
        message("  Receiver: ", receiver_type, " (", sr$n_receivers, ")")
        message("  Radii: ", paste(radii, collapse = ", "), " um")
        message("  Expression threshold: ", round(sr$threshold, 3))
    }

    # Standardize expression (gene-wise)
    expr_norm <- standardize_matrix(data$expr)

    # Extract coordinates and expression
    sender_coords <- data$coords[sr$sender_idx, , drop = FALSE]
    receiver_coords <- data$coords[sr$receiver_idx, , drop = FALSE]

    factor_expr <- expr_norm[factor_gene, sr$sender_idx]
    receiver_expr_norm <- expr_norm[, sr$receiver_idx, drop = FALSE]

    # Compute for each radius
    n_genes <- data$n_genes
    n_radii <- length(radii)
    I_matrix <- matrix(NA_real_, n_genes, n_radii)

    for (r in seq_along(radii)) {
        radius <- radii[r]
        inner_radius <- if (r == 1) 0 else radii[r - 1]

        if (verbose) message("  Radius ", radius, " um...")

        # Build Gaussian weight matrix
        W <- create_weights_sc(
            sender_coords, receiver_coords,
            radius = radius,
            inner_radius = inner_radius
        )

        if (length(W@x) < min_connections) {
            warning("Insufficient connections at radius ", radius)
            next
        }

        # Compute spatial lag for all genes
        lag_G <- compute_spatial_lag_batch(W, t(receiver_expr_norm))

        # Compute metric
        I_matrix[, r] <- compute_metric_batch(
            factor_expr, lag_G, metric
        )
    }

    # Build result
    result <- data.frame(gene = data$gene_names, stringsAsFactors = FALSE)
    colnames(I_matrix) <- paste0("I_", radii)
    result <- cbind(result, I_matrix)
    result$I_r1 <- I_matrix[, 1]
    result$I_max <- apply(I_matrix, 1, function(x) {
        if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
    })

    # Delta I
    if (compute_delta && n_radii > 1) {
        delta_results <- apply(I_matrix, 1, function(curve) {
            compute_delta_i(curve, radii)
        })
        result$delta_I <- sapply(delta_results, `[[`, "delta_I")
        result$delta_I_sign <- sapply(delta_results, `[[`, "sign")
    }

    # Permutation test
    if (n_perm > 0) {
        if (verbose) message("  Permutation test (", n_perm, ")...")

        W1 <- create_weights_sc(sender_coords, receiver_coords, radii[1])

        perm_result <- batch_permutation_test(
            z_f = factor_expr,
            Z_g = t(receiver_expr_norm),
            W = W1,
            metric = metric,
            n_perm = n_perm
        )

        result$p_value <- perm_result$p_value
        result$p_adj <- adjust_pvalues(perm_result$p_value, "BH")
        result$z_score <- perm_result$z_score
    }

    structure(
        result,
        class = c("ScSignature", "data.frame"),
        factor_gene = factor_gene,
        sender_type = sender_type,
        receiver_type = receiver_type,
        metric = metric,
        radii = radii,
        quantile_prob = quantile_prob,
        n_senders = sr$n_senders,
        n_receivers = sr$n_receivers
    )
}

#' Print Method for ScSignature
#'
#' @param x An ScSignature object.
#' @param ... Additional arguments (ignored).
#' @return Invisibly returns the input object.
#' @export
print.ScSignature <- function(x, ...) {
    cat("ScSignature object (Single-Cell Spatial)\n")
    cat("  Factor gene:", attr(x, "factor_gene"), "\n")
    cat("  Sender:", attr(x, "sender_type"), "(", attr(x, "n_senders"), "cells )\n")
    cat("  Receiver:", attr(x, "receiver_type"), "(", attr(x, "n_receivers"), "cells )\n")
    cat("  Metric:", attr(x, "metric"), "\n")
    cat("  Radii:", paste(attr(x, "radii"), collapse = ", "), "um\n")
    cat("  Genes tested:", nrow(x), "\n")

    if ("p_value" %in% names(x)) {
        n_sig <- sum(x$p_adj < 0.05, na.rm = TRUE)
        cat("  Significant (p_adj < 0.05):", n_sig, "\n")
    }

    cat("\nTop 5 genes by first radius metric:\n")
    top_idx <- order(x$I_r1, decreasing = TRUE)[1:min(5, nrow(x))]
    top_df <- as.data.frame(x[top_idx, c("gene", "I_r1")])
    print(top_df, row.names = FALSE)

    invisible(x)
}
