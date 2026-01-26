#' Compute Cell Type Pair Specific Spatial Analysis
#'
#' Main entry point for computing spatial correlation metrics (Moran's I, I_ND,
#' Delta I) across cell type sender-receiver pairs at multiple radii. Designed
#' for large single-cell spatial datasets with streaming HDF5 output.
#'
#' @param data SCData list with expr, coords, cell_types, gene_names
#' @param radii Numeric vector of distance radii (default: seq(20, 500, by=20))
#' @param pairs Character matrix or NULL. Cell type pairs to analyze.
#'   If NULL (default), analyzes all pairs. Format: 2-column matrix with
#'   "sender" and "receiver" columns.
#' @param min_cells Minimum cells per type to include pair (default: 10)
#' @param output_dir Directory for HDF5 output. If NULL, returns results in memory.
#' @param compute_moran Compute pairwise Moran's I (default: TRUE)
#' @param compute_ind Compute I_ND factor-gene metric (default: FALSE)
#' @param factor_genes Character vector of factor genes for I_ND
#' @param method Computation method: "streaming" (default) or "matrix"
#' @param sigma_factor Sigma = radius * sigma_factor (default: 1/3)
#' @param verbose Print progress messages (default: TRUE)
#'
#' @return If output_dir is NULL, returns list with:
#' \describe{
#'   \item{delta_I}{List of delta_I matrices per pair}
#'   \item{I_max}{List of I_max matrices per pair}
#'   \item{pair_stats}{Data frame with pair statistics}
#' }
#' If output_dir is provided, saves to HDF5 and returns path.
#'
#' @details
#' Workflow for each cell type pair:
#' \enumerate{
#'   \item Extract sender and receiver cells
#'   \item For each radius, compute Moran's I matrix (genes x genes)
#'   \item Compute Delta I from I(r) curves
#'   \item Stream results to HDF5 or accumulate in memory
#' }
#'
#' Memory usage: ~100 MB per pair (vs 200 GB if all pairs in memory)
#'
#' @examples
#' \dontrun{
#' # Load CosMx data
#' data <- load_data_cosmx("expr.csv", "meta.csv")
#'
#' # Compute for all cell type pairs
#' result <- compute_celltype_pair_analysis(data, radii = seq(50, 500, by=50))
#'
#' # Specific pairs only
#' pairs <- matrix(c("Macrophage", "Fibroblast",
#'                   "Macrophage", "Epithelial"), ncol = 2, byrow = TRUE)
#' result <- compute_celltype_pair_analysis(data, pairs = pairs)
#'
#' # With HDF5 output for large analysis
#' compute_celltype_pair_analysis(data, output_dir = "results/")
#' }
#'
#' @export
compute_celltype_pair_analysis <- function(
    data,
    radii = seq(20, 500, by = 20),
    pairs = NULL,
    min_cells = 10,
    output_dir = NULL,
    compute_moran = TRUE,
    compute_ind = FALSE,
    factor_genes = NULL,
    method = c("streaming", "matrix"),
    sigma_factor = 1/3,
    verbose = TRUE
) {
    method <- match.arg(method)

    # Validate inputs
    if (!inherits(data, "SCData") && !is.list(data)) {
        stop("data must be an SCData object or list with expr, coords, cell_types")
    }

    required_fields <- c("expr", "coords", "cell_types", "gene_names")
    missing <- setdiff(required_fields, names(data))
    if (length(missing) > 0) {
        stop("data is missing required fields: ", paste(missing, collapse = ", "))
    }

    if (compute_ind && is.null(factor_genes)) {
        stop("factor_genes must be provided when compute_ind = TRUE")
    }

    # Get unique cell types
    cell_types_unique <- sort(unique(data$cell_types))
    n_types <- length(cell_types_unique)

    if (verbose) {
        message("Cell Type Pair Analysis")
        message("  ", data$n_genes, " genes, ", data$n_cells, " cells")
        message("  ", n_types, " cell types")
        message("  ", length(radii), " radii: ", radii[1], " to ", radii[length(radii)])
    }

    # Generate pairs if not provided
    if (is.null(pairs)) {
        pairs <- expand.grid(
            sender = cell_types_unique,
            receiver = cell_types_unique,
            stringsAsFactors = FALSE
        )
        pairs <- as.matrix(pairs)
    }
    n_pairs <- nrow(pairs)

    if (verbose) {
        message("  ", n_pairs, " cell type pairs to analyze")
    }

    # Setup HDF5 output if requested
    use_hdf5 <- !is.null(output_dir)
    if (use_hdf5) {
        if (!requireNamespace("rhdf5", quietly = TRUE)) {
            stop("rhdf5 package required for HDF5 output. Install with:\n",
                 "BiocManager::install('rhdf5')")
        }
        if (!dir.exists(output_dir)) {
            dir.create(output_dir, recursive = TRUE)
        }
        h5_path <- file.path(output_dir, "celltype_pair_analysis.h5")
        .init_hdf5_output(h5_path, n_pairs, data$n_genes, length(radii),
                          pairs, radii, data$gene_names)
    }

    # Results storage (if not using HDF5)
    results <- list(
        delta_I = list(),
        I_max = list(),
        I_curves = list()
    )

    pair_stats <- data.frame(
        pair_idx = integer(),
        sender = character(),
        receiver = character(),
        n_sender = integer(),
        n_receiver = integer(),
        mean_delta_I = numeric(),
        max_delta_I = numeric(),
        stringsAsFactors = FALSE
    )

    # Process each pair
    for (p in seq_len(n_pairs)) {
        sender_type <- pairs[p, 1]
        receiver_type <- pairs[p, 2]
        pair_name <- paste0(sender_type, "->", receiver_type)

        if (verbose) {
            message("\nPair ", p, "/", n_pairs, ": ", pair_name)
        }

        # Get cell indices
        sender_idx <- which(data$cell_types == sender_type)
        receiver_idx <- which(data$cell_types == receiver_type)

        n_sender <- length(sender_idx)
        n_receiver <- length(receiver_idx)

        if (verbose) {
            message("  Senders: ", n_sender, ", Receivers: ", n_receiver)
        }

        # Skip if too few cells
        if (n_sender < min_cells || n_receiver < min_cells) {
            if (verbose) {
                message("  Skipping: insufficient cells (min = ", min_cells, ")")
            }
            next
        }

        # Extract data subsets
        sender_data <- data$expr[, sender_idx, drop = FALSE]
        receiver_data <- data$expr[, receiver_idx, drop = FALSE]
        sender_coords <- data$coords[sender_idx, , drop = FALSE]
        receiver_coords <- data$coords[receiver_idx, , drop = FALSE]

        # Compute Moran's I for all radii
        if (compute_moran) {
            pair_result <- compute_celltype_pair_moran_cpp(
                sender_data = sender_data,
                receiver_data = receiver_data,
                sender_coords = as.matrix(sender_coords),
                receiver_coords = as.matrix(receiver_coords),
                radii = radii,
                sigma_factor = sigma_factor,
                verbose = verbose
            )

            # Store results
            if (use_hdf5) {
                .write_pair_to_hdf5(h5_path, p, pair_result)
            } else {
                results$delta_I[[pair_name]] <- pair_result$delta_I
                results$I_max[[pair_name]] <- pair_result$I_max
                results$I_curves[[pair_name]] <- pair_result$I_curves
            }

            # Compute pair statistics
            delta_vals <- as.vector(pair_result$delta_I)
            delta_vals <- delta_vals[is.finite(delta_vals)]

            pair_stats <- rbind(pair_stats, data.frame(
                pair_idx = p,
                sender = sender_type,
                receiver = receiver_type,
                n_sender = n_sender,
                n_receiver = n_receiver,
                mean_delta_I = mean(abs(delta_vals), na.rm = TRUE),
                max_delta_I = max(abs(delta_vals), na.rm = TRUE),
                stringsAsFactors = FALSE
            ))
        }

        # Force garbage collection to free memory
        gc(verbose = FALSE)
    }

    # Finalize output
    if (use_hdf5) {
        if (verbose) message("\nResults saved to: ", h5_path)

        # Save summary RDS
        summary_path <- file.path(output_dir, "pair_stats.rds")
        saveRDS(pair_stats, summary_path)

        return(invisible(list(
            h5_path = h5_path,
            pair_stats = pair_stats
        )))
    } else {
        results$pair_stats <- pair_stats
        results$gene_names <- data$gene_names
        results$radii <- radii
        results$pairs <- pairs

        class(results) <- c("CellTypePairResult", "list")
        return(results)
    }
}


#' Compute Directional Pairwise Moran's I for Cell Type Pair
#'
#' Lower-level function to compute Moran's I for a single cell type pair
#' at a single radius.
#'
#' @param sender_data Gene expression for sender cells (genes x n_sender)
#' @param receiver_data Gene expression for receiver cells (genes x n_receiver)
#' @param sender_coords Sender cell coordinates (n_sender x 2)
#' @param receiver_coords Receiver cell coordinates (n_receiver x 2)
#' @param radius Distance radius
#' @param sigma Gaussian sigma (default: radius/3)
#' @param method Computation method: "streaming" or "matrix"
#' @param normalize Normalize data before computation (default: TRUE)
#' @param verbose Print progress messages (default: FALSE)
#'
#' @return List with moran matrix, weight_sum, n_edges
#'
#' @export
pairwise_moran_celltype_pair <- function(
    sender_data,
    receiver_data,
    sender_coords,
    receiver_coords,
    radius,
    sigma = NULL,
    method = c("streaming", "matrix"),
    normalize = TRUE,
    verbose = FALSE
) {
    method <- match.arg(method)

    sender_data <- as.matrix(sender_data)
    receiver_data <- as.matrix(receiver_data)
    sender_coords <- as.matrix(sender_coords)
    receiver_coords <- as.matrix(receiver_coords)

    if (is.null(sigma)) {
        sigma <- radius / 3
    }

    if (method == "streaming") {
        result <- pairwise_moran_directional_streaming_cpp(
            sender_data = sender_data,
            receiver_data = receiver_data,
            sender_coords = sender_coords,
            receiver_coords = receiver_coords,
            radius = radius,
            sigma = sigma,
            normalize_data = normalize,
            verbose = verbose
        )
    } else {
        # Matrix method: build explicit weight matrix
        W_result <- create_directional_weights_sc_cpp(
            sender_coords = sender_coords,
            receiver_coords = receiver_coords,
            radius = radius,
            sigma = sigma,
            use_kdtree = TRUE,
            verbose = verbose
        )

        # Normalize data if requested
        if (normalize) {
            sender_data <- standardize_matrix(sender_data)
            receiver_data <- standardize_matrix(receiver_data)
        }

        result <- pairwise_moran_directional_cpp(
            sender_data = sender_data,
            receiver_data = receiver_data,
            W = W_result$W,
            verbose = verbose
        )
        result$n_edges <- W_result$n_edges
    }

    # Add gene names if available
    if (!is.null(rownames(sender_data))) {
        rownames(result$moran) <- rownames(sender_data)
        colnames(result$moran) <- rownames(sender_data)
    }

    result
}


#' Create Directional Weight Matrix for Cell Type Pair
#'
#' Creates a sparse weight matrix from sender to receiver cells.
#' Supports annular (ring) weights via inner_radius parameter.
#'
#' @param sender_coords Sender cell coordinates (n_sender x 2)
#' @param receiver_coords Receiver cell coordinates (n_receiver x 2)
#' @param radius Distance radius (outer radius)
#' @param inner_radius Inner radius for annular weights (default: 0, circular)
#' @param sigma Gaussian sigma (default: radius/3)
#' @param verbose Print progress messages (default: FALSE)
#'
#' @return List with W (sparse matrix), weight_sum, n_edges
#'
#' @export
create_directional_weights <- function(
    sender_coords,
    receiver_coords,
    radius,
    inner_radius = 0,
    sigma = NULL,
    verbose = FALSE
) {
    sender_coords <- as.matrix(sender_coords)
    receiver_coords <- as.matrix(receiver_coords)

    if (is.null(sigma)) {
        sigma <- radius / 3
    }

    create_directional_weights_sc_cpp(
        sender_coords = sender_coords,
        receiver_coords = receiver_coords,
        radius = radius,
        inner_radius = inner_radius,
        sigma = sigma,
        use_kdtree = nrow(receiver_coords) > 1000,
        verbose = verbose
    )
}


#' Print Method for CellTypePairResult
#'
#' @param x CellTypePairResult object
#' @param ... Additional arguments (ignored)
#' @return Invisibly returns the input object
#' @export
print.CellTypePairResult <- function(x, ...) {
    cat("CellTypePairResult object\n")
    cat("  Genes:", length(x$gene_names), "\n")
    cat("  Radii:", length(x$radii), "(", x$radii[1], "-", x$radii[length(x$radii)], ")\n")
    cat("  Pairs analyzed:", nrow(x$pair_stats), "\n")

    if (nrow(x$pair_stats) > 0) {
        cat("\nTop 5 pairs by max delta I:\n")
        top_idx <- order(-x$pair_stats$max_delta_I)[1:min(5, nrow(x$pair_stats))]
        print(x$pair_stats[top_idx, c("sender", "receiver", "n_sender",
                                       "n_receiver", "max_delta_I")],
              row.names = FALSE)
    }

    invisible(x)
}


# =============================================================================
# HDF5 I/O Helper Functions
# =============================================================================

#' Initialize HDF5 Output File
#'
#' @keywords internal
.init_hdf5_output <- function(h5_path, n_pairs, n_genes, n_radii,
                               pairs, radii, gene_names) {
    # Delete existing file
    if (file.exists(h5_path)) {
        file.remove(h5_path)
    }

    rhdf5::h5createFile(h5_path)

    # Create groups
    rhdf5::h5createGroup(h5_path, "moran")
    rhdf5::h5createGroup(h5_path, "metadata")

    # Create datasets with chunking for efficient streaming
    # Delta I: pairs x genes x genes
    rhdf5::h5createDataset(
        h5_path, "moran/delta_I",
        dims = c(n_pairs, n_genes, n_genes),
        chunk = c(1, n_genes, n_genes),  # Chunk by pair
        level = 6
    )

    # I_max: pairs x genes x genes
    rhdf5::h5createDataset(
        h5_path, "moran/I_max",
        dims = c(n_pairs, n_genes, n_genes),
        chunk = c(1, n_genes, n_genes),
        level = 6
    )

    # Argmax: pairs x genes x genes
    rhdf5::h5createDataset(
        h5_path, "moran/argmax",
        dims = c(n_pairs, n_genes, n_genes),
        chunk = c(1, n_genes, n_genes),
        level = 6
    )

    # I_curves: pairs x genes x genes x radii (optional, large)
    # Commented out by default to save space
    # rhdf5::h5createDataset(
    #     h5_path, "moran/I_curves",
    #     dims = c(n_pairs, n_genes, n_genes, n_radii),
    #     chunk = c(1, n_genes, n_genes, n_radii),
    #     level = 6
    # )

    # Metadata
    rhdf5::h5write(pairs[, 1], h5_path, "metadata/sender_types")
    rhdf5::h5write(pairs[, 2], h5_path, "metadata/receiver_types")
    rhdf5::h5write(radii, h5_path, "metadata/radii")
    rhdf5::h5write(gene_names, h5_path, "metadata/gene_names")

    invisible(h5_path)
}


#' Write Pair Results to HDF5
#'
#' @keywords internal
.write_pair_to_hdf5 <- function(h5_path, pair_idx, result) {
    # Write delta_I
    rhdf5::h5write(
        result$delta_I,
        h5_path,
        "moran/delta_I",
        index = list(pair_idx, NULL, NULL)
    )

    # Write I_max
    rhdf5::h5write(
        result$I_max,
        h5_path,
        "moran/I_max",
        index = list(pair_idx, NULL, NULL)
    )

    # Write argmax
    rhdf5::h5write(
        result$argmax,
        h5_path,
        "moran/argmax",
        index = list(pair_idx, NULL, NULL)
    )

    invisible(NULL)
}


#' Read Cell Type Pair Results from HDF5
#'
#' @param h5_path Path to HDF5 file
#' @param pair_idx Pair index or indices to read (default: all)
#' @param what What to read: "delta_I", "I_max", "argmax", or "all"
#'
#' @return List with requested data
#'
#' @export
read_celltype_pair_hdf5 <- function(h5_path, pair_idx = NULL, what = "all") {
    if (!requireNamespace("rhdf5", quietly = TRUE)) {
        stop("rhdf5 package required")
    }

    # Read metadata
    gene_names <- rhdf5::h5read(h5_path, "metadata/gene_names")
    radii <- rhdf5::h5read(h5_path, "metadata/radii")
    sender_types <- rhdf5::h5read(h5_path, "metadata/sender_types")
    receiver_types <- rhdf5::h5read(h5_path, "metadata/receiver_types")

    result <- list(
        gene_names = gene_names,
        radii = radii,
        pairs = cbind(sender = sender_types, receiver = receiver_types)
    )

    # Read data
    if (is.null(pair_idx)) {
        # Read all
        if (what %in% c("all", "delta_I")) {
            result$delta_I <- rhdf5::h5read(h5_path, "moran/delta_I")
        }
        if (what %in% c("all", "I_max")) {
            result$I_max <- rhdf5::h5read(h5_path, "moran/I_max")
        }
        if (what %in% c("all", "argmax")) {
            result$argmax <- rhdf5::h5read(h5_path, "moran/argmax")
        }
    } else {
        # Read specific pairs
        if (what %in% c("all", "delta_I")) {
            result$delta_I <- rhdf5::h5read(h5_path, "moran/delta_I",
                                            index = list(pair_idx, NULL, NULL))
        }
        if (what %in% c("all", "I_max")) {
            result$I_max <- rhdf5::h5read(h5_path, "moran/I_max",
                                          index = list(pair_idx, NULL, NULL))
        }
        if (what %in% c("all", "argmax")) {
            result$argmax <- rhdf5::h5read(h5_path, "moran/argmax",
                                           index = list(pair_idx, NULL, NULL))
        }
        result$pairs <- result$pairs[pair_idx, , drop = FALSE]
    }

    result
}


#' Extract Top Gene Pairs by Delta I
#'
#' @param result CellTypePairResult object or HDF5 path
#' @param n_top Number of top pairs to return per cell type pair
#' @param min_delta_I Minimum absolute delta I threshold
#'
#' @return Data frame with top gene pairs
#'
#' @export
extract_top_delta_i <- function(result, n_top = 100, min_delta_I = 0.1) {
    if (is.character(result)) {
        result <- read_celltype_pair_hdf5(result)
    }

    top_pairs <- data.frame()

    pair_names <- if (!is.null(result$pair_stats)) {
        paste0(result$pair_stats$sender, "->", result$pair_stats$receiver)
    } else {
        paste0(result$pairs[, 1], "->", result$pairs[, 2])
    }

    for (i in seq_along(result$delta_I)) {
        delta_mat <- result$delta_I[[i]]
        if (is.null(delta_mat)) next

        pair_name <- if (is.list(result$delta_I)) {
            names(result$delta_I)[i]
        } else {
            pair_names[i]
        }

        # Extract upper triangle (avoid duplicates for symmetric)
        n_genes <- nrow(delta_mat)
        gene_names <- result$gene_names

        for (g1 in seq_len(n_genes)) {
            for (g2 in seq_len(n_genes)) {
                delta_val <- delta_mat[g1, g2]
                if (is.finite(delta_val) && abs(delta_val) >= min_delta_I) {
                    top_pairs <- rbind(top_pairs, data.frame(
                        pair = pair_name,
                        gene1 = gene_names[g1],
                        gene2 = gene_names[g2],
                        delta_I = delta_val,
                        stringsAsFactors = FALSE
                    ))
                }
            }
        }
    }

    # Sort and limit
    if (nrow(top_pairs) > 0) {
        top_pairs <- top_pairs[order(-abs(top_pairs$delta_I)), ]
        if (nrow(top_pairs) > n_top) {
            top_pairs <- top_pairs[seq_len(n_top), ]
        }
    }

    top_pairs
}


#' Apply FDR Correction to Spatial Analysis Results
#'
#' Applies Benjamini-Hochberg FDR correction to p-values from
#' permutation testing or other significance tests.
#'
#' @param results Data frame with 'pvalue' column, or CellTypePairResult object
#' @param alpha FDR threshold (default: 0.05)
#' @param method P-value adjustment method (default: "BH" for Benjamini-Hochberg)
#'
#' @return Results with added columns: p_adj (adjusted p-values), significant (logical)
#'
#' @details
#' Supports multiple input types:
#' \itemize{
#'   \item CellTypePairResult objects with a pvalue component
#'   \item Data frames with a 'pvalue' column
#' }
#'
#' Available methods (passed to \code{stats::p.adjust}):
#' \itemize{
#'   \item "BH" or "fdr": Benjamini-Hochberg (default, controls FDR)
#'   \item "bonferroni": Bonferroni (controls FWER, conservative)
#'   \item "holm": Holm-Bonferroni (controls FWER)
#'   \item "BY": Benjamini-Yekutieli (controls FDR under dependence)
#' }
#'
#' @examples
#' \dontrun{
#' # Correct p-values from permutation test
#' results <- data.frame(
#'     gene = paste0("Gene", 1:100),
#'     pvalue = runif(100, 0, 1)
#' )
#' corrected <- apply_fdr_correction(results, alpha = 0.05)
#'
#' # View significant results
#' significant_genes <- corrected[corrected$significant, ]
#' }
#'
#' @export
apply_fdr_correction <- function(results, alpha = 0.05, method = "BH") {
    if (inherits(results, "CellTypePairResult")) {
        # Handle CellTypePairResult object
        if (is.null(results$pvalue)) {
            warning("No p-values found in results. Run permutation test first.")
            return(results)
        }
        results$p_adj <- stats::p.adjust(results$pvalue, method = method)
        results$significant <- results$p_adj < alpha
    } else if (is.data.frame(results)) {
        # Handle data frame
        if (!"pvalue" %in% names(results)) {
            stop("Results must have 'pvalue' column")
        }
        results$p_adj <- stats::p.adjust(results$pvalue, method = method)
        results$significant <- results$p_adj < alpha
    } else {
        stop("results must be a CellTypePairResult or data.frame with 'pvalue' column")
    }

    n_total <- if (is.data.frame(results)) nrow(results) else length(results$pvalue)
    n_sig <- sum(results$significant, na.rm = TRUE)
    message(sprintf("FDR correction (method=%s, alpha=%.2f): %d/%d significant",
                    method, alpha, n_sig, n_total))

    results
}
