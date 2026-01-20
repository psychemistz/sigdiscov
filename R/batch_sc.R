#' Compute I_ND for All Cell Type Pairs
#'
#' Batch computation of I_ND for all cell type pairs across multiple radii.
#' Matches Python genomewide_interaction output format (4D array).
#'
#' @param data SCData list (from load_data_cosmx or similar)
#' @param factor_genes Character vector. Factor gene names (NULL = all genes)
#' @param radii Numeric vector. Distance radii in micrometers
#'   (default: c(10, 20, 30, 50, 100, 200, 300, 500))
#' @param quantile_prob Numeric. Quantile threshold (default: 0.25 = top 75 percent)
#' @param exclude_self Logical. Exclude same-type pairs (default: TRUE)
#' @param min_cells Integer. Minimum cells per type (default: 10)
#' @param min_connections Integer. Minimum neighbors (default: 10)
#' @param output_file Character. Optional HDF5 output path (requires rhdf5)
#' @param verbose Logical. Print progress (default: TRUE)
#'
#' @return List with components:
#'   \itemize{
#'     \item IND_array: 4D array (pairs x radii x factors x genes)
#'     \item pairs: Character vector of pair names ("sender->receiver")
#'     \item radii: Numeric vector of radii
#'     \item factor_genes: Character vector of factor gene names
#'     \item target_genes: Character vector of all gene names
#'     \item cell_types: Character vector of valid cell types
#'   }
#'
#' @details
#' This function computes the full I_ND matrix for all cell type pairs,
#' producing output compatible with Python analysis pipelines. The 4D array
#' can be indexed as: \code{IND_array[pair, radius, factor, gene]}
#'
#' @examples
#' \dontrun{
#' # Load data
#' data <- load_data_cosmx("expression.csv", "metadata.csv")
#'
#' # Compute for specific factors
#' result <- compute_batch_sc(
#'     data,
#'     factor_genes = c("IL1B", "IL6", "TGFB1", "IFNG"),
#'     output_file = "IND_matrix.h5"
#' )
#'
#' # Access specific result
#' # result$IND_array[pair_idx, radius_idx, factor_idx, gene_idx]
#' }
#'
#' @seealso
#'   \code{\link{compute_signature_sc}}, \code{\link{save_hdf5_sc}},
#'   \code{\link{load_hdf5_sc}}
#'
#' @export
compute_batch_sc <- function(data,
                              factor_genes = NULL,
                              radii = c(10, 20, 30, 50, 100, 200, 300, 500),
                              quantile_prob = 0.25,
                              exclude_self = TRUE,
                              min_cells = 10,
                              min_connections = 10,
                              output_file = NULL,
                              verbose = TRUE) {

    # Validate data
    if (!inherits(data, "SCData") && !is.list(data)) {
        stop("data must be an SCData object or list with expr, coords, cell_types")
    }

    # Determine factors
    if (is.null(factor_genes)) {
        factor_genes <- data$gene_names
    }
    factor_idx <- match(factor_genes, data$gene_names)
    valid <- !is.na(factor_idx)
    if (sum(valid) == 0) {
        stop("No valid factor genes found in expression matrix")
    }
    factor_idx <- factor_idx[valid]
    factor_genes <- factor_genes[valid]

    # Get valid cell type pairs
    pairs <- get_celltype_pairs(data$cell_types, min_cells, exclude_self)

    n_pairs <- nrow(pairs)
    n_radii <- length(radii)
    n_factors <- length(factor_genes)
    n_genes <- data$n_genes

    if (verbose) {
        message("Batch computation (matching Python format):")
        message("  Cell types: ", length(unique(c(pairs$sender, pairs$receiver))))
        message("  Pairs: ", n_pairs)
        message("  Radii: ", n_radii)
        message("  Factors: ", n_factors)
        message("  Genes: ", n_genes)
    }

    # Call C++ batch function
    result_list <- batch_compute_all_pairs_cpp(
        expr_matrix = data$expr,
        coords = data$coords,
        cell_types = data$cell_types,
        factor_indices = factor_idx - 1L,  # Convert to 0-based
        radii = radii,
        pairs = pairs,
        quantile_prob = quantile_prob,
        min_cells = min_cells,
        min_connections = min_connections,
        verbose = verbose
    )

    # Convert to 4D array: (pairs, radii, factors, genes)
    IND_array <- array(NA_real_, dim = c(n_pairs, n_radii, n_factors, n_genes))
    for (r in seq_len(n_radii)) {
        IND_array[, r, , ] <- result_list[[r]]
    }

    # Get valid cell types
    type_counts <- table(data$cell_types)
    valid_types <- names(type_counts[type_counts >= min_cells])

    output <- list(
        IND_array = IND_array,
        pairs = paste0(pairs$sender, "->", pairs$receiver),
        pairs_df = pairs,
        radii = radii,
        factor_genes = factor_genes,
        target_genes = data$gene_names,
        cell_types = valid_types
    )

    # Save to HDF5 if requested
    if (!is.null(output_file)) {
        save_hdf5_sc(output, output_file, verbose = verbose)
    }

    output
}

#' Save Results to HDF5 (Python-Compatible Format)
#'
#' Save batch computation results to HDF5 file format compatible with
#' Python analysis pipelines.
#'
#' @param result Output from compute_batch_sc
#' @param file Path to HDF5 file
#' @param verbose Logical. Print progress (default: TRUE)
#' @return NULL (invisibly)
#' @export
save_hdf5_sc <- function(result, file, verbose = TRUE) {
    if (!requireNamespace("rhdf5", quietly = TRUE)) {
        stop("Package 'rhdf5' required. Install from Bioconductor: ",
             "BiocManager::install('rhdf5')")
    }

    if (file.exists(file)) file.remove(file)
    rhdf5::h5createFile(file)

    # Main data (matching Python name)
    rhdf5::h5write(result$IND_array, file, "morans_i_matrices")

    # Metadata
    rhdf5::h5write(result$pairs, file, "pairs")
    rhdf5::h5write(result$radii, file, "radii")
    rhdf5::h5write(result$factor_genes, file, "factor_genes")
    rhdf5::h5write(result$target_genes, file, "target_genes")
    rhdf5::h5write(result$cell_types, file, "cell_types")

    rhdf5::H5close()

    if (verbose) {
        message("Saved to: ", file)
        message("  Shape: ", paste(dim(result$IND_array), collapse = " x "))
    }

    invisible(NULL)
}

#' Load Results from HDF5
#'
#' Load batch computation results from HDF5 file.
#'
#' @param file Path to HDF5 file
#' @return List with IND_array and metadata
#' @export
load_hdf5_sc <- function(file) {
    if (!requireNamespace("rhdf5", quietly = TRUE)) {
        stop("Package 'rhdf5' required. Install from Bioconductor: ",
             "BiocManager::install('rhdf5')")
    }

    if (!file.exists(file)) {
        stop("File not found: ", file)
    }

    result <- list(
        IND_array = rhdf5::h5read(file, "morans_i_matrices"),
        pairs = rhdf5::h5read(file, "pairs"),
        radii = rhdf5::h5read(file, "radii"),
        factor_genes = rhdf5::h5read(file, "factor_genes"),
        target_genes = rhdf5::h5read(file, "target_genes"),
        cell_types = rhdf5::h5read(file, "cell_types")
    )

    rhdf5::H5close()
    result
}

#' Extract Signature for Specific Pair
#'
#' Extract results for a specific cell type pair from batch results.
#'
#' @param result Output from compute_batch_sc
#' @param sender_type Sender cell type
#' @param receiver_type Receiver cell type
#' @param factor_gene Factor gene name
#' @return Data frame with gene, I values at each radius, etc.
#' @export
extract_pair_signature <- function(result, sender_type, receiver_type, factor_gene) {
    pair_name <- paste0(sender_type, "->", receiver_type)
    pair_idx <- which(result$pairs == pair_name)
    if (length(pair_idx) == 0) {
        stop("Pair '", pair_name, "' not found in results.\n",
             "Available pairs: ", paste(head(result$pairs), collapse = ", "), "...")
    }

    factor_idx <- which(result$factor_genes == factor_gene)
    if (length(factor_idx) == 0) {
        stop("Factor gene '", factor_gene, "' not found in results.\n",
             "Available factors: ", paste(head(result$factor_genes), collapse = ", "), "...")
    }

    # Extract I values: (n_radii, n_genes)
    I_matrix <- result$IND_array[pair_idx, , factor_idx, ]

    df <- data.frame(gene = result$target_genes, stringsAsFactors = FALSE)
    for (r in seq_along(result$radii)) {
        df[[paste0("I_", result$radii[r])]] <- I_matrix[r, ]
    }
    df$I_r1 <- I_matrix[1, ]
    df$I_max <- apply(I_matrix, 2, function(x) {
        if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
    })

    structure(
        df,
        class = c("ScSignature", "data.frame"),
        factor_gene = factor_gene,
        sender_type = sender_type,
        receiver_type = receiver_type,
        radii = result$radii
    )
}
