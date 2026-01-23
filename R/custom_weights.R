#' @title Custom Weight Matrix Support
#' @description Functions for loading and using custom spatial weight matrices.
#' @name custom_weights
NULL

#' Load Custom Weight Matrix from File
#'
#' Loads a spatial weight matrix from a file. Supports multiple formats
#' including dense TSV, sparse COO, and sparse TSV.
#'
#' @param file Path to weight matrix file.
#' @param format Format of the file: "auto" (default), "dense", "sparse_coo",
#'   or "sparse_tsv". If "auto", the format is detected from file content.
#' @param spot_names Optional vector of spot names to use. If NULL, names are
#'   read from the file.
#'
#' @return A square matrix (dense or sparse Matrix) with spot names as
#'   row and column names.
#'
#' @details
#' Supported formats:
#' \describe{
#'   \item{dense}{Full matrix with spot names as header and first column.
#'     Tab-separated values.}
#'   \item{sparse_coo}{Coordinate format with three columns: row_idx, col_idx,
#'     weight. Indices are 1-based or spot names.}
#'   \item{sparse_tsv}{Three columns: spot1, spot2, weight. Spot names as
#'     strings.}
#' }
#'
#' The matrix is symmetrized if not already symmetric.
#'
#' @examples
#' \dontrun{
#' # Load dense matrix
#' W <- load_weight_matrix("weights_dense.tsv", format = "dense")
#'
#' # Load sparse matrix
#' W <- load_weight_matrix("weights_sparse.tsv", format = "sparse_tsv")
#'
#' # Use with pairwise_moran_custom
#' result <- pairwise_moran_custom(data, W)
#' }
#'
#' @export
load_weight_matrix <- function(file, format = "auto", spot_names = NULL) {
    stopifnot(file.exists(file))

    # Read first few lines to detect format
    sample_lines <- readLines(file, n = 5)
    first_line <- sample_lines[1]

    # Detect delimiter (tab or comma)
    if (grepl("\t", first_line)) {
        sep <- "\t"
    } else if (grepl(",", first_line)) {
        sep <- ","
    } else {
        sep <- "\t"  # default
    }

    fields <- strsplit(first_line, sep)[[1]]
    n_fields <- length(fields)

    # Auto-detect format
    if (format == "auto") {
        if (n_fields == 3) {
            # Could be sparse_coo or sparse_tsv
            # Check if first field looks like a coordinate
            if (grepl("^[0-9]+$", fields[1]) ||
                grepl("x", fields[1], ignore.case = TRUE)) {
                # Check second line to distinguish
                if (length(sample_lines) > 1) {
                    second_fields <- strsplit(sample_lines[2], sep)[[1]]
                    if (all(grepl("^[0-9]+$", second_fields[1:2]))) {
                        format <- "sparse_coo"
                    } else {
                        format <- "sparse_tsv"
                    }
                } else {
                    format <- "sparse_tsv"
                }
            } else {
                format <- "sparse_tsv"
            }
        } else if (n_fields > 10) {
            format <- "dense"
        } else {
            # Assume sparse TSV for safety
            format <- "sparse_tsv"
        }
        message(sprintf("Auto-detected format: %s", format))
    }

    # Load based on format
    if (format == "dense") {
        W <- load_dense_weight_matrix_internal(file, sep)
    } else if (format == "sparse_coo") {
        W <- load_sparse_coo_matrix_internal(file, sep, spot_names)
    } else if (format == "sparse_tsv") {
        W <- load_sparse_tsv_matrix_internal(file, sep)
    } else {
        stop("Unknown format: ", format)
    }

    # Ensure symmetry
    if (!isSymmetric(as.matrix(W), tol = 1e-10)) {
        message("Symmetrizing weight matrix...")
        W <- (W + t(W)) / 2
    }

    W
}


#' Load Dense Weight Matrix
#' @keywords internal
load_dense_weight_matrix_internal <- function(file, sep) {
    df <- read.table(file, header = TRUE, row.names = 1, sep = sep,
                     stringsAsFactors = FALSE, check.names = FALSE)
    W <- as.matrix(df)

    # Ensure numeric
    storage.mode(W) <- "double"

    W
}


#' Load Sparse COO Matrix
#' @keywords internal
load_sparse_coo_matrix_internal <- function(file, sep, spot_names = NULL) {
    df <- read.table(file, header = TRUE, sep = sep, stringsAsFactors = FALSE)

    if (ncol(df) < 3) {
        stop("Sparse COO format requires 3 columns: row, col, weight")
    }

    row_idx <- df[[1]]
    col_idx <- df[[2]]
    weights <- as.numeric(df[[3]])

    # Determine if indices are numeric or names
    if (all(grepl("^[0-9]+$", row_idx))) {
        # Numeric indices
        row_idx <- as.integer(row_idx)
        col_idx <- as.integer(col_idx)
        n <- max(c(row_idx, col_idx))

        if (is.null(spot_names)) {
            spot_names <- paste0("spot_", seq_len(n))
        }
    } else {
        # Spot names
        all_spots <- unique(c(row_idx, col_idx))
        if (is.null(spot_names)) {
            spot_names <- sort(all_spots)
        }
        n <- length(spot_names)

        spot_idx <- setNames(seq_along(spot_names), spot_names)
        row_idx <- spot_idx[row_idx]
        col_idx <- spot_idx[col_idx]
    }

    # Create sparse matrix
    W <- Matrix::sparseMatrix(
        i = row_idx,
        j = col_idx,
        x = weights,
        dims = c(n, n),
        dimnames = list(spot_names, spot_names)
    )

    W
}


#' Load Sparse TSV Matrix
#' @keywords internal
load_sparse_tsv_matrix_internal <- function(file, sep) {
    df <- read.table(file, header = TRUE, sep = sep, stringsAsFactors = FALSE)

    if (ncol(df) < 3) {
        stop("Sparse TSV format requires 3 columns: spot1, spot2, weight")
    }

    spot1 <- as.character(df[[1]])
    spot2 <- as.character(df[[2]])
    weights <- as.numeric(df[[3]])

    # Get unique spots
    all_spots <- unique(c(spot1, spot2))
    n <- length(all_spots)
    spot_idx <- setNames(seq_along(all_spots), all_spots)

    # Create sparse matrix
    W <- Matrix::sparseMatrix(
        i = spot_idx[spot1],
        j = spot_idx[spot2],
        x = weights,
        dims = c(n, n),
        dimnames = list(all_spots, all_spots)
    )

    W
}


#' Compute Moran's I with Custom Weight Matrix
#'
#' Computes pairwise Moran's I using a user-provided spatial weight matrix
#' instead of computing weights from coordinates.
#'
#' @param data Gene expression matrix (genes x spots). Column names should
#'   match the row/column names of the weight matrix.
#' @param weight_matrix Square weight matrix. Can be dense matrix or sparse
#'   Matrix object.
#' @param normalize Logical. If TRUE, divide each weight by S0 before
#'   computation. Default: FALSE.
#' @param mode Calculation mode: "paired" (all pairs), "single" (univariate
#'   only). Default: "paired".
#' @param verbose Print progress messages. Default: TRUE.
#'
#' @return A list containing:
#'   \item{moran}{Moran's I matrix (if mode = "paired") or named vector
#'     (if mode = "single")}
#'   \item{gene_names}{Names of genes}
#'   \item{S0}{Sum of weights used}
#'
#' @details
#' This function allows using arbitrary spatial weight matrices, enabling:
#' \itemize{
#'   \item Non-standard spatial relationships
#'   \item Alternative distance metrics
#'   \item Prior biological knowledge about connectivity
#' }
#'
#' The weight matrix should typically be symmetric and have zeros on the
#' diagonal (no self-connections).
#'
#' @examples
#' \dontrun{
#' # Load custom weights
#' W <- load_weight_matrix("custom_weights.tsv")
#'
#' # Compute Moran's I
#' result <- pairwise_moran_custom(data, W)
#'
#' # Get significant pairs
#' high_corr <- which(result$moran > 0.5 & lower.tri(result$moran),
#'                    arr.ind = TRUE)
#' }
#'
#' @export
pairwise_moran_custom <- function(
    data,
    weight_matrix,
    normalize = FALSE,
    mode = "paired",
    verbose = TRUE
) {
    # Convert data to matrix if needed
    if (!is.matrix(data)) {
        data <- as.matrix(data)
    }

    # Check if weight matrix is sparse
    is_sparse <- methods::is(weight_matrix, "sparseMatrix")

    # Validate dimensions
    n_spots <- ncol(data)

    if (nrow(weight_matrix) != ncol(weight_matrix)) {
        stop("Weight matrix must be square")
    }
    if (n_spots != nrow(weight_matrix)) {
        stop(sprintf(
            "Number of spots in data (%d) must match weight matrix dimension (%d)",
            n_spots, nrow(weight_matrix)
        ))
    }

    # Check for matching spot names if available
    data_spots <- colnames(data)
    weight_spots <- rownames(weight_matrix)

    if (!is.null(data_spots) && !is.null(weight_spots)) {
        if (!all(data_spots %in% weight_spots)) {
            warning("Some spot names in data not found in weight matrix")
        }
        # Reorder weight matrix to match data
        common_spots <- intersect(data_spots, weight_spots)
        if (length(common_spots) < n_spots) {
            stop(sprintf(
                "Only %d spots match between data and weight matrix",
                length(common_spots)
            ))
        }
        weight_matrix <- weight_matrix[data_spots, data_spots]
    }

    # Normalize if requested
    if (normalize) {
        S0 <- sum(weight_matrix)
        weight_matrix <- weight_matrix / S0
    }

    # Convert mode to integer for C++
    mode_int <- switch(mode,
        "paired" = 0L,
        "single" = 1L,
        "first" = 2L,
        0L  # default to paired
    )

    # Call C++ function (handles z-normalization with population SD internally)
    if (is_sparse) {
        cpp_result <- cpp_pairwise_moran_custom_sparse(data, weight_matrix, mode_int, verbose)
    } else {
        W <- as.matrix(weight_matrix)
        cpp_result <- cpp_pairwise_moran_custom(data, W, mode_int, verbose)
    }

    # Format result
    gene_names <- rownames(data)
    moran <- cpp_result$moran

    if (mode == "paired" && !is.null(gene_names)) {
        rownames(moran) <- gene_names
        colnames(moran) <- gene_names
    } else if (!is.null(gene_names)) {
        names(moran) <- gene_names
    }

    list(
        moran = moran,
        gene_names = gene_names,
        weight_sum = cpp_result$weight_sum
    )
}


#' Save Weight Matrix to File
#'
#' Saves a weight matrix to a file in various formats.
#'
#' @param weight_matrix Weight matrix to save.
#' @param file Output file path.
#' @param format Output format: "dense", "sparse_tsv", or "sparse_coo".
#'   Default: "sparse_tsv".
#' @param threshold Minimum weight to include in sparse formats. Default: 1e-10.
#'
#' @export
save_weight_matrix <- function(
    weight_matrix,
    file,
    format = "sparse_tsv",
    threshold = 1e-10
) {
    W <- as.matrix(weight_matrix)
    spot_names <- rownames(W)
    if (is.null(spot_names)) {
        spot_names <- paste0("spot_", seq_len(nrow(W)))
    }

    if (format == "dense") {
        df <- as.data.frame(W)
        rownames(df) <- spot_names
        colnames(df) <- spot_names
        write.table(df, file, sep = "\t", quote = FALSE, row.names = TRUE)

    } else if (format == "sparse_tsv") {
        # Find non-zero entries (upper triangle only for symmetric)
        idx <- which(W > threshold & upper.tri(W, diag = FALSE), arr.ind = TRUE)

        df <- data.frame(
            spot1 = spot_names[idx[, 1]],
            spot2 = spot_names[idx[, 2]],
            weight = W[idx],
            stringsAsFactors = FALSE
        )

        write.table(df, file, sep = "\t", quote = FALSE, row.names = FALSE)

    } else if (format == "sparse_coo") {
        idx <- which(W > threshold & upper.tri(W, diag = FALSE), arr.ind = TRUE)

        df <- data.frame(
            row = idx[, 1],
            col = idx[, 2],
            weight = W[idx]
        )

        write.table(df, file, sep = "\t", quote = FALSE, row.names = FALSE)

    } else {
        stop("Unknown format: ", format)
    }

    invisible(file)
}


#' Create Weight Matrix from Distance Matrix
#'
#' Creates a spatial weight matrix from a pre-computed distance matrix
#' using RBF kernel.
#'
#' @param dist_matrix Distance matrix (spots x spots).
#' @param sigma RBF kernel bandwidth. Default: 100.
#' @param max_dist Maximum distance to include. Weights for larger distances
#'   are set to zero. Default: Inf (no cutoff).
#' @param self_weight Weight for self-connections (diagonal). Default: 0.
#'
#' @return Weight matrix with same dimensions as input.
#'
#' @export
weight_matrix_from_distance <- function(
    dist_matrix,
    sigma = 100,
    max_dist = Inf,
    self_weight = 0
) {
    D <- as.matrix(dist_matrix)

    # Apply RBF kernel
    W <- exp(-(D^2) / (2 * sigma^2))

    # Apply distance cutoff
    W[D > max_dist] <- 0

    # Set diagonal
    diag(W) <- self_weight

    # Copy names
    if (!is.null(rownames(D))) {
        rownames(W) <- rownames(D)
        colnames(W) <- colnames(D)
    }

    W
}
