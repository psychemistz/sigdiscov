#' Create Gaussian Weight Matrix for Single-Cell Data
#'
#' Matches Python implementation: sigma = radius / 3
#'
#' @param sender_coords Sender cell coordinates (n x 2 matrix)
#' @param receiver_coords Receiver cell coordinates (m x 2 matrix)
#' @param radius Outer radius (micrometers)
#' @param inner_radius Inner radius for ring (default: 0)
#' @param sigma Gaussian sigma (default: radius/3)
#' @param min_weight Minimum weight threshold (default: 1e-6)
#' @return Sparse weight matrix (n x m), row-normalized
#' @export
create_weights_sc <- function(sender_coords, receiver_coords,
                               radius, inner_radius = 0, sigma = NULL,
                               min_weight = 1e-6) {
    if (is.null(sigma)) {
        sigma <- radius / 3  # Match Python default
    }

    create_gaussian_weights_cpp(
        as.matrix(sender_coords),
        as.matrix(receiver_coords),
        radius,
        inner_radius,
        sigma,
        min_weight
    )
}

#' Create Gaussian Ring Weight Matrix
#'
#' Weight matrix for annular region between inner and outer radius.
#'
#' @param coords Cell coordinates (n x 2 matrix)
#' @param outer_radius Outer radius (micrometers)
#' @param inner_radius Inner radius (micrometers)
#' @param sigma Gaussian sigma (default: outer_radius/3)
#' @return Sparse weight matrix (n x n), row-normalized
#' @export
create_ring_weights_sc <- function(coords, outer_radius, inner_radius,
                                    sigma = NULL) {
    if (is.null(sigma)) {
        sigma <- outer_radius / 3
    }

    create_gaussian_ring_weights_cpp(
        as.matrix(coords),
        outer_radius,
        inner_radius,
        sigma
    )
}


#' Create Scalable Weight Matrix for Large Single-Cell Datasets
#'
#' Uses native C++ KD-tree (nanoflann) for efficient neighbor search to create
#' sparse weight matrices for datasets with >10,000 cells. Falls back to
#' brute-force for smaller datasets.
#'
#' @param coords Cell coordinates (n x 2 matrix)
#' @param radius Radius for neighbor search
#' @param sigma Gaussian sigma (default: radius/3)
#' @param max_neighbors Maximum neighbors per cell (default: 200)
#' @param threshold Cell count threshold for using KD-tree (default: 10000)
#' @param verbose Print progress messages (default: TRUE)
#' @return List with W (sparse weight matrix), weight_sum, and n_edges
#'
#' @details
#' For datasets with more than \code{threshold} cells, uses nanoflann C++
#' KD-tree for O(n log n) neighbor search instead of O(n²) brute-force.
#' This enables analysis of datasets with hundreds of thousands of cells.
#' No external R package dependencies required.
#'
#' @examples
#' \dontrun{
#' # Large CosMx dataset (~100k cells)
#' coords <- as.matrix(meta[, c("x", "y")])
#' W_result <- create_weights_sc_large(coords, radius = 50, sigma = 20)
#' }
#'
#' @export
create_weights_sc_large <- function(coords, radius, sigma = NULL,
                                     max_neighbors = 200, threshold = 10000,
                                     verbose = TRUE) {
    coords <- as.matrix(coords)
    n <- nrow(coords)

    if (is.null(sigma)) {
        sigma <- radius / 3
    }

    # For small datasets, use direct method
    if (n <= threshold) {
        if (verbose) message("Using direct method (n <= ", threshold, ")")
        W <- create_weights_sc(coords, coords, radius, sigma = sigma)
        return(list(W = W, weight_sum = sum(W), n_edges = length(W@x)))
    }

    # For large datasets, use native C++ KD-tree
    if (verbose) message("Using native C++ KD-tree for neighbor search (n = ", n, ")")

    # Find neighbors and build weight matrix using nanoflann
    if (verbose) message("Finding neighbors within radius ", radius, "...")
    t1 <- Sys.time()
    nn <- find_neighbors_radius_cpp(coords, radius, max_neighbors)
    t2 <- Sys.time()
    if (verbose) message("Neighbor search: ", round(difftime(t2, t1, units = "secs"), 2), " sec")

    # Create weight matrix from neighbors
    if (verbose) message("Building sparse weight matrix...")
    result <- create_weights_from_neighbors_cpp(
        nn$nn.idx,
        nn$nn.dists,
        sigma,
        radius,
        row_normalize = TRUE,
        self_weight = 0.0
    )

    if (verbose) {
        message("W non-zeros: ", result$n_edges)
        message("Avg neighbors/cell: ", round(result$n_edges / n, 1))
    }

    result
}


#' Compute Pairwise Moran's I for Large Single-Cell Datasets
#'
#' Scalable computation of pairwise Moran's I using native C++ KD-tree
#' (nanoflann) for neighbor search. Handles datasets with hundreds of
#' thousands of cells. No external R package dependencies required.
#'
#' @param data Gene expression matrix (genes x cells)
#' @param coords Cell coordinates (n x 2 matrix)
#' @param radius Radius for neighbor search
#' @param sigma Gaussian sigma (default: radius/3)
#' @param method Method to use: "streaming" (default, no neighbor limit) or
#'        "matrix" (builds explicit W matrix with neighbor limit)
#' @param max_neighbors Maximum neighbors per cell when method="matrix" (default: 200)
#' @param verbose Print progress messages (default: TRUE)
#' @return List with moran matrix, weight_sum, and n_edges
#'
#' @details
#' Two methods are available:
#' \itemize{
#'   \item \code{"streaming"} (default): Computes spatial lag on-the-fly without
#'     storing the weight matrix. No neighbor limit, memory efficient, fastest.
#'   \item \code{"matrix"}: Builds explicit sparse weight matrix. Uses max_neighbors
#'     limit to control memory usage.
#' }
#'
#' Both methods use nanoflann (header-only C++ KD-tree) for O(n log n) neighbor
#' search. No external R package dependencies required.
#'
#' @examples
#' \dontrun{
#' # CosMx dataset (~100k cells) - streaming (recommended)
#' result <- pairwise_moran_sc_large(vst_data, coords, radius = 0.1)
#'
#' # With explicit W matrix and neighbor limit
#' result <- pairwise_moran_sc_large(vst_data, coords, radius = 0.1,
#'                                    method = "matrix", max_neighbors = 200)
#' }
#'
#' @export
pairwise_moran_sc_large <- function(data, coords, radius, sigma = NULL,
                                     method = c("auto", "dense", "streaming", "matrix"),
                                     max_neighbors = 200, chunk_size = 500,
                                     verbose = TRUE) {
    data <- as.matrix(data)
    coords <- as.matrix(coords)
    n <- nrow(coords)
    method <- match.arg(method)

    if (ncol(data) != n) {
        stop("Number of columns in data (", ncol(data),
             ") must match number of rows in coords (", n, ")")
    }

    if (is.null(sigma)) {
        sigma <- radius / 3
    }

    # Auto-select method based on expected density
    if (method == "auto") {
        # Estimate neighbors by sampling
        coord_range <- max(diff(range(coords[,1])), diff(range(coords[,2])))
        density_estimate <- (pi * radius^2) / (coord_range^2) * n

        if (verbose) {
            message(sprintf("Estimated neighbors/cell: %.0f", density_estimate))
        }

        # Use dense BLAS when > 500 neighbors (about 0.5% density for 98k cells)
        method <- if (density_estimate > 500) "dense" else "streaming"
        if (verbose) message(sprintf("Auto-selected method: %s", method))
    }

    # Choose implementation based on method
    if (method == "dense") {
        # Dense BLAS: chunked matrix operations, best for large radii
        result <- pairwise_moran_dense_cpp(
            data,
            coords,
            radius,
            sigma,
            chunk_size,
            verbose
        )
    } else if (method == "streaming") {
        # Streaming: no W matrix storage, good for small radii
        result <- pairwise_moran_streaming_cpp(
            data,
            coords,
            radius,
            sigma,
            verbose
        )
    } else {
        # Matrix: builds explicit W with neighbor limit
        result <- pairwise_moran_native_cpp(
            data,
            coords,
            radius,
            sigma,
            max_neighbors,
            verbose
        )
    }

    # Add gene names
    if (!is.null(rownames(data))) {
        rownames(result$moran) <- rownames(data)
        colnames(result$moran) <- rownames(data)
    }

    result
}
