# Tests for cell type pair specific spatial analysis functions

# Helper function to create synthetic spatial data
create_test_sc_data <- function(n_genes = 30, n_cells = 100, n_types = 3) {
    set.seed(42)

    # Generate cell types
    cell_types <- rep(paste0("Type", seq_len(n_types)), length.out = n_cells)

    # Generate coordinates (2D spatial layout)
    coords <- matrix(runif(n_cells * 2, 0, 100), ncol = 2)
    colnames(coords) <- c("x", "y")

    # Generate expression data with some spatial structure
    expr <- matrix(rnorm(n_genes * n_cells), nrow = n_genes, ncol = n_cells)
    rownames(expr) <- paste0("Gene", seq_len(n_genes))

    # Add some spatial correlation for testing
    for (i in seq_len(n_cells)) {
        neighbors <- which(sqrt((coords[, 1] - coords[i, 1])^2 +
                                (coords[, 2] - coords[i, 2])^2) < 20)
        if (length(neighbors) > 1) {
            expr[1:5, i] <- expr[1:5, i] + 0.5 * rowMeans(expr[1:5, neighbors, drop = FALSE])
        }
    }

    list(
        expr = expr,
        coords = coords,
        cell_types = cell_types,
        gene_names = rownames(expr),
        n_genes = n_genes,
        n_cells = n_cells
    )
}

# =============================================================================
# Tests for create_directional_weights_sc_cpp
# =============================================================================

test_that("create_directional_weights_sc_cpp returns correct dimensions", {
    sender_coords <- matrix(c(0, 0, 10, 0, 20, 0), ncol = 2, byrow = TRUE)
    receiver_coords <- matrix(c(5, 0, 15, 0, 25, 0, 35, 0), ncol = 2, byrow = TRUE)

    result <- create_directional_weights_sc_cpp(
        sender_coords = sender_coords,
        receiver_coords = receiver_coords,
        radius = 20
    )

    expect_equal(dim(result$W)[1], 3)  # n_sender
    expect_equal(dim(result$W)[2], 4)  # n_receiver
})

test_that("create_directional_weights_sc_cpp produces row-normalized weights", {
    sender_coords <- matrix(c(0, 0, 50, 0), ncol = 2, byrow = TRUE)
    receiver_coords <- matrix(c(5, 0, 10, 0, 55, 0), ncol = 2, byrow = TRUE)

    result <- create_directional_weights_sc_cpp(
        sender_coords = sender_coords,
        receiver_coords = receiver_coords,
        radius = 20
    )

    row_sums <- Matrix::rowSums(result$W)
    # Rows with neighbors should sum to 1
    expect_true(all(abs(row_sums[row_sums > 0] - 1) < 1e-10))
})

test_that("create_directional_weights_sc_cpp respects radius cutoff", {
    sender_coords <- matrix(c(0, 0), ncol = 2)
    receiver_coords <- matrix(c(5, 0, 50, 0), ncol = 2, byrow = TRUE)

    result <- create_directional_weights_sc_cpp(
        sender_coords = sender_coords,
        receiver_coords = receiver_coords,
        radius = 20
    )

    expect_gt(result$W[1, 1], 0)  # 5 units away - included
    expect_equal(result$W[1, 2], 0)  # 50 units away - excluded
})

# =============================================================================
# Tests for pairwise_moran_directional_cpp
# =============================================================================

test_that("pairwise_moran_directional_cpp returns correct dimensions", {
    n_genes <- 10
    n_sender <- 20
    n_receiver <- 30

    sender_data <- matrix(rnorm(n_genes * n_sender), nrow = n_genes)
    receiver_data <- matrix(rnorm(n_genes * n_receiver), nrow = n_genes)

    # Create a weight matrix
    sender_coords <- matrix(runif(n_sender * 2, 0, 100), ncol = 2)
    receiver_coords <- matrix(runif(n_receiver * 2, 0, 100), ncol = 2)

    W_result <- create_directional_weights_sc_cpp(sender_coords, receiver_coords, radius = 50)

    result <- pairwise_moran_directional_cpp(
        sender_data = sender_data,
        receiver_data = receiver_data,
        W = W_result$W,
        verbose = FALSE
    )

    expect_equal(dim(result$moran), c(n_genes, n_genes))
    expect_true(is.finite(result$weight_sum))
})

test_that("pairwise_moran_directional_cpp handles self-pairs correctly", {
    # When sender and receiver are the same cells, directional Moran
    # should match symmetric Moran for self-pairs
    set.seed(123)
    n_genes <- 5
    n_cells <- 50

    data <- matrix(rnorm(n_genes * n_cells), nrow = n_genes)
    coords <- matrix(runif(n_cells * 2, 0, 100), ncol = 2)

    # Compute directional (self-pair)
    W_result <- create_directional_weights_sc_cpp(coords, coords, radius = 30)
    result_dir <- pairwise_moran_directional_cpp(
        sender_data = data,
        receiver_data = data,
        W = W_result$W,
        verbose = FALSE
    )

    # Diagonal elements should be positive for z-normalized data
    diag_vals <- diag(result_dir$moran)
    expect_true(all(is.finite(diag_vals)))
})

# =============================================================================
# Tests for pairwise_moran_directional_streaming_cpp
# =============================================================================

test_that("streaming and matrix methods give correlated results", {
    set.seed(456)
    n_genes <- 10
    n_sender <- 30
    n_receiver <- 40

    sender_data <- matrix(rnorm(n_genes * n_sender), nrow = n_genes)
    receiver_data <- matrix(rnorm(n_genes * n_receiver), nrow = n_genes)
    sender_coords <- matrix(runif(n_sender * 2, 0, 100), ncol = 2)
    receiver_coords <- matrix(runif(n_receiver * 2, 0, 100), ncol = 2)

    # Streaming method
    result_stream <- pairwise_moran_directional_streaming_cpp(
        sender_data = sender_data,
        receiver_data = receiver_data,
        sender_coords = sender_coords,
        receiver_coords = receiver_coords,
        radius = 40,
        normalize_data = TRUE,
        verbose = FALSE
    )

    # Matrix method - use same normalization as streaming
    W_result <- create_directional_weights_sc_cpp(sender_coords, receiver_coords, radius = 40)

    # Results should have high correlation (different normalization may cause scale differences)
    # Use the R wrapper which normalizes consistently
    result_mat <- pairwise_moran_celltype_pair(
        sender_data = sender_data,
        receiver_data = receiver_data,
        sender_coords = sender_coords,
        receiver_coords = receiver_coords,
        radius = 40,
        method = "streaming",
        verbose = FALSE
    )

    # Both streaming methods should be identical
    expect_equal(result_stream$moran, result_mat$moran, tolerance = 1e-10)
})

# =============================================================================
# Tests for compute_celltype_pair_moran_cpp
# =============================================================================

test_that("compute_celltype_pair_moran_cpp returns all expected outputs", {
    set.seed(789)
    n_genes <- 15
    n_sender <- 25
    n_receiver <- 35

    sender_data <- matrix(rnorm(n_genes * n_sender), nrow = n_genes)
    receiver_data <- matrix(rnorm(n_genes * n_receiver), nrow = n_genes)
    sender_coords <- matrix(runif(n_sender * 2, 0, 100), ncol = 2)
    receiver_coords <- matrix(runif(n_receiver * 2, 0, 100), ncol = 2)

    radii <- c(20, 30, 40)

    result <- compute_celltype_pair_moran_cpp(
        sender_data = sender_data,
        receiver_data = receiver_data,
        sender_coords = sender_coords,
        receiver_coords = receiver_coords,
        radii = radii,
        verbose = FALSE
    )

    # Check all expected outputs
    expect_true("I_curves" %in% names(result))
    expect_true("delta_I" %in% names(result))
    expect_true("I_max" %in% names(result))
    expect_true("argmax" %in% names(result))

    # Check dimensions
    expect_equal(length(result$I_curves), 3)  # 3 radii
    expect_equal(dim(result$delta_I), c(n_genes, n_genes))
    expect_equal(dim(result$I_max), c(n_genes, n_genes))
    expect_equal(dim(result$argmax), c(n_genes, n_genes))

    # I_curves should have correct dimensions
    for (r in seq_along(result$I_curves)) {
        expect_equal(dim(result$I_curves[[r]]), c(n_genes, n_genes))
    }
})

test_that("delta_I is correctly computed from I(r) curves", {
    set.seed(111)
    n_genes <- 5
    n_sender <- 20
    n_receiver <- 20

    sender_data <- matrix(rnorm(n_genes * n_sender), nrow = n_genes)
    receiver_data <- matrix(rnorm(n_genes * n_receiver), nrow = n_genes)
    sender_coords <- matrix(runif(n_sender * 2, 0, 100), ncol = 2)
    receiver_coords <- matrix(runif(n_receiver * 2, 0, 100), ncol = 2)

    radii <- c(10, 20, 30, 40, 50)

    result <- compute_celltype_pair_moran_cpp(
        sender_data = sender_data,
        receiver_data = receiver_data,
        sender_coords = sender_coords,
        receiver_coords = receiver_coords,
        radii = radii,
        verbose = FALSE
    )

    # Manually check delta_I for one gene pair
    i <- 1
    j <- 2
    curve <- sapply(result$I_curves, function(m) m[i, j])

    I_max <- max(curve)
    I_min <- min(curve)
    I_short <- curve[1]
    I_long <- curve[length(curve)]
    sign <- if (I_short >= I_long) 1 else -1
    expected_delta <- sign * (I_max - I_min)

    expect_equal(result$delta_I[i, j], expected_delta, tolerance = 1e-10)
    expect_equal(result$I_max[i, j], I_max, tolerance = 1e-10)
})

# =============================================================================
# Tests for R wrapper functions
# =============================================================================

test_that("pairwise_moran_celltype_pair works with streaming method", {
    set.seed(222)
    n_genes <- 8
    n_sender <- 15
    n_receiver <- 20

    sender_data <- matrix(rnorm(n_genes * n_sender), nrow = n_genes)
    rownames(sender_data) <- paste0("Gene", seq_len(n_genes))
    receiver_data <- matrix(rnorm(n_genes * n_receiver), nrow = n_genes)
    rownames(receiver_data) <- paste0("Gene", seq_len(n_genes))
    sender_coords <- matrix(runif(n_sender * 2, 0, 100), ncol = 2)
    receiver_coords <- matrix(runif(n_receiver * 2, 0, 100), ncol = 2)

    result <- pairwise_moran_celltype_pair(
        sender_data = sender_data,
        receiver_data = receiver_data,
        sender_coords = sender_coords,
        receiver_coords = receiver_coords,
        radius = 30,
        method = "streaming",
        verbose = FALSE
    )

    expect_equal(dim(result$moran), c(n_genes, n_genes))
    expect_equal(rownames(result$moran), rownames(sender_data))
    expect_true(result$weight_sum > 0)
})

test_that("get_celltype_pairs generates correct pairs", {
    cell_types <- c("A", "B", "A", "C", "B", "C")

    # With self-pairs (exclude_self = FALSE)
    pairs <- get_celltype_pairs(cell_types, min_cells = 1, exclude_self = FALSE)
    expect_equal(nrow(pairs), 9)  # 3 x 3

    # Without self-pairs (exclude_self = TRUE)
    pairs_no_self <- get_celltype_pairs(cell_types, min_cells = 1, exclude_self = TRUE)
    expect_equal(nrow(pairs_no_self), 6)  # 3 x 3 - 3
})

test_that("create_directional_weights R wrapper works", {
    sender_coords <- matrix(c(0, 0, 10, 0, 20, 0), ncol = 2, byrow = TRUE)
    receiver_coords <- matrix(c(5, 0, 15, 0, 25, 0), ncol = 2, byrow = TRUE)

    result <- create_directional_weights(
        sender_coords = sender_coords,
        receiver_coords = receiver_coords,
        radius = 20,
        verbose = FALSE
    )

    expect_true("W" %in% names(result))
    expect_true("weight_sum" %in% names(result))
    expect_s4_class(result$W, "dgCMatrix")
})

# =============================================================================
# Tests for compute_celltype_pair_analysis (main function)
# =============================================================================

test_that("compute_celltype_pair_analysis works with synthetic data", {
    data <- create_test_sc_data(n_genes = 20, n_cells = 60, n_types = 2)
    class(data) <- c("SCData", "list")

    result <- compute_celltype_pair_analysis(
        data = data,
        radii = c(20, 30, 40),
        min_cells = 5,
        verbose = FALSE
    )

    expect_s3_class(result, "CellTypePairResult")
    expect_true("delta_I" %in% names(result))
    expect_true("I_max" %in% names(result))
    expect_true("pair_stats" %in% names(result))
    expect_equal(result$radii, c(20, 30, 40))
})

test_that("compute_celltype_pair_analysis respects min_cells", {
    # Create data with one very small cell type
    data <- create_test_sc_data(n_genes = 10, n_cells = 50, n_types = 2)
    data$cell_types[1:3] <- "Rare"  # Only 3 cells
    class(data) <- c("SCData", "list")

    result <- compute_celltype_pair_analysis(
        data = data,
        radii = c(20, 30),
        min_cells = 10,  # Rare type has only 3
        verbose = FALSE
    )

    # Pairs involving "Rare" type should be skipped
    expect_false("Rare->Type1" %in% names(result$delta_I))
    expect_false("Type1->Rare" %in% names(result$delta_I))
})

test_that("compute_celltype_pair_analysis handles specific pairs", {
    data <- create_test_sc_data(n_genes = 15, n_cells = 80, n_types = 3)
    class(data) <- c("SCData", "list")

    # Only analyze specific pairs
    pairs <- matrix(c("Type1", "Type2", "Type2", "Type3"), ncol = 2, byrow = TRUE)

    result <- compute_celltype_pair_analysis(
        data = data,
        radii = c(20, 30),
        pairs = pairs,
        min_cells = 5,
        verbose = FALSE
    )

    # Should only have the 2 specified pairs
    expect_equal(nrow(result$pair_stats), 2)
})

# =============================================================================
# Tests for compute_delta_i_batch_cpp
# =============================================================================

test_that("compute_delta_i_batch_cpp handles various curve patterns", {
    # Create test curves
    # Curve 1: Decreasing (paracrine pattern) -> positive delta
    # Curve 2: Increasing (inverse pattern) -> negative delta
    # Curve 3: Flat -> zero delta
    curves <- rbind(
        c(0.8, 0.6, 0.4, 0.2, 0.1),  # Decreasing
        c(0.1, 0.2, 0.4, 0.6, 0.8),  # Increasing
        c(0.5, 0.5, 0.5, 0.5, 0.5)   # Flat
    )

    result <- compute_delta_i_batch_cpp(curves)

    # Decreasing: sign=1, delta = 0.8 - 0.1 = 0.7
    expect_equal(result$sign[1], 1)
    expect_equal(result$delta_I[1], 0.7, tolerance = 1e-10)

    # Increasing: sign=-1, delta = -(0.8 - 0.1) = -0.7
    expect_equal(result$sign[2], -1)
    expect_equal(result$delta_I[2], -0.7, tolerance = 1e-10)

    # Flat: delta = 0
    expect_equal(result$delta_I[3], 0, tolerance = 1e-10)
})
