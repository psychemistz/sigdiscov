# Tests for Visium weight matrix functions
# Note: create_weights_visium returns a list with W (weight matrix) and weight_sum

test_that("create_weights_visium produces valid weight matrix", {
    # Use grid coordinates (row, col integers) for gaussian weights
    coords <- data.frame(row = c(1, 2, 1, 2), col = c(1, 1, 2, 2))
    result <- create_weights_visium(coords, radius = 150)

    # Should return a list
    expect_true(is.list(result))
    expect_true("W" %in% names(result))
    expect_true("weight_sum" %in% names(result))

    # W should be sparse matrix
    expect_true(inherits(result$W, "sparseMatrix"))
})

test_that("create_weights_visium respects radius", {
    # Use grid coordinates (row, col integers) for gaussian weights
    coords <- data.frame(row = c(1, 1, 1), col = c(1, 2, 3))

    # Small radius: only immediate neighbors (radius 100 = 1 grid unit)
    result_small <- create_weights_visium(coords, radius = 100)
    expect_true(sum(result_small$W[1, ] > 0) >= 1)  # Spot 1 has at least 1 neighbor

    # Large radius: all neighbors (radius 300 = 3 grid units)
    result_large <- create_weights_visium(coords, radius = 300)
    expect_true(sum(result_large$W[1, ] > 0) >= 2)  # Spot 1 has 2 neighbors
})

test_that("create_weights_visium excludes self by default", {
    # Use grid coordinates (row, col integers)
    coords <- data.frame(row = c(1, 1), col = c(1, 2))
    result <- create_weights_visium(coords, radius = 200)

    # Diagonal should be 0
    expect_equal(diag(as.matrix(result$W)), c(0, 0))
})

test_that("create_weights_visium includes self when requested", {
    # Use grid coordinates (row, col integers)
    coords <- data.frame(row = c(1, 1), col = c(1, 2))
    result <- create_weights_visium(coords, radius = 200, include_self = TRUE)

    # Diagonal should be non-zero
    expect_true(all(diag(as.matrix(result$W)) > 0))
})

test_that("create_ring_weights_visium excludes inner region", {
    coords <- matrix(c(0, 0, 50, 0, 150, 0), ncol = 2, byrow = TRUE)

    W <- create_ring_weights_visium(coords, inner_radius = 100, outer_radius = 200)

    # Spot 1 should only connect to spot 3 (distance 150), not spot 2 (distance 50)
    expect_equal(W[1, 2], 0)
    expect_true(W[1, 3] > 0)
})

test_that("create_directional_weights_visium has correct dimensions", {
    sender_coords <- matrix(c(0, 0, 100, 0), ncol = 2, byrow = TRUE)
    receiver_coords <- matrix(c(0, 100, 100, 100, 200, 100), ncol = 2, byrow = TRUE)

    W <- create_directional_weights_visium(sender_coords, receiver_coords, radius = 150)

    # Should be n_senders x n_receivers
    expect_equal(nrow(W), 2)
    expect_equal(ncol(W), 3)
})

test_that("create_directional_weights_visium is row-normalized", {
    sender_coords <- matrix(runif(20), 10, 2) * 100
    receiver_coords <- matrix(runif(40), 20, 2) * 100

    W <- create_directional_weights_visium(sender_coords, receiver_coords, radius = 50)

    row_sums <- Matrix::rowSums(W)
    expect_true(all(abs(row_sums - 1) < 1e-10 | row_sums == 0))
})
