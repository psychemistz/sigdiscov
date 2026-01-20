# Tests for Visium weight matrix functions

test_that("create_weights_visium produces row-normalized matrix", {
    coords <- matrix(c(0, 0, 100, 0, 0, 100, 100, 100), ncol = 2, byrow = TRUE)
    W <- create_weights_visium(coords, radius = 150)

    # Should be sparse matrix
    expect_true(inherits(W, "sparseMatrix"))

    # Each row should sum to 1 (or 0 if no neighbors)
    row_sums <- Matrix::rowSums(W)
    expect_true(all(abs(row_sums - 1) < 1e-10 | row_sums == 0))
})

test_that("create_weights_visium respects radius", {
    coords <- matrix(c(0, 0, 100, 0, 200, 0), ncol = 2, byrow = TRUE)

    # Small radius: only immediate neighbors
    W_small <- create_weights_visium(coords, radius = 110)
    expect_equal(sum(W_small[1, ] > 0), 1)  # Spot 1 has 1 neighbor

    # Large radius: all neighbors
    W_large <- create_weights_visium(coords, radius = 250)
    expect_equal(sum(W_large[1, ] > 0), 2)  # Spot 1 has 2 neighbors
})

test_that("create_weights_visium excludes self by default", {
    coords <- matrix(c(0, 0, 100, 0), ncol = 2, byrow = TRUE)
    W <- create_weights_visium(coords, radius = 200)

    # Diagonal should be 0
    expect_equal(diag(as.matrix(W)), c(0, 0))
})

test_that("create_weights_visium includes self when requested", {
    coords <- matrix(c(0, 0, 100, 0), ncol = 2, byrow = TRUE)
    W <- create_weights_visium(coords, radius = 200, include_self = TRUE)

    # Diagonal should be non-zero
    expect_true(all(diag(as.matrix(W)) > 0))
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
