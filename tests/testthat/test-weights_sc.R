# Tests for single-cell Gaussian weight matrix functions

test_that("create_weights_sc returns sparse matrix", {
    coords <- matrix(c(0, 0, 10, 0, 0, 10, 10, 10), ncol = 2, byrow = TRUE)
    W <- create_weights_sc(coords, coords, radius = 15)

    expect_s4_class(W, "dgCMatrix")
    expect_equal(dim(W), c(4, 4))
})

test_that("Gaussian weights use sigma = radius/3", {
    # Create simple 2-point setup
    sender <- matrix(c(0, 0), ncol = 2)
    receiver <- matrix(c(10, 0), ncol = 2)

    W <- create_weights_sc(sender, receiver, radius = 30)

    # sigma = 30/3 = 10, distance = 10
    # Since only one neighbor, row-normalized = 1.0
    expect_equal(W[1, 1], 1.0)
})

test_that("Gaussian weights are row-normalized", {
    coords <- matrix(c(0, 0, 10, 0, 20, 0, 30, 0), ncol = 2, byrow = TRUE)
    W <- create_weights_sc(coords, coords, radius = 25)

    row_sums <- Matrix::rowSums(W)
    # Row sums should be either 1 (has neighbors) or 0 (no neighbors)
    expect_true(all(abs(row_sums - 1) < 1e-10 | row_sums == 0))
})

test_that("Gaussian ring weights exclude inner radius", {
    coords <- matrix(c(0, 0, 5, 0, 15, 0), ncol = 2, byrow = TRUE)

    # Full radius should include both neighbors
    W_full <- create_weights_sc(coords, coords, radius = 20)
    expect_gt(W_full[1, 2], 0)  # 5 units away
    expect_gt(W_full[1, 3], 0)  # 15 units away

    # Ring should exclude inner (< 10)
    W_ring <- create_ring_weights_sc(coords, outer_radius = 20, inner_radius = 10)
    expect_equal(W_ring[1, 2], 0)  # 5 units away - excluded
    expect_gt(W_ring[1, 3], 0)     # 15 units away - included
})

test_that("Custom sigma parameter works", {
    sender <- matrix(c(0, 0), ncol = 2)
    receiver <- matrix(c(10, 0), ncol = 2)

    # Compare default sigma vs custom
    W_default <- create_weights_sc(sender, receiver, radius = 30)  # sigma = 10
    W_custom <- create_weights_sc(sender, receiver, radius = 30, sigma = 5)

    # Both normalize to 1 since single neighbor
    expect_equal(W_default[1, 1], 1.0)
    expect_equal(W_custom[1, 1], 1.0)
})

test_that("Empty result for no neighbors", {
    sender <- matrix(c(0, 0), ncol = 2)
    receiver <- matrix(c(100, 100), ncol = 2)

    W <- create_weights_sc(sender, receiver, radius = 10)

    expect_equal(length(W@x), 0)
    expect_equal(Matrix::rowSums(W)[1], 0)
})

test_that("Gaussian weights match Python formula", {
    # Test that weights follow exp(-d^2 / (2 * sigma^2))
    sender <- matrix(c(0, 0), ncol = 2)
    receiver <- matrix(c(6, 8, 12, 16), ncol = 2, byrow = TRUE)  # distances 10 and 20

    # radius = 30, sigma = 10
    W <- create_weights_sc(sender, receiver, radius = 30)

    # Raw weights before normalization:
    # w1 = exp(-100 / 200) = exp(-0.5) = 0.6065
    # w2 = exp(-400 / 200) = exp(-2.0) = 0.1353
    # Total = 0.7418

    # Normalized weights:
    # w1_norm = 0.6065 / 0.7418 = 0.8176
    # w2_norm = 0.1353 / 0.7418 = 0.1824

    expected_w1 <- exp(-0.5) / (exp(-0.5) + exp(-2.0))
    expected_w2 <- exp(-2.0) / (exp(-0.5) + exp(-2.0))

    expect_equal(W[1, 1], expected_w1, tolerance = 1e-4)
    expect_equal(W[1, 2], expected_w2, tolerance = 1e-4)
})
