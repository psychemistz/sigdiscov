# Tests for utility functions

test_that("standardize produces mean 0 and population sd 1", {
    x <- c(1, 2, 3, 4, 5)
    z <- standardize(x)

    # Population SD: sqrt(sum((x - mean)^2) / n)
    pop_sd <- sqrt(sum((z - mean(z))^2) / length(z))

    expect_equal(mean(z), 0, tolerance = 1e-10)
    expect_equal(pop_sd, 1, tolerance = 1e-10)
})

test_that("standardize handles constant vectors", {
    x <- rep(5, 10)
    z <- standardize(x)

    # Should be all zeros (or very close due to numerical stability)
    expect_true(all(abs(z) < 1))
})

test_that("standardize_matrix standardizes row-wise", {
    X <- matrix(1:20, nrow = 4, ncol = 5)
    X_norm <- standardize_matrix(X)

    # Each row should have mean ~0 and population sd ~1
    row_means <- rowMeans(X_norm)
    row_pop_sds <- apply(X_norm, 1, function(r) sqrt(sum((r - mean(r))^2) / length(r)))

    expect_true(all(abs(row_means) < 1e-10))
    expect_true(all(abs(row_pop_sds - 1) < 1e-10))
})

test_that("adjust_pvalues applies correction", {
    p <- c(0.01, 0.02, 0.03, 0.04, 0.05)

    p_adj <- adjust_pvalues(p, method = "BH")

    # BH adjusted p-values should be >= raw p-values
    expect_true(all(p_adj >= p))
})

test_that("sparse_row_normalize normalizes rows", {
    library(Matrix)
    W <- sparseMatrix(
        i = c(1, 1, 2, 2, 3),
        j = c(2, 3, 1, 3, 1),
        x = c(2, 4, 1, 2, 3),
        dims = c(3, 3)
    )

    W_norm <- sparse_row_normalize(W)

    row_sums <- rowSums(W_norm)

    # Rows with entries should sum to 1
    expect_equal(row_sums[1], 1, tolerance = 1e-10)
    expect_equal(row_sums[2], 1, tolerance = 1e-10)
    expect_equal(row_sums[3], 1, tolerance = 1e-10)
})
