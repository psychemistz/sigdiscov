# Tests for metric computation functions

test_that("compute_moran_from_lag computes correctly", {
    # Perfect positive autocorrelation
    z_f <- c(1, 1, -1, -1)
    lag_g <- c(1, 1, -1, -1)
    I <- compute_moran_from_lag(z_f, lag_g)
    expect_equal(I, 1, tolerance = 1e-10)

    # Perfect negative autocorrelation
    z_f <- c(1, -1, 1, -1)
    lag_g <- c(-1, 1, -1, 1)
    I <- compute_moran_from_lag(z_f, lag_g)
    expect_equal(I, -1, tolerance = 1e-10)

    # No correlation (orthogonal)
    z_f <- c(1, 0, -1, 0)
    lag_g <- c(0, 1, 0, -1)
    I <- compute_moran_from_lag(z_f, lag_g)
    expect_equal(I, 0, tolerance = 1e-10)
})

test_that("compute_ind_from_lag is bounded [-1, 1]", {
    set.seed(42)
    for (i in 1:100) {
        z_f <- rnorm(50)
        lag_g <- rnorm(50)
        I <- compute_ind_from_lag(z_f, lag_g)
        expect_true(is.na(I) || (I >= -1 && I <= 1))
    }
})

test_that("compute_ind_from_lag returns NA for zero vectors", {
    z_f <- rep(0, 10)
    lag_g <- rnorm(10)
    expect_true(is.na(compute_ind_from_lag(z_f, lag_g)))

    z_f <- rnorm(10)
    lag_g <- rep(0, 10)
    expect_true(is.na(compute_ind_from_lag(z_f, lag_g)))
})

test_that("compute_metric_batch returns correct length", {
    n <- 100
    n_genes <- 50
    z_f <- rnorm(n)
    lag_G <- matrix(rnorm(n * n_genes), n, n_genes)

    result_moran <- compute_metric_batch(z_f, lag_G, "moran")
    expect_equal(length(result_moran), n_genes)

    result_ind <- compute_metric_batch(z_f, lag_G, "ind")
    expect_equal(length(result_ind), n_genes)
})

test_that("compute_spatial_lag produces correct output", {
    library(Matrix)
    # Simple 3x3 weight matrix
    W <- sparseMatrix(
        i = c(1, 1, 2, 2, 3, 3),
        j = c(2, 3, 1, 3, 1, 2),
        x = c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
        dims = c(3, 3)
    )
    z <- c(1, 2, 3)

    lag <- compute_spatial_lag(W, z)

    # lag[1] = 0.5*2 + 0.5*3 = 2.5
    # lag[2] = 0.5*1 + 0.5*3 = 2.0
    # lag[3] = 0.5*1 + 0.5*2 = 1.5
    expect_equal(lag, c(2.5, 2.0, 1.5))
})
