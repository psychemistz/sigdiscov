# Tests for delta I computation

test_that("compute_delta_i computes positive delta for decreasing signal", {
    # Paracrine pattern: I decreases with distance
    I_curve <- c(0.8, 0.6, 0.4, 0.2, 0.1)
    radii <- c(100, 200, 300, 400, 500)

    result <- compute_delta_i(I_curve, radii)

    expect_true(result$delta_I > 0)
    expect_equal(result$sign, 1L)
    expect_equal(result$I_max, 0.8)
    expect_equal(result$I_min, 0.1)
})

test_that("compute_delta_i computes negative delta for increasing signal", {
    # Inverse pattern: I increases with distance
    I_curve <- c(0.1, 0.2, 0.4, 0.6, 0.8)
    radii <- c(100, 200, 300, 400, 500)

    result <- compute_delta_i(I_curve, radii)

    expect_true(result$delta_I < 0)
    expect_equal(result$sign, -1L)
})

test_that("compute_delta_i handles all-NA input", {
    I_curve <- rep(NA_real_, 5)

    result <- compute_delta_i(I_curve)

    expect_true(is.na(result$delta_I))
    expect_true(is.na(result$sign))
})

test_that("compute_delta_i returns argmax correctly", {
    I_curve <- c(0.1, 0.5, 0.9, 0.3, 0.2)  # Max at position 3
    result <- compute_delta_i(I_curve)

    expect_equal(result$argmax, 3)
})

test_that("compute_delta_i smoothing works", {
    # Need a longer curve for smoothing to take effect
    I_curve <- c(0.8, 0.2, 0.9, 0.1, 0.7, 0.3, 0.6)  # Noisy with max in middle

    result_unsmoothed <- compute_delta_i(I_curve, smooth = FALSE)
    result_smoothed <- compute_delta_i(I_curve, smooth = TRUE, smooth_window = 3)

    # Unsmoothed I_max should be 0.9 (index 3), smoothed should be lower
    expect_equal(result_unsmoothed$I_max, 0.9)
    expect_true(result_smoothed$I_max < 0.9)
})
