#' Compute Delta I from I(r) Curve
#'
#' Computes the signed delta I statistic from a vector of I values at
#' different distance radii. Delta I captures the distance-dependent
#' change in spatial correlation.
#'
#' @param I_curve Numeric vector. I values at different radii.
#' @param radii Numeric vector. Distance radii (optional, for reference).
#' @param smooth Logical. Apply moving average smoothing (default: FALSE).
#' @param smooth_window Integer. Window size for smoothing (default: 5).
#'
#' @return A list with components:
#' \describe{
#'   \item{delta_I}{Signed delta I value}
#'   \item{sign}{Direction: +1 if short-range > long-range, -1 otherwise}
#'   \item{I_max}{Maximum I value across radii}
#'   \item{I_min}{Minimum I value across radii}
#'   \item{I_short}{I value at first (shortest) radius}
#'   \item{I_long}{I value at last (longest) radius}
#'   \item{argmax}{Index of radius with maximum I}
#' }
#'
#' @details
#' The signed delta I is computed as:
#' \deqn{\Delta I = sign(I_{r1} - I_{rn}) \times (I_{max} - I_{min})}
#'
#' Interpretation:
#' \itemize{
#'   \item Positive delta I: Signal decreases with distance (paracrine pattern)
#'   \item Negative delta I: Signal increases with distance (inverse pattern)
#' }
#'
#' @examples
#' # Paracrine-like pattern (decreasing with distance)
#' I_curve <- c(0.8, 0.6, 0.4, 0.2, 0.1)
#' radii <- c(100, 200, 300, 400, 500)
#' result <- compute_delta_i(I_curve, radii)
#' result$delta_I  # Positive
#'
#' # Inverse pattern (increasing with distance)
#' I_curve <- c(0.1, 0.2, 0.4, 0.6, 0.8)
#' result <- compute_delta_i(I_curve, radii)
#' result$delta_I  # Negative
#'
#' @export
compute_delta_i <- function(I_curve, radii = NULL, smooth = FALSE,
                            smooth_window = 5) {
    # Handle all-NA case
    if (all(is.na(I_curve))) {
        return(list(
            delta_I = NA_real_,
            sign = NA_integer_,
            I_max = NA_real_,
            I_min = NA_real_,
            I_short = NA_real_,
            I_long = NA_real_,
            argmax = NA_integer_
        ))
    }

    # Smooth if requested
    if (smooth && length(I_curve) >= smooth_window) {
        I_curve <- .moving_average(I_curve, smooth_window)
    }

    # Compute statistics
    I_max <- max(I_curve, na.rm = TRUE)
    I_min <- min(I_curve, na.rm = TRUE)
    I_short <- I_curve[1]
    I_long <- I_curve[length(I_curve)]
    argmax <- which.max(I_curve)

    # Compute signed delta I
    sign <- ifelse(I_short >= I_long, 1L, -1L)
    delta_I <- sign * (I_max - I_min)

    list(
        delta_I = delta_I,
        sign = sign,
        I_max = I_max,
        I_min = I_min,
        I_short = I_short,
        I_long = I_long,
        argmax = argmax
    )
}

#' Simple Moving Average Smoothing
#'
#' @param x Numeric vector to smooth.
#' @param window Integer window size (must be odd for symmetric smoothing).
#'
#' @return Smoothed numeric vector.
#'
#' @keywords internal
.moving_average <- function(x, window) {
    n <- length(x)
    half <- floor(window / 2)
    result <- x

    for (i in (half + 1):(n - half)) {
        result[i] <- mean(x[(i - half):(i + half)], na.rm = TRUE)
    }

    result
}
