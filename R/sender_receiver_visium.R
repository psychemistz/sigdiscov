#' Split Spots by Expression Level
#'
#' Divides spots into sender and receiver groups based on factor gene
#' expression level. Senders are high-expressing spots (top percentile),
#' receivers are the remaining spots.
#'
#' @param expr Numeric vector. Expression values for the factor gene across
#'   all spots.
#' @param percentile Numeric. Percentile threshold (default: 75, meaning top
#'   25 percent of expressors become senders).
#' @param use_nonzero Logical. Compute threshold using only non-zero values
#'   (default: TRUE). This avoids inflation by zero-expressing spots.
#' @param min_spots Integer. Minimum required spots per group (default: 30).
#'   Warns if either group is smaller.
#'
#' @return A list with components:
#' \describe{
#'   \item{sender_idx}{Integer vector. Indices of sender spots (1-based).}
#'   \item{receiver_idx}{Integer vector. Indices of receiver spots (1-based).}
#'   \item{threshold}{Numeric. Expression threshold used for splitting.}
#'   \item{percentile}{Numeric. Percentile used (echoed from input).}
#'   \item{n_senders}{Integer. Number of sender spots.}
#'   \item{n_receivers}{Integer. Number of receiver spots.}
#' }
#'
#' @details
#' For Visium data, sender/receiver definition is expression-based because
#' spots contain mixed cell populations (1-10 cells per spot) and we cannot
#' use cell type annotations.
#'
#' The default strategy:
#' \enumerate{
#'   \item Compute 75th percentile of non-zero expression values
#'   \item Spots with expression >= threshold become senders (top ~25 percent)
#'   \item Remaining spots become receivers
#' }
#'
#' Using non-zero values for threshold computation avoids the situation where
#' many zero-expressing spots would inflate the receiver group.
#'
#' @examples
#' # Expression with some zeros
#' expr <- c(0, 0, 1, 2, 3, 4, 5, 6, 7, 10)
#' split <- split_by_expression(expr, percentile = 75)
#' split$n_senders
#' split$threshold
#'
#' @seealso \code{\link{create_directional_weights_visium}}
#'
#' @export
split_by_expression <- function(expr, percentile = 75, use_nonzero = TRUE,
                                 min_spots = 30) {
    # Validate input
    if (!is.numeric(expr)) {
        stop("expr must be a numeric vector")
    }

    if (percentile < 0 || percentile > 100) {
        stop("percentile must be between 0 and 100")
    }

    # Compute threshold
    if (use_nonzero) {
        expr_for_threshold <- expr[expr > 0]
        if (length(expr_for_threshold) == 0) {
            stop("No non-zero expression values for this gene")
        }
    } else {
        expr_for_threshold <- expr
    }

    threshold <- stats::quantile(expr_for_threshold, percentile / 100)

    # Split spots
    sender_idx <- which(expr >= threshold)
    receiver_idx <- which(expr < threshold)

    # Validate group sizes
    if (length(sender_idx) < min_spots) {
        warning("Only ", length(sender_idx), " senders (minimum: ", min_spots, ")")
    }
    if (length(receiver_idx) < min_spots) {
        warning("Only ", length(receiver_idx), " receivers (minimum: ", min_spots, ")")
    }

    list(
        sender_idx = sender_idx,
        receiver_idx = receiver_idx,
        threshold = as.numeric(threshold),
        percentile = percentile,
        n_senders = length(sender_idx),
        n_receivers = length(receiver_idx)
    )
}
