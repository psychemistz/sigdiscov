#' Define Senders and Receivers by Cell Type
#'
#' Get cell indices for sender and receiver cell types.
#'
#' @param cell_types Character vector of cell type annotations
#' @param sender_type Sender cell type name
#' @param receiver_type Receiver cell type name
#' @param min_cells Minimum cells required (default: 10)
#' @return List with sender_type_idx, receiver_idx, counts
#' @export
split_by_celltype <- function(cell_types, sender_type, receiver_type,
                               min_cells = 10) {
    sender_type_idx <- which(cell_types == sender_type)
    receiver_idx <- which(cell_types == receiver_type)

    if (length(sender_type_idx) < min_cells) {
        stop("Insufficient cells of type '", sender_type, "': ",
             length(sender_type_idx), " (minimum: ", min_cells, ")")
    }
    if (length(receiver_idx) < min_cells) {
        stop("Insufficient cells of type '", receiver_type, "': ",
             length(receiver_idx), " (minimum: ", min_cells, ")")
    }

    list(
        sender_type_idx = sender_type_idx,
        receiver_idx = receiver_idx,
        n_sender_type = length(sender_type_idx),
        n_receivers = length(receiver_idx),
        sender_type = sender_type,
        receiver_type = receiver_type
    )
}

#' Filter Senders by Expression (Within Cell Type)
#'
#' Apply expression quantile filter to keep high-expression cells.
#' Matches Python: keep cells with expr > quantile(expr, quantile_prob)
#'
#' @param expr Expression vector (for cells of sender type)
#' @param quantile_prob Quantile threshold (default: 0.25 = keep top 75 percent)
#' @return List with high_expr_idx (relative to input), threshold, n_high
#' @export
filter_high_expression <- function(expr, quantile_prob = 0.25) {
    threshold <- quantile(expr, quantile_prob)
    high_expr_idx <- which(expr > threshold)

    list(
        high_expr_idx = high_expr_idx,
        threshold = as.numeric(threshold),
        n_high = length(high_expr_idx),
        quantile_prob = quantile_prob
    )
}

#' Combined Sender/Receiver Definition
#'
#' Define senders (cell type + expression filter) and receivers (cell type only).
#' Matches Python genomewide_interaction implementation.
#'
#' @param cell_types Character vector of cell type annotations
#' @param sender_type Sender cell type name
#' @param receiver_type Receiver cell type name
#' @param factor_expr RAW factor expression vector (all cells). If NULL, no
#'   expression filter is applied.
#' @param quantile_prob Quantile for expression filter (default: 0.25 = top 75 percent)
#' @param min_cells Minimum cells (default: 10)
#' @return List with sender_idx, receiver_idx, threshold, counts, etc.
#' @export
define_sender_receiver_sc <- function(cell_types, sender_type, receiver_type,
                                       factor_expr = NULL, quantile_prob = 0.25,
                                       min_cells = 10) {

    # Cell type split
    ct_split <- split_by_celltype(cell_types, sender_type, receiver_type, min_cells)

    # Apply expression filter within sender type
    if (!is.null(factor_expr)) {
        sender_expr <- factor_expr[ct_split$sender_type_idx]
        expr_filter <- filter_high_expression(sender_expr, quantile_prob)

        sender_idx <- ct_split$sender_type_idx[expr_filter$high_expr_idx]
        threshold <- expr_filter$threshold
    } else {
        sender_idx <- ct_split$sender_type_idx
        threshold <- NA_real_
    }

    if (length(sender_idx) < min_cells) {
        stop("Insufficient high-expression senders: ", length(sender_idx),
             " (minimum: ", min_cells, ")")
    }

    list(
        sender_idx = sender_idx,
        receiver_idx = ct_split$receiver_idx,
        sender_type_idx = ct_split$sender_type_idx,
        n_senders = length(sender_idx),
        n_receivers = ct_split$n_receivers,
        n_sender_type = ct_split$n_sender_type,
        sender_type = sender_type,
        receiver_type = receiver_type,
        threshold = threshold,
        quantile_prob = quantile_prob
    )
}

#' Get All Valid Cell Type Pairs
#'
#' Generate all sender-receiver cell type pairs with sufficient cells.
#'
#' @param cell_types Character vector of cell type annotations
#' @param min_cells Minimum cells per type (default: 10)
#' @param exclude_self Exclude same-type pairs (default: TRUE)
#' @return Data frame with sender and receiver columns
#' @export
get_celltype_pairs <- function(cell_types, min_cells = 10, exclude_self = TRUE) {
    type_counts <- table(cell_types)
    valid_types <- names(type_counts[type_counts >= min_cells])

    pairs <- expand.grid(sender = valid_types, receiver = valid_types,
                         stringsAsFactors = FALSE)

    if (exclude_self) {
        pairs <- pairs[pairs$sender != pairs$receiver, ]
    }

    pairs
}
