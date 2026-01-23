#' Create Weight Matrix for Visium
#'
#' Creates a spatial weight matrix for Visium data. Supports circular RBF weights
#' (SpaCET-compatible, default), grid-based Gaussian weights, and binary weights.
#'
#' @param coords Spot coordinates. For type="circular", should have x/y physical
#'   coordinates (e.g., micrometers). For type="grid", should have 'row' and 'col'
#'   grid coordinates.
#' @param radius Numeric. Distance cutoff in coordinate units. Default: 200.
#' @param type Character. Weight type: "circular" (default, SpaCET-compatible),
#'   "grid" (legacy hexagonal), or "binary" (row-normalized).
#' @param sigma Numeric. RBF kernel bandwidth for circular weights. Default: 100.
#' @param include_self Logical. Include self-connections (default: FALSE).
#' @param sparse Logical. Return sparse matrix (default: TRUE for circular).
#' @param platform Character. Platform type for grid weights: "visium" or "old".
#'
#' @return A list with components:
#'   \item{W}{Weight matrix (sparse or dense)}
#'   \item{weight_sum}{Sum of all weights (S0 for Moran's I normalization)}
#'
#' @details
#' Three weight types are supported:
#'
#' **Circular (default)**: Uses Euclidean distance with RBF kernel and circular
#' cutoff. Matches SpaCET's calWeights behavior.
#' \deqn{w_{ij} = \exp(-d_{ij}^2 / (2\sigma^2)) \text{ if } d_{ij} \leq r}
#'
#' **Grid**: Uses hexagonal grid coordinates with Gaussian distance decay.
#' Legacy method for backward compatibility.
#'
#' **Binary**: All spots within radius receive equal weight (row-normalized).
#' \eqn{w_{ij} = 1/n_i} if \eqn{d(i,j) \le r}, else 0.
#'
#' @examples
#' \dontrun{
#' # Circular weights (default, SpaCET-compatible)
#' coords <- data.frame(x = spot_x_um, y = spot_y_um)
#' W_result <- create_weights_visium(coords, radius = 200, sigma = 100)
#'
#' # Grid-based weights (legacy)
#' grid_coords <- parse_spot_names(colnames(data))
#' W_result <- create_weights_visium(grid_coords, radius = 300, type = "grid")
#' }
#'
#' @seealso \code{\link{pairwise_moran}}, \code{\link{create_ring_weights_visium}}
#'
#' @export
create_weights_visium <- function(coords, radius = 200,
                                   type = c("circular", "grid", "binary"),
                                   sigma = 100,
                                   include_self = FALSE,
                                   sparse = TRUE,
                                   platform = c("visium", "old")) {
    type <- match.arg(type)
    platform <- match.arg(platform)

    # Circular RBF weights (default, SpaCET-compatible)
    if (type == "circular") {
        # Parse coordinates
        if (is.data.frame(coords)) {
            if ("coordinate_x_um" %in% names(coords) && "coordinate_y_um" %in% names(coords)) {
                coord_mat <- as.matrix(coords[, c("coordinate_x_um", "coordinate_y_um")])
            } else if ("x" %in% names(coords) && "y" %in% names(coords)) {
                coord_mat <- as.matrix(coords[, c("x", "y")])
            } else {
                coord_mat <- as.matrix(coords[, 1:2])
            }
        } else {
            coord_mat <- as.matrix(coords)
        }

        if (ncol(coord_mat) != 2) {
            stop("coords must have exactly 2 columns (x, y coordinates)")
        }

        if (sparse) {
            result <- create_circular_weights_cpp(coord_mat, radius, sigma, include_self)
        } else {
            result <- create_circular_weights_dense_cpp(coord_mat, radius, sigma, include_self)
        }
        return(result)
    }

    # Binary weights (row-normalized)
    if (type == "binary") {
        W <- create_binary_weights_cpp(as.matrix(coords), radius, include_self)
        weight_sum <- sum(W)
        return(list(W = W, weight_sum = weight_sum))
    }

    # Grid-based Gaussian weights (legacy)
    if (is.data.frame(coords)) {
        if ("row" %in% names(coords) && "col" %in% names(coords)) {
            spot_row <- as.integer(coords$row)
            spot_col <- as.integer(coords$col)
        } else {
            spot_row <- as.integer(coords[[1]])
            spot_col <- as.integer(coords[[2]])
        }
    } else if (is.matrix(coords)) {
        spot_row <- as.integer(coords[, 1])
        spot_col <- as.integer(coords[, 2])
    } else {
        stop("coords must be a data frame or matrix with row/col coordinates")
    }

    max_radius <- as.integer(ceiling(radius / 100))
    if (max_radius < 1) max_radius <- 1L

    platform_int <- if (platform == "visium") 0L else 1L
    distance <- cpp_create_distance(max_radius, platform_int)

    if (!include_self) {
        distance[1, 1] <- 0
    }

    W_result <- cpp_create_weight_matrix_sparse(spot_row, spot_col, distance,
                                                 max_radius, include_self)

    return(list(W = W_result$W, weight_sum = W_result$weight_sum))
}

#' Create Ring Weight Matrix for Visium
#'
#' Creates a weight matrix where only spots in an annulus (ring) between
#' inner and outer radii receive weight. Useful for excluding immediate
#' neighbors or analyzing specific distance bands.
#'
#' @param coords Numeric matrix. Spot coordinates (n x 2).
#' @param inner_radius Numeric. Inner radius (exclusive).
#' @param outer_radius Numeric. Outer radius (inclusive).
#'
#' @return Sparse ring weight matrix (n x n), row-normalized.
#'
#' @details
#' Weight is assigned to spots where inner_radius < d(i,j) <= outer_radius.
#' This is useful for:
#' \itemize{
#'   \item Analyzing paracrine signaling at specific distance ranges
#'   \item Excluding immediate neighbors to focus on longer-range effects
#'   \item Computing spatial correlation at multiple distance bands
#' }
#'
#' @examples
#' coords <- matrix(runif(200), 100, 2) * 1000
#' W_ring <- create_ring_weights_visium(coords, inner_radius = 100, outer_radius = 200)
#'
#' @seealso \code{\link{create_weights_visium}}
#'
#' @export
create_ring_weights_visium <- function(coords, inner_radius, outer_radius) {
    if (inner_radius >= outer_radius) {
        stop("inner_radius must be less than outer_radius")
    }
    create_ring_weights_cpp(as.matrix(coords), inner_radius, outer_radius)
}

#' Create Directional Weight Matrix (Sender to Receiver)
#'
#' Creates a weight matrix for directional analysis where senders and
#' receivers are different sets of spots.
#'
#' @param sender_coords Numeric matrix. Sender spot coordinates (n_s x 2).
#' @param receiver_coords Numeric matrix. Receiver spot coordinates (n_r x 2).
#' @param radius Numeric. Distance threshold.
#'
#' @return Sparse weight matrix (n_senders x n_receivers), row-normalized.
#'
#' @details
#' For directional mode, the weight matrix is not square. Each sender spot
#' gets connections to receiver spots within the radius, with row normalization
#' so that weights sum to 1 per sender.
#'
#' This is the key function for expression-based sender/receiver analysis
#' in Visium data.
#'
#' @examples
#' # Senders at top, receivers at bottom
#' sender_coords <- matrix(c(0, 500, 100, 500, 200, 500), ncol = 2, byrow = TRUE)
#' receiver_coords <- matrix(c(0, 0, 100, 0, 200, 0), ncol = 2, byrow = TRUE)
#' W <- create_directional_weights_visium(sender_coords, receiver_coords, radius = 600)
#'
#' @seealso \code{\link{split_by_expression}}
#'
#' @export
create_directional_weights_visium <- function(sender_coords, receiver_coords,
                                               radius) {
    create_directional_weights_cpp(
        as.matrix(sender_coords),
        as.matrix(receiver_coords),
        radius
    )
}


#' @rdname create_weights_visium
#' @keywords internal
create_circular_weights_visium <- function(coords, radius, sigma = 100,
                                            include_self = FALSE, sparse = TRUE) {
    create_weights_visium(coords, radius, type = "circular", sigma = sigma,
                          include_self = include_self, sparse = sparse)
}


#' @rdname pairwise_moran
#' @keywords internal
pairwise_moran_circular <- function(data, coords, radius, sigma = 100,
                                     include_self = FALSE, mode = "paired",
                                     sparse = TRUE, verbose = TRUE) {
    # Create circular weight matrix (sparse by default for memory efficiency)
    if (verbose) message("Creating circular RBF weight matrix...")
    W_result <- create_weights_visium(coords, radius, type = "circular",
                                       sigma = sigma, include_self = include_self,
                                       sparse = sparse)

    if (W_result$weight_sum == 0) {
        stop("Weight matrix has no non-zero weights. Try increasing radius.")
    }

    if (verbose) {
        message(sprintf("  Weight sum (S0): %.4f", W_result$weight_sum))
        if (sparse) {
            message(sprintf("  W non-zeros: %d", Matrix::nnzero(W_result$W)))
        } else {
            message(sprintf("  Non-zero weights: %d", sum(W_result$W > 0)))
        }
    }

    # Compute Moran's I using custom weights (handles sparse/dense automatically)
    result <- pairwise_moran_custom(data, W_result$W, normalize = FALSE,
                                     mode = mode, verbose = verbose)
    result$W <- W_result$W
    result
}
