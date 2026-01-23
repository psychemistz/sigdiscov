#' Create Weight Matrix for Visium
#'
#' Creates a spatial weight matrix for Visium data. Supports both Gaussian
#' distance decay weights (consistent with pairwise_moran) and binary weights.
#'
#' @param coords Spot coordinates. For type="gaussian", should be a data frame
#'   with 'row' and 'col' columns (grid coordinates). For type="binary", can be
#'   any n x 2 matrix of physical coordinates.
#' @param radius Numeric. For type="binary", distance threshold in coordinate units.
#'   For type="gaussian", converted to max_radius in grid units (radius / 100).
#' @param type Character. Weight type: "gaussian" (default) or "binary".
#' @param include_self Logical. Include self-connections (default: FALSE).
#' @param platform Character. Platform type for gaussian weights: "visium" (default)
#'   or "old". Only used when type = "gaussian".
#'
#' @return A list with components:
#'   \item{W}{Sparse weight matrix (n x n)}
#'   \item{weight_sum}{Sum of all weights (for normalization)}
#'
#' @details
#' Two weight types are supported:
#'
#' **Gaussian (default)**: Uses the same grid-based Gaussian distance decay as
#' \code{pairwise_moran()}. Requires grid coordinates (row, col integers).
#' Uses: \eqn{w_{ij} = exp(-d_{ij}^2 / (2 * 100^2))} where d is computed using
#' Visium hexagonal grid geometry. This is recommended for consistency with
#' pairwise Moran's I analysis.
#'
#' **Binary**: All spots within radius receive equal weight (row-normalized).
#' \eqn{w_{ij} = 1/n_i} if \eqn{d(i,j) \le r}, else 0.
#'
#' @examples
#' \dontrun{
#' # Parse grid coordinates from spot names
#' spot_coords <- parse_spot_names(colnames(data))
#'
#' # Gaussian weights (default, matches pairwise_moran)
#' W_result <- create_weights_visium(spot_coords, radius = 300, type = "gaussian")
#'
#' # Binary weights (legacy)
#' W_result <- create_weights_visium(spot_coords, radius = 300, type = "binary")
#' }
#'
#' @seealso \code{\link{pairwise_moran}}, \code{\link{create_ring_weights_visium}},
#'   \code{\link{create_directional_weights_visium}}
#'
#' @export
create_weights_visium <- function(coords, radius, type = c("gaussian", "binary"),
                                   include_self = FALSE, platform = c("visium", "old")) {
    type <- match.arg(type)
    platform <- match.arg(platform)

    if (type == "binary") {
        # Legacy binary weights (row-normalized)
        W <- create_binary_weights_cpp(as.matrix(coords), radius, include_self)
        weight_sum <- sum(W)
        return(list(W = W, weight_sum = weight_sum))
    }

    # Gaussian distance decay weights (same as pairwise_moran)
    # Parse coordinates - need row/col integers for grid-based computation
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

    # Convert radius to max_radius (grid units)
    # Visium center-to-center distance is ~100 um
    max_radius <- as.integer(ceiling(radius / 100))
    if (max_radius < 1) max_radius <- 1L

    # Use the same functions as pairwise_moran
    platform_int <- if (platform == "visium") 0L else 1L
    distance <- cpp_create_distance(max_radius, platform_int)

    if (!include_self) {
        distance[1, 1] <- 0
    }

    # Create weight matrix (same as pairwise_moran uses internally)
    # Use sparse version for memory efficiency and compatibility with permutation tests
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


#' Create Circular RBF Weight Matrix for Visium
#'
#' Creates a spatial weight matrix using RBF kernel with circular (Euclidean)
#' distance cutoff. This matches the SpaCET calWeights behavior for spatial
#' correlation analysis.
#'
#' @param coords Spot coordinates. Can be:
#'   - A matrix with columns for x/y coordinates (in physical units like micrometers)
#'   - A data frame with columns 'x' and 'y', or 'coordinate_x_um' and 'coordinate_y_um'
#' @param radius Distance cutoff in the same units as coords. Only spots within
#'   this radius are connected. For Visium, typical values are 100-500 um.
#' @param sigma RBF kernel bandwidth. Default: 100. Controls distance decay rate.
#' @param include_self Include self-connections (diagonal). Default: FALSE.
#' @param sparse Return sparse matrix (TRUE) or dense matrix (FALSE). Default: TRUE.
#'
#' @return A list with components:
#'   \item{W}{Weight matrix (sparse or dense)}
#'   \item{weight_sum}{Sum of all weights (S0 for Moran's I normalization)}
#'
#' @details
#' The weight between spots i and j is computed as:
#' \deqn{w_{ij} = \exp\left(-\frac{d_{ij}^2}{2\sigma^2}\right) \text{ if } d_{ij} \leq r}
#'
#' where \eqn{d_{ij}} is the Euclidean distance.
#'
#' This is the recommended weight function for spatial correlation analysis as it:
#' - Uses circular (Euclidean) distance instead of rectangular grid
#' - Applies continuous distance decay via RBF kernel
#' - Is consistent with SpaCET's approach
#'
#' @examples
#' \dontrun{
#' # Load Visium data with physical coordinates
#' coords <- data.frame(
#'   x = visium_data$coords$coordinate_x_um,
#'   y = visium_data$coords$coordinate_y_um
#' )
#'
#' # Create weight matrix with 200um radius
#' W_result <- create_circular_weights_visium(coords, radius = 200, sigma = 100)
#'
#' # Use for pairwise Moran's I
#' result <- pairwise_moran_custom(data, W_result$W)
#' }
#'
#' @seealso \code{\link{create_weights_visium}}, \code{\link{pairwise_moran_custom}}
#'
#' @export
create_circular_weights_visium <- function(coords,
                                            radius,
                                            sigma = 100,
                                            include_self = FALSE,
                                            sparse = TRUE) {
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

    # Call C++ function
    if (sparse) {
        result <- create_circular_weights_cpp(coord_mat, radius, sigma, include_self)
    } else {
        result <- create_circular_weights_dense_cpp(coord_mat, radius, sigma, include_self)
    }

    result
}


#' Compute Pairwise Moran's I with Circular Weights
#'
#' Convenience function that creates circular RBF weights and computes
#' pairwise Moran's I in one step. Matches the SpaCET workflow.
#'
#' @param data Gene expression matrix (genes x spots). Will be z-normalized.
#' @param coords Spot coordinates (see \code{\link{create_circular_weights_visium}}).
#' @param radius Distance cutoff for weight matrix.
#' @param sigma RBF kernel bandwidth. Default: 100.
#' @param include_self Include self-connections. Default: FALSE.
#' @param mode Calculation mode: "paired" (all pairs) or "single" (univariate).
#' @param verbose Print progress. Default: TRUE.
#'
#' @return A list containing:
#'   \item{moran}{Moran's I matrix (if mode = "paired") or vector (if mode = "single")}
#'   \item{gene_names}{Names of genes}
#'   \item{weight_sum}{Sum of weights used (S0)}
#'   \item{W}{The weight matrix used}
#'
#' @examples
#' \dontrun{
#' # One-step pairwise Moran's I with circular weights
#' result <- pairwise_moran_circular(
#'   data_vst,
#'   coords,
#'   radius = 200,
#'   sigma = 100
#' )
#' }
#'
#' @seealso \code{\link{create_circular_weights_visium}}, \code{\link{pairwise_moran_custom}}
#'
#' @export
pairwise_moran_circular <- function(data,
                                     coords,
                                     radius,
                                     sigma = 100,
                                     include_self = FALSE,
                                     mode = "paired",
                                     verbose = TRUE) {
    # Create circular weight matrix
    if (verbose) message("Creating circular RBF weight matrix...")
    W_result <- create_circular_weights_visium(
        coords, radius, sigma, include_self, sparse = FALSE
    )

    if (W_result$weight_sum == 0) {
        stop("Weight matrix has no non-zero weights. Try increasing radius.")
    }

    if (verbose) {
        message(sprintf("  Weight sum (S0): %.4f", W_result$weight_sum))
        message(sprintf("  Non-zero weights: %d", sum(W_result$W > 0)))
    }

    # Compute Moran's I using custom weights
    result <- pairwise_moran_custom(
        data, W_result$W,
        normalize = FALSE,
        mode = mode,
        verbose = verbose
    )

    # Add weight matrix to result
    result$W <- W_result$W

    result
}
