#' Create Gaussian Weight Matrix for Single-Cell Data
#'
#' Matches Python implementation: sigma = radius / 3
#'
#' @param sender_coords Sender cell coordinates (n x 2 matrix)
#' @param receiver_coords Receiver cell coordinates (m x 2 matrix)
#' @param radius Outer radius (micrometers)
#' @param inner_radius Inner radius for ring (default: 0)
#' @param sigma Gaussian sigma (default: radius/3)
#' @param min_weight Minimum weight threshold (default: 1e-6)
#' @return Sparse weight matrix (n x m), row-normalized
#' @export
create_weights_sc <- function(sender_coords, receiver_coords,
                               radius, inner_radius = 0, sigma = NULL,
                               min_weight = 1e-6) {
    if (is.null(sigma)) {
        sigma <- radius / 3  # Match Python default
    }

    create_gaussian_weights_cpp(
        as.matrix(sender_coords),
        as.matrix(receiver_coords),
        radius,
        inner_radius,
        sigma,
        min_weight
    )
}

#' Create Gaussian Ring Weight Matrix
#'
#' Weight matrix for annular region between inner and outer radius.
#'
#' @param coords Cell coordinates (n x 2 matrix)
#' @param outer_radius Outer radius (micrometers)
#' @param inner_radius Inner radius (micrometers)
#' @param sigma Gaussian sigma (default: outer_radius/3)
#' @return Sparse weight matrix (n x n), row-normalized
#' @export
create_ring_weights_sc <- function(coords, outer_radius, inner_radius,
                                    sigma = NULL) {
    if (is.null(sigma)) {
        sigma <- outer_radius / 3
    }

    create_gaussian_ring_weights_cpp(
        as.matrix(coords),
        outer_radius,
        inner_radius,
        sigma
    )
}
