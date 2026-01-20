#' Create Binary Weight Matrix for Visium
#'
#' Creates a row-normalized binary spatial weight matrix where all spots
#' within the specified radius receive equal weight.
#'
#' @param coords Numeric matrix. Spot coordinates (n x 2).
#' @param radius Numeric. Distance threshold for neighbor definition.
#' @param include_self Logical. Include self-connections (default: FALSE).
#'
#' @return Sparse weight matrix (n x n), row-normalized so each row sums to 1.
#'
#' @details
#' For Visium data, binary weights are preferred because:
#' \itemize{
#'   \item Grid structure has regular spacing (~100 um center-to-center)
#'   \item Spot density is relatively low (~4K spots per sample)
#'   \item Distance-decay effects are less pronounced at this resolution
#' }
#'
#' The weight formula is:
#' \deqn{w_{ij} = \frac{1}{n_i} \text{ if } d(i,j) \leq r, \text{ else } 0}
#'
#' @examples
#' # Create simple coordinates
#' coords <- matrix(c(0, 0, 100, 0, 0, 100, 100, 100), ncol = 2, byrow = TRUE)
#' W <- create_weights_visium(coords, radius = 150)
#' Matrix::rowSums(W)  # Each row sums to 1
#'
#' @seealso \code{\link{create_ring_weights_visium}},
#'   \code{\link{create_directional_weights_visium}}
#'
#' @export
create_weights_visium <- function(coords, radius, include_self = FALSE) {
    create_binary_weights_cpp(as.matrix(coords), radius, include_self)
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
