#' Pairwise Moran's I for Spatial Transcriptomics
#'
#' Compute pairwise Moran's I statistics between all genes in spatial
#' transcriptomics data.
#'
#' @param data A matrix or sparse matrix (dgCMatrix) of gene expression values.
#'   Rows are genes, columns are spots.
#' @param spot_coords A data frame or matrix with spot coordinates. Should have
#'   columns named 'row' and 'col', or be a 2-column matrix.
#' @param max_radius Maximum grid radius for computing distances. Default: 5.
#' @param platform Platform type: "visium" (default) or "old".
#' @param same_spot Whether to consider the same spot in computation. Default: FALSE.
#' @param mode Computation mode: "paired" (all pairs), "first" (pairs with first gene),
#'   or "single" (diagonal only). Default: "paired".
#' @param verbose Print progress messages. Default: TRUE.
#'
#' @return A matrix of pairwise Moran's I values. For mode="paired", returns a
#'   symmetric matrix where element [i,j] is the Moran's I between genes i and j.
#'
#' @examples
#' \dontrun{
#' # Create example data
#' data <- matrix(rnorm(1000 * 100), nrow = 1000, ncol = 100)
#' spot_coords <- data.frame(
#'   row = rep(1:10, each = 10),
#'   col = rep(1:10, times = 10)
#' )
#'
#' # Compute pairwise Moran's I
#' result <- pairwise_moran(data, spot_coords, max_radius = 3)
#' }
#'
#' @export
pairwise_moran <- function(data,
                           spot_coords,
                           max_radius = 5,
                           platform = c("visium", "old"),
                           same_spot = FALSE,
                           mode = c("paired", "first", "single"),
                           verbose = TRUE) {

  platform <- match.arg(platform)
  mode <- match.arg(mode)

  # Convert platform to integer
  platform_int <- if (platform == "visium") 0L else 1L

  # Convert mode to booleans
  paired_genes <- mode != "single"
  all_genes <- mode == "paired"

  # Parse spot coordinates
  if (is.data.frame(spot_coords)) {
    if ("row" %in% names(spot_coords) && "col" %in% names(spot_coords)) {
      spot_row <- as.integer(spot_coords$row)
      spot_col <- as.integer(spot_coords$col)
    } else {
      spot_row <- as.integer(spot_coords[[1]])
      spot_col <- as.integer(spot_coords[[2]])
    }
  } else if (is.matrix(spot_coords)) {
    spot_row <- as.integer(spot_coords[, 1])
    spot_col <- as.integer(spot_coords[, 2])
  } else {
    stop("spot_coords must be a data frame or matrix")
  }

  # Validate dimensions
  if (ncol(data) != length(spot_row)) {
    stop("Number of columns in data must match number of spots in spot_coords")
  }

  # Convert to matrix if needed and call C++ function
  if (methods::is(data, "sparseMatrix")) {
    # Sparse matrix path
    data <- methods::as(data, "dgCMatrix")
    result <- cpp_compute_moran_full_sparse(
      data,
      spot_row, spot_col,
      as.integer(max_radius),
      platform_int,
      same_spot,
      paired_genes,
      all_genes,
      verbose
    )
  } else {
    # Dense matrix path
    data <- as.matrix(data)
    result <- cpp_compute_moran_full(
      data,
      spot_row, spot_col,
      as.integer(max_radius),
      platform_int,
      same_spot,
      paired_genes,
      all_genes,
      verbose
    )
  }

  # Set row/column names if available
  if (!is.null(rownames(data))) {
    if (mode == "paired") {
      rownames(result) <- rownames(data)
      colnames(result) <- rownames(data)
    } else {
      names(result) <- rownames(data)
    }
  }

  return(result)
}


#' Compute Moran's I for a Single Gene Pair
#'
#' Compute Moran's I between two specific genes.
#'
#' @param x Numeric vector of expression values for gene 1.
#' @param y Numeric vector of expression values for gene 2.
#' @param W Weight matrix (spots x spots).
#' @param normalize Whether to z-normalize the input vectors. Default: TRUE.
#'
#' @return A single numeric value (Moran's I).
#'
#' @export
moran_I <- function(x, y, W, normalize = TRUE) {
  if (length(x) != length(y)) {
    stop("x and y must have the same length")
  }

  if (nrow(W) != length(x) || ncol(W) != length(x)) {
    stop("W must be a square matrix with dimensions matching x and y")
  }

  # Z-normalize if requested
  if (normalize) {
    x <- (x - mean(x)) / sd(x)
    y <- (y - mean(y)) / sd(y)
  }

  # Compute Moran's I: (x' W y) / sum(W)
  weight_sum <- sum(W)
  if (weight_sum == 0) {
    warning("Weight sum is zero")
    return(0)
  }

  result <- as.numeric(t(x) %*% W %*% y) / weight_sum
  return(result)
}


#' Create Spatial Weight Matrix
#'
#' Create a weight matrix based on spatial distances between spots.
#'
#' @param spot_coords A data frame or matrix with spot coordinates.
#' @param max_radius Maximum grid radius. Default: 5.
#' @param platform Platform type: "visium" or "old". Default: "visium".
#' @param same_spot Whether to include same-spot weights. Default: TRUE.
#'
#' @return A list with components:
#'   \item{W}{The weight matrix (spots x spots)}
#'   \item{weight_sum}{Sum of all weights}
#'
#' @export
create_weight_matrix <- function(spot_coords,
                                  max_radius = 5,
                                  platform = c("visium", "old"),
                                  same_spot = TRUE) {

  platform <- match.arg(platform)
  platform_int <- if (platform == "visium") 0L else 1L

  # Parse spot coordinates
  if (is.data.frame(spot_coords)) {
    if ("row" %in% names(spot_coords) && "col" %in% names(spot_coords)) {
      spot_row <- as.integer(spot_coords$row)
      spot_col <- as.integer(spot_coords$col)
    } else {
      spot_row <- as.integer(spot_coords[[1]])
      spot_col <- as.integer(spot_coords[[2]])
    }
  } else if (is.matrix(spot_coords)) {
    spot_row <- as.integer(spot_coords[, 1])
    spot_col <- as.integer(spot_coords[, 2])
  } else {
    stop("spot_coords must be a data frame or matrix")
  }

  # Create distance lookup
  distance <- cpp_create_distance(as.integer(max_radius), platform_int)
  if (!same_spot) {
    distance[1, 1] <- 0
  }

  # Create weight matrix
  result <- cpp_create_weight_matrix(
    spot_row, spot_col, distance,
    as.integer(max_radius), same_spot
  )

  return(result)
}


#' Parse Spot Names to Coordinates
#'
#' Convert spot names in format "ROWxCOL" to coordinate data frame.
#'
#' @param spot_names Character vector of spot names (e.g., "5x10", "3x7").
#'
#' @return A data frame with columns 'row' and 'col'.
#'
#' @examples
#' spots <- c("5x10", "5x12", "6x11")
#' coords <- parse_spot_names(spots)
#'
#' @export
parse_spot_names <- function(spot_names) {
  parts <- strsplit(spot_names, "x")
  coords <- data.frame(
    row = as.integer(sapply(parts, `[`, 1)),
    col = as.integer(sapply(parts, `[`, 2))
  )
  return(coords)
}
