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

#' Standardize Vector (Z-score)
#'
#' Standardizes a numeric vector to have mean 0 and standard deviation 1.
#'
#' @param x Numeric vector.
#' @param na.rm Logical. Remove NAs before computing mean/sd (default: TRUE).
#'
#' @return Standardized numeric vector.
#'
#' @details
#' Computes: \code{(x - mean(x)) / sd(x)}
#'
#' A small constant (1e-10) is added to the denominator to avoid
#' division by zero for constant vectors.
#'
#' @examples
#' x <- c(1, 2, 3, 4, 5)
#' standardize(x)
#' # Returns approximately: -1.26, -0.63, 0, 0.63, 1.26
#'
#' @importFrom stats sd
#' @export
standardize <- function(x, na.rm = TRUE) {
    m <- mean(x, na.rm = na.rm)
    s <- sd(x, na.rm = na.rm)
    (x - m) / (s + 1e-10)
}

#' Standardize Matrix (Gene-wise)
#'
#' Standardizes each row (gene) of a matrix to have mean 0 and sd 1.
#'
#' @param X Numeric matrix (genes x observations).
#'
#' @return Gene-wise standardized matrix with same dimensions.
#'
#' @details
#' Each row is standardized independently. This is the standard
#' preprocessing for spatial correlation analysis.
#'
#' @examples
#' X <- matrix(rnorm(100), 10, 10)
#' X_norm <- standardize_matrix(X)
#' rowMeans(X_norm)  # All approximately 0
#' apply(X_norm, 1, sd)  # All approximately 1
#'
#' @export
standardize_matrix <- function(X) {
    t(scale(t(X)))
}

#' Adjust P-values for Multiple Testing
#'
#' Wrapper for \code{\link[stats]{p.adjust}} for multiple testing correction.
#'
#' @param p Numeric vector of p-values.
#' @param method Character. Adjustment method (default: "BH" for
#'   Benjamini-Hochberg).
#'
#' @return Numeric vector of adjusted p-values.
#'
#' @details
#' Available methods include:
#' \itemize{
#'   \item "BH" or "fdr": Benjamini-Hochberg (controls FDR)
#'   \item "bonferroni": Bonferroni correction (controls FWER)
#'   \item "holm": Holm-Bonferroni
#'   \item "hochberg": Hochberg
#'   \item "hommel": Hommel
#'   \item "BY": Benjamini-Yekutieli
#' }
#'
#' @examples
#' p <- c(0.001, 0.01, 0.05, 0.1, 0.5)
#' adjust_pvalues(p, "BH")
#'
#' @export
adjust_pvalues <- function(p, method = "BH") {
    stats::p.adjust(p, method = method)
}

#' Row-Normalize Sparse Matrix
#'
#' Normalizes each row of a sparse matrix to sum to 1 (or 0 for empty rows).
#'
#' @param W Sparse matrix (Matrix package format).
#'
#' @return Row-normalized sparse matrix.
#'
#' @details
#' Row normalization ensures that spatial weights represent proper
#' averaging: each row sums to 1 (for non-isolated observations) or 0
#' (for isolated observations with no neighbors).
#'
#' @examples
#' library(Matrix)
#' W <- sparseMatrix(i = c(1,1,2,2,3,3), j = c(2,3,1,3,1,2),
#'                   x = c(1,2,1,1,3,1), dims = c(3,3))
#' W_norm <- sparse_row_normalize(W)
#' Matrix::rowSums(W_norm)  # All 1s
#'
#' @importFrom Matrix rowSums
#' @export
sparse_row_normalize <- function(W) {
    sparse_row_normalize_cpp(W)
}
