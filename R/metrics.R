#' Compute Bivariate Moran's I from Pre-computed Spatial Lag
#'
#' Calculates Moran's I statistic given a factor expression vector and
#' pre-computed spatial lag of gene expression.
#'
#' @param z_f Numeric vector. Standardized factor expression (length n).
#' @param lag_g Numeric vector. Spatial lag of gene expression (W * z_g),
#'   same length as z_f.
#'
#' @return Scalar Moran's I value (unbounded).
#'
#' @details
#' The bivariate Moran's I is computed as:
#' \deqn{I = \frac{z_f' \cdot lag_g}{n}}
#'
#' where:
#' \itemize{
#'   \item \code{z_f} is the standardized factor expression
#'   \item \code{lag_g} is the pre-computed spatial lag (W * z_g)
#'   \item \code{n} is the number of observations
#' }
#'
#' @family core
#'
#' @examples
#' # Simple example
#' z_f <- c(1, -1, 1, -1)
#' lag_g <- c(0.5, -0.5, 0.5, -0.5)
#' compute_moran_from_lag(z_f, lag_g)
#'
#' @export
compute_moran_from_lag <- function(z_f, lag_g) {
    compute_moran_from_lag_cpp(z_f, lag_g)
}

#' Compute I_ND (Cosine Similarity) from Pre-computed Spatial Lag
#'
#' Calculates the normalized directional Moran's I (I_ND) using cosine
#' similarity between factor expression and spatial lag of gene expression.
#'
#' @param z_f Numeric vector. Standardized factor expression (length n).
#' @param lag_g Numeric vector. Spatial lag of gene expression (W * z_g).
#'
#' @return Scalar I_ND value bounded between -1 and 1, or NA if norms are too small.
#'
#' @details
#' The I_ND is computed as cosine similarity:
#' \deqn{I_{ND} = \frac{z_f' \cdot lag_g}{\|z_f\| \cdot \|lag_g\|}}
#'
#' This metric is bounded between -1 and 1, making it interpretable as a
#' correlation coefficient.
#'
#' @family core
#'
#' @examples
#' z_f <- rnorm(100)
#' lag_g <- rnorm(100)
#' compute_ind_from_lag(z_f, lag_g)
#'
#' @export
compute_ind_from_lag <- function(z_f, lag_g) {
    compute_ind_from_lag_cpp(z_f, lag_g)
}

#' Batch Compute Metrics for All Genes
#'
#' Efficiently computes Moran's I or I_ND for one factor against all genes
#' using pre-computed spatial lag matrix.
#'
#' @param z_f Numeric vector. Standardized factor expression (length n).
#' @param lag_G Numeric matrix. Spatial lag matrix (n x n_genes), where each
#'   column is W * z_g for gene g.
#' @param metric Character. Either "moran" or "ind".
#'
#' @return Numeric vector of metric values (length n_genes).
#'
#' @family core
#'
#' @examples
#' n <- 100
#' n_genes <- 50
#' z_f <- rnorm(n)
#' lag_G <- matrix(rnorm(n * n_genes), n, n_genes)
#' compute_metric_batch(z_f, lag_G, "ind")
#'
#' @export
compute_metric_batch <- function(z_f, lag_G, metric = c("moran", "ind")) {
    metric <- match.arg(metric)
    compute_metric_batch_cpp(z_f, lag_G, metric)
}

#' Compute Spatial Lag (W * z)
#'
#' Computes the spatial lag of an expression vector using a weight matrix.
#'
#' @param W Sparse weight matrix (n x m), typically row-normalized.
#' @param z Numeric vector. Expression values (length m).
#'
#' @return Numeric vector. Spatial lag (length n).
#'
#' @details
#' The spatial lag represents the weighted average of neighboring values:
#' \deqn{lag_i = \sum_j w_{ij} z_j}
#'
#' @family core
#'
#' @examples
#' library(Matrix)
#' W <- sparseMatrix(i = c(1,1,2,2), j = c(2,3,1,3), x = c(0.5,0.5,0.5,0.5))
#' z <- c(1, 2, 3)
#' compute_spatial_lag(W, z)
#'
#' @export
compute_spatial_lag <- function(W, z) {
    as.vector(W %*% z)
}

#' Batch Compute Spatial Lags (W * Z)
#'
#' Computes spatial lags for all genes at once using matrix multiplication.
#'
#' @param W Sparse weight matrix (n x m), row-normalized.
#' @param Z Numeric matrix. Expression matrix (m x n_genes), genes in columns.
#'
#' @return Numeric matrix. Spatial lag matrix (n x n_genes).
#'
#' @family core
#'
#' @examples
#' library(Matrix)
#' n <- 50
#' m <- 50
#' n_genes <- 100
#' W <- rsparsematrix(n, m, density = 0.1)
#' W <- W / rowSums(abs(W))
#' Z <- matrix(rnorm(m * n_genes), m, n_genes)
#' lag_G <- compute_spatial_lag_batch(W, Z)
#'
#' @export
compute_spatial_lag_batch <- function(W, Z) {
    as.matrix(W %*% Z)
}
