#' @title Permutation Tests for Spatial Correlation
#' @description Functions for permutation-based significance testing of Moran's I
#' and I_ND spatial correlation metrics.
#' @name permutation
#' @importFrom stats p.adjust quantile
NULL


#' Permutation Test for Spatial Correlation
#'
#' Performs permutation test for Moran's I or I_ND to assess statistical
#' significance of spatial association between factor and gene expression.
#'
#' @param expr_matrix Gene expression matrix (genes x cells/spots).
#' @param cell_meta Data frame with cell metadata (for single-cell mode).
#'   Must contain: cell_id, x, y, cell_type.
#' @param spot_coords Coordinate matrix or data frame for Visium (spots x 2).
#' @param factor_gene Name or index of factor gene.
#' @param sender_type Sender cell type (for single-cell directional mode).
#' @param receiver_type Receiver cell type (for single-cell directional mode).
#' @param radius Distance threshold for weight matrix. Default: 50 for single-cell, 200 for Visium.
#' @param mode "bivariate" (all cells/spots) or "directional" (sender to receiver).
#' @param metric "moran" (Moran's I) or "ind" (I_ND, default).
#' @param n_perm Number of permutations. Default: 999.
#' @param seed Random seed for reproducibility.
#' @param adjust_method P-value adjustment method. Default: "BH" (Benjamini-Hochberg).
#' @param verbose Print progress. Default: TRUE.
#'
#' @return A data.frame with columns:
#'   \item{gene}{Gene name}
#'   \item{I_obs}{Observed statistic (Moran's I or I_ND)}
#'   \item{p_value}{Permutation p-value}
#'   \item{p_adj}{Adjusted p-value (FDR)}
#'   \item{z_score}{Standardized effect size}
#'
#' @details
#' \strong{Permutation Strategy:}
#' The test permutes the factor expression while keeping the spatial lag fixed.
#' This is efficient because:
#' \enumerate{
#'   \item Spatial lag (W * z_g) is computed once per gene
#'   \item Permutation involves only cheap vector operations
#'   \item Same permutation indices are reused across all genes
#' }
#'
#' \strong{P-value Calculation:}
#' \deqn{p = \frac{1 + \sum I(|I^{(b)}| \geq |I_{obs}|)}{B + 1}}
#'
#' The "+1" provides continuity correction and ensures p > 0.
#'
#' @examples
#' \dontrun{
#' # Visium: Bivariate Moran's I test
#' result <- permutation_test_spatial(
#'     expr_matrix = expr_vst,
#'     spot_coords = coords,
#'     factor_gene = "IL1B",
#'     mode = "bivariate",
#'     metric = "moran",
#'     n_perm = 999
#' )
#'
#' # CosMx: Directional I_ND test (cell-type specific)
#' result <- permutation_test_spatial(
#'     expr_matrix = cosmx_expr,
#'     cell_meta = meta,
#'     factor_gene = "IL1B",
#'     sender_type = "Macrophage",
#'     receiver_type = "Fibroblast",
#'     mode = "directional",
#'     metric = "ind",
#'     n_perm = 999
#' )
#'
#' # Significant genes (FDR < 0.05)
#' sig_genes <- result[result$p_adj < 0.05, ]
#' }
#'
#' @export
permutation_test_spatial <- function(expr_matrix,
                                     cell_meta = NULL,
                                     spot_coords = NULL,
                                     factor_gene,
                                     sender_type = NULL,
                                     receiver_type = NULL,
                                     radius = NULL,
                                     mode = c("directional", "bivariate"),
                                     metric = c("ind", "moran"),
                                     n_perm = 999,
                                     seed = NULL,
                                     adjust_method = "BH",
                                     verbose = TRUE) {

  mode <- match.arg(mode)
  metric <- match.arg(metric)

  # Convert sparse matrix if needed
  if (inherits(expr_matrix, "dgCMatrix") || inherits(expr_matrix, "dgTMatrix")) {
    expr_matrix <- as.matrix(expr_matrix)
  }
  if (!is.matrix(expr_matrix)) {
    expr_matrix <- as.matrix(expr_matrix)
  }

  gene_names <- rownames(expr_matrix)
  if (is.null(gene_names)) {
    gene_names <- paste0("Gene_", seq_len(nrow(expr_matrix)))
  }

  # Get factor index
  if (is.character(factor_gene)) {
    factor_idx <- which(gene_names == factor_gene)
    if (length(factor_idx) == 0) {
      stop("Factor gene '", factor_gene, "' not found")
    }
    factor_idx <- factor_idx[1] - 1L  # 0-based
  } else {
    factor_idx <- as.integer(factor_gene) - 1L
  }

  # Determine if single-cell or Visium mode
  is_singlecell <- !is.null(cell_meta) && !is.null(sender_type) && !is.null(receiver_type)

  if (is_singlecell) {
    # Single-cell mode (CosMx, Xenium, etc.)
    required_cols <- c("cell_id", "x", "y", "cell_type")
    missing_cols <- setdiff(required_cols, colnames(cell_meta))
    if (length(missing_cols) > 0) {
      stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
    }

    # Match cell order
    if (!all(colnames(expr_matrix) == cell_meta$cell_id)) {
      cell_meta <- cell_meta[match(colnames(expr_matrix), cell_meta$cell_id), ]
    }

    coords <- as.matrix(cell_meta[, c("x", "y")])
    cell_types <- as.character(cell_meta$cell_type)

    # Get sender/receiver indices (0-based)
    sender_idx <- which(cell_types == sender_type) - 1L
    receiver_idx <- which(cell_types == receiver_type) - 1L

    if (is.null(radius)) radius <- 50  # Default for single-cell: 50 um

    if (verbose) {
      message("Single-cell permutation test:")
      message("  Factor: ", gene_names[factor_idx + 1])
      message("  Triplet: ", sender_type, " -> ", receiver_type)
      message("  Senders: ", length(sender_idx), " cells")
      message("  Receivers: ", length(receiver_idx), " cells")
    }

  } else if (!is.null(spot_coords)) {
    # Visium mode
    coords <- as.matrix(spot_coords)
    if (ncol(coords) < 2) {
      stop("spot_coords must have at least 2 columns (x, y or row, col)")
    }
    coords <- coords[, 1:2]

    n_spots <- ncol(expr_matrix)

    if (mode == "directional") {
      # Expression-based sender/receiver split (top 25% = senders)
      factor_expr <- expr_matrix[factor_idx + 1, ]
      nonzero_expr <- factor_expr[factor_expr > 0]
      if (length(nonzero_expr) == 0) {
        stop("Factor gene has no non-zero expression")
      }
      threshold <- quantile(nonzero_expr, 0.75)
      sender_idx <- which(factor_expr >= threshold) - 1L
      receiver_idx <- which(factor_expr < threshold) - 1L
    } else {
      # Bivariate: all spots
      sender_idx <- seq(0, n_spots - 1)
      receiver_idx <- seq(0, n_spots - 1)
    }

    if (is.null(radius)) radius <- 200  # Default for Visium

    if (verbose) {
      message("Visium permutation test:")
      message("  Factor: ", gene_names[factor_idx + 1])
      message("  Mode: ", mode)
      message("  Senders: ", length(sender_idx), " spots")
      message("  Receivers: ", length(receiver_idx), " spots")
    }

  } else {
    stop("Must provide either cell_meta (single-cell) or spot_coords (Visium)")
  }

  if (length(sender_idx) < 10 || length(receiver_idx) < 10) {
    stop("Insufficient sender or receiver cells/spots for permutation test")
  }

  # Build weight matrix
  if (verbose) {
    message("  Building weight matrix (radius = ", radius, ")...")
  }

  W <- .build_weight_matrix_perm(
    coords = coords,
    sender_idx = sender_idx,
    receiver_idx = receiver_idx,
    r_outer = radius
  )

  # Run permutation test
  if (verbose) {
    message("  Running ", n_perm, " permutations on ", nrow(expr_matrix), " genes...")
  }

  seed_int <- if (is.null(seed)) -1L else as.integer(seed)

  result <- cpp_batch_permutation_test(
    expr_matrix = expr_matrix,
    gene_names = gene_names,
    factor_idx = factor_idx,
    W = W,
    sender_idx = as.integer(sender_idx),
    receiver_idx = as.integer(receiver_idx),
    n_perm = as.integer(n_perm),
    metric = if (metric == "moran") 0L else 1L,
    seed = seed_int,
    verbose = verbose
  )

  # Adjust p-values
  result$p_adj <- p.adjust(result$p_value, method = adjust_method)

  # Add metadata
  attr(result, "factor_gene") <- gene_names[factor_idx + 1]
  attr(result, "mode") <- mode
  attr(result, "metric") <- metric
  attr(result, "n_perm") <- n_perm
  attr(result, "radius") <- radius

  return(result)
}


#' Permutation Test for Single-Cell ST with Signature Computation
#'
#' Combines I_ND signature computation with permutation testing.
#' Computes I_ND at multiple radii and tests significance at the first radius.
#'
#' @param expr_matrix Gene expression matrix (genes x cells).
#' @param cell_meta Data frame with cell metadata.
#' @param factor_gene Name or index of factor gene.
#' @param sender_type Sender cell type.
#' @param receiver_type Receiver cell type.
#' @param radii Distance bins in micrometers. Default: seq(20, 100, 20).
#' @param n_perm Number of permutations. Default: 999.
#' @param seed Random seed.
#' @param adjust_method P-value adjustment method. Default: "BH".
#' @param verbose Print progress. Default: TRUE.
#'
#' @return A data.frame with columns:
#'   \item{gene}{Gene name}
#'   \item{I_ND_r1}{I_ND at first radius}
#'   \item{delta_I}{Signed delta I (spatial decay signature)}
#'   \item{p_value}{Permutation p-value for I_ND_r1}
#'   \item{p_adj}{Adjusted p-value}
#'   \item{z_score}{Standardized effect size}
#'
#' @examples
#' \dontrun{
#' result <- permutation_test_IND_sc(
#'     expr_matrix = cosmx_expr,
#'     cell_meta = meta,
#'     factor_gene = "IL1B",
#'     sender_type = "Macrophage",
#'     receiver_type = "Fibroblast",
#'     radii = seq(20, 100, 20),
#'     n_perm = 999
#' )
#'
#' # Significant responders (positive I_ND, FDR < 0.05)
#' responders <- result[result$p_adj < 0.05 & result$I_ND_r1 > 0, ]
#' }
#'
#' @export
permutation_test_IND_sc <- function(expr_matrix,
                                    cell_meta,
                                    factor_gene,
                                    sender_type,
                                    receiver_type,
                                    radii = seq(20, 100, 20),
                                    n_perm = 999,
                                    seed = NULL,
                                    adjust_method = "BH",
                                    verbose = TRUE) {

  # Convert sparse matrix if needed
  if (inherits(expr_matrix, "dgCMatrix") || inherits(expr_matrix, "dgTMatrix")) {
    expr_matrix <- as.matrix(expr_matrix)
  }
  if (!is.matrix(expr_matrix)) {
    expr_matrix <- as.matrix(expr_matrix)
  }

  # Validate metadata
  required_cols <- c("cell_id", "x", "y", "cell_type")
  missing_cols <- setdiff(required_cols, colnames(cell_meta))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # Match cell order
  if (!all(colnames(expr_matrix) == cell_meta$cell_id)) {
    cell_meta <- cell_meta[match(colnames(expr_matrix), cell_meta$cell_id), ]
  }

  coords <- as.matrix(cell_meta[, c("x", "y")])
  cell_types <- as.character(cell_meta$cell_type)

  gene_names <- rownames(expr_matrix)
  if (is.null(gene_names)) {
    gene_names <- paste0("Gene_", seq_len(nrow(expr_matrix)))
  }

  # Get factor index
  if (is.character(factor_gene)) {
    factor_idx <- which(gene_names == factor_gene)
    if (length(factor_idx) == 0) {
      stop("Factor gene '", factor_gene, "' not found")
    }
    factor_idx <- factor_idx[1] - 1L
    factor_name <- factor_gene
  } else {
    factor_idx <- as.integer(factor_gene) - 1L
    factor_name <- gene_names[factor_idx + 1]
  }

  # Get sender/receiver indices (0-based)
  sender_idx <- which(cell_types == sender_type) - 1L
  receiver_idx <- which(cell_types == receiver_type) - 1L

  if (length(sender_idx) < 10) {
    stop("Insufficient sender cells: ", length(sender_idx))
  }
  if (length(receiver_idx) < 10) {
    stop("Insufficient receiver cells: ", length(receiver_idx))
  }

  if (verbose) {
    message("Permutation test for I_ND signatures:")
    message("  Factor: ", factor_name)
    message("  Triplet: ", sender_type, " -> ", receiver_type)
    message("  Senders: ", length(sender_idx), " cells")
    message("  Receivers: ", length(receiver_idx), " cells")
    message("  Radii: ", paste(radii, collapse = ", "), " um")
    message("  Permutations: ", n_perm)
  }

  seed_int <- if (is.null(seed)) -1L else as.integer(seed)

  result <- cpp_batch_permutation_test_radii(
    expr_matrix = expr_matrix,
    gene_names = gene_names,
    factor_idx = factor_idx,
    coords = coords,
    sender_idx = as.integer(sender_idx),
    receiver_idx = as.integer(receiver_idx),
    radii = as.numeric(radii),
    n_perm = as.integer(n_perm),
    seed = seed_int,
    verbose = verbose
  )

  # Adjust p-values
  result$p_adj <- p.adjust(result$p_value, method = adjust_method)

  # Reorder columns
  result <- result[, c("gene", "I_ND_r1", "delta_I", "p_value", "p_adj", "z_score")]

  # Add metadata
  attr(result, "factor_gene") <- factor_name
  attr(result, "sender_type") <- sender_type
  attr(result, "receiver_type") <- receiver_type
  attr(result, "radii") <- radii
  attr(result, "n_perm") <- n_perm

  return(result)
}


#' Vectorized Multi-Factor Permutation Test
#'
#' Highly optimized permutation test for multiple factor genes using
#' matrix multiplication. Much faster than running separate tests.
#'
#' @param expr_matrix Gene expression matrix (genes x cells).
#' @param cell_meta Data frame with cell metadata (cell_id, x, y, cell_type).
#' @param factor_genes Character vector of factor gene names.
#' @param sender_type Sender cell type.
#' @param receiver_type Receiver cell type.
#' @param radius Distance threshold for weight matrix. Default: 50.
#' @param n_perm Number of permutations. Default: 999.
#' @param seed Random seed.
#' @param adjust_method P-value adjustment method. Default: "BH".
#' @param verbose Print progress. Default: TRUE.
#'
#' @return A named list of data.frames, one per factor gene.
#'
#' @details
#' \strong{Optimization Strategy:}
#' \enumerate{
#'   \item Precompute spatial lags for ALL genes once: LAG = W x Z_receiver
#'   \item Generate permutation indices once (shared across factors)
#'   \item For each factor, use BLAS matrix multiply: I_perm = Z_perm' x LAG
#' }
#'
#' This reduces complexity from O(n_factors x n_genes x n_perm x n_senders)
#' to O(n_genes x n_senders + n_factors x n_perm x n_genes).
#'
#' \strong{Speedup:} Typically 5-20x faster than separate calls for 5+ factors.
#'
#' @examples
#' \dontrun{
#' # Test multiple cytokines efficiently
#' results <- permutation_test_multi_factor(
#'     expr_matrix = expr,
#'     cell_meta = meta,
#'     factor_genes = c("IL1B", "TGFB1", "IFNG", "TNF", "IL6"),
#'     sender_type = "Macrophage",
#'     receiver_type = "Fibroblast",
#'     n_perm = 999
#' )
#'
#' # Combine results with FDR correction across all tests
#' all_results <- do.call(rbind, lapply(names(results), function(f) {
#'     df <- results[[f]]
#'     df$factor <- f
#'     df
#' }))
#' all_results$p_adj_global <- p.adjust(all_results$p_value, method = "BH")
#' }
#'
#' @export
permutation_test_multi_factor <- function(expr_matrix,
                                          cell_meta,
                                          factor_genes,
                                          sender_type,
                                          receiver_type,
                                          radius = 50,
                                          n_perm = 999,
                                          seed = NULL,
                                          adjust_method = "BH",
                                          verbose = TRUE) {

  # Convert sparse matrix if needed
  if (inherits(expr_matrix, "dgCMatrix") || inherits(expr_matrix, "dgTMatrix")) {
    expr_matrix <- as.matrix(expr_matrix)
  }
  if (!is.matrix(expr_matrix)) {
    expr_matrix <- as.matrix(expr_matrix)
  }

  # Validate metadata
  required_cols <- c("cell_id", "x", "y", "cell_type")
  missing_cols <- setdiff(required_cols, colnames(cell_meta))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # Match cell order
  if (!all(colnames(expr_matrix) == cell_meta$cell_id)) {
    cell_meta <- cell_meta[match(colnames(expr_matrix), cell_meta$cell_id), ]
  }

  coords <- as.matrix(cell_meta[, c("x", "y")])
  cell_types <- as.character(cell_meta$cell_type)

  gene_names <- rownames(expr_matrix)
  if (is.null(gene_names)) {
    gene_names <- paste0("Gene_", seq_len(nrow(expr_matrix)))
  }

  # Validate factor genes
  missing_factors <- setdiff(factor_genes, gene_names)
  if (length(missing_factors) > 0) {
    warning("Factor genes not found: ", paste(missing_factors, collapse = ", "))
    factor_genes <- intersect(factor_genes, gene_names)
  }
  if (length(factor_genes) == 0) {
    stop("No valid factor genes found")
  }

  # Get factor indices (0-based)
  factor_indices <- match(factor_genes, gene_names) - 1L

  # Get sender/receiver indices (0-based)
  sender_idx <- which(cell_types == sender_type) - 1L
  receiver_idx <- which(cell_types == receiver_type) - 1L

  if (length(sender_idx) < 10) {
    stop("Insufficient sender cells: ", length(sender_idx))
  }
  if (length(receiver_idx) < 10) {
    stop("Insufficient receiver cells: ", length(receiver_idx))
  }

  if (verbose) {
    message("Vectorized multi-factor permutation test:")
    message("  Factors: ", paste(factor_genes, collapse = ", "))
    message("  Triplet: ", sender_type, " -> ", receiver_type)
    message("  Senders: ", length(sender_idx), " cells")
    message("  Receivers: ", length(receiver_idx), " cells")
    message("  Genes: ", nrow(expr_matrix))
    message("  Permutations: ", n_perm)
  }

  # Build weight matrix
  if (verbose) {
    message("Building weight matrix (radius = ", radius, ")...")
  }

  W <- .build_weight_matrix_perm(
    coords = coords,
    sender_idx = sender_idx,
    receiver_idx = receiver_idx,
    r_outer = radius
  )

  seed_int <- if (is.null(seed)) -1L else as.integer(seed)

  # Run vectorized permutation test
  results <- cpp_batch_permutation_multi_factor_vectorized(
    expr_matrix = expr_matrix,
    gene_names = gene_names,
    factor_indices = as.integer(factor_indices),
    factor_names = factor_genes,
    W = W,
    sender_idx = as.integer(sender_idx),
    receiver_idx = as.integer(receiver_idx),
    n_perm = as.integer(n_perm),
    seed = seed_int,
    verbose = verbose
  )

  # Add adjusted p-values to each result
  for (f in names(results)) {
    results[[f]]$p_adj <- p.adjust(results[[f]]$p_value, method = adjust_method)
  }

  # Add metadata
  attr(results, "sender_type") <- sender_type
  attr(results, "receiver_type") <- receiver_type
  attr(results, "radius") <- radius
  attr(results, "n_perm") <- n_perm

  return(results)
}


#' Build Weight Matrix for Permutation Test
#'
#' Internal function to build sparse weight matrix.
#'
#' @param coords Coordinate matrix (N x 2)
#' @param sender_idx Sender indices (0-based)
#' @param receiver_idx Receiver indices (0-based)
#' @param r_outer Outer radius
#'
#' @return Sparse weight matrix
#'
#' @keywords internal
.build_weight_matrix_perm <- function(coords, sender_idx, receiver_idx, r_outer) {
  n_senders <- length(sender_idx)
  n_receivers <- length(receiver_idx)

  # Use C++ function if available, otherwise build in R
  if (exists("cpp_create_ring_weight_matrix_sc")) {
    return(cpp_create_ring_weight_matrix_sc(
      coords = coords,
      sender_idx = as.integer(sender_idx),
      receiver_idx = as.integer(receiver_idx),
      r_inner = 0,
      r_outer = r_outer,
      row_normalize = TRUE
    ))
  }

  # Fallback: build in R (slower)
  row_idx <- integer(0)
  col_idx <- integer(0)
  values <- numeric(0)
  row_sums <- numeric(n_senders)

  for (i in seq_len(n_senders)) {
    si <- sender_idx[i] + 1  # 1-based for R
    xi <- coords[si, 1]
    yi <- coords[si, 2]

    for (j in seq_len(n_receivers)) {
      rj <- receiver_idx[j] + 1
      if (si == rj) next

      dx <- xi - coords[rj, 1]
      dy <- yi - coords[rj, 2]
      dist <- sqrt(dx^2 + dy^2)

      if (dist < r_outer) {
        row_sums[i] <- row_sums[i] + 1
      }
    }
  }

  for (i in seq_len(n_senders)) {
    if (row_sums[i] == 0) next

    si <- sender_idx[i] + 1
    xi <- coords[si, 1]
    yi <- coords[si, 2]

    for (j in seq_len(n_receivers)) {
      rj <- receiver_idx[j] + 1
      if (si == rj) next

      dx <- xi - coords[rj, 1]
      dy <- yi - coords[rj, 2]
      dist <- sqrt(dx^2 + dy^2)

      if (dist < r_outer) {
        row_idx <- c(row_idx, i)
        col_idx <- c(col_idx, j)
        values <- c(values, 1 / row_sums[i])
      }
    }
  }

  Matrix::sparseMatrix(
    i = row_idx,
    j = col_idx,
    x = values,
    dims = c(n_senders, n_receivers)
  )
}
