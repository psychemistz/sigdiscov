#' @title Single-Cell Spatial Transcriptomics Analysis
#' @description Functions for computing cell-type-aware I_ND (Normalized Directional
#' Moran's I) signatures for single-cell spatial transcriptomics data (CosMx, Xenium,
#' MERFISH, etc.).
#' @name singlecell
NULL


#' Compute I_ND Signatures for Single-Cell Spatial Transcriptomics
#'
#' Computes cell-type-aware I_ND (Normalized Directional Moran's I) signatures
#' for spatial transcriptomics data with single-cell resolution.
#'
#' @param expr_matrix Gene expression matrix (genes x cells). Sparse matrix (dgCMatrix) supported.
#' @param cell_meta Data frame with cell metadata. Must contain columns:
#'   \itemize{
#'     \item cell_id: unique identifier matching column names of expr_matrix
#'     \item x: x-coordinate in micrometers
#'     \item y: y-coordinate in micrometers
#'     \item cell_type: cell type annotation
#'   }
#' @param factor_gene Name or row index of the factor gene.
#' @param sender_type Cell type for sender cells (e.g., "Macrophage").
#' @param receiver_type Cell type for receiver cells (e.g., "Fibroblast").
#' @param radii Distance bins in micrometers. Default: seq(20, 200, 20).
#' @param min_cells Minimum cells per type required. Default: 30.
#' @param smooth_window Savitzky-Golay smoothing window size. Default: 5.
#' @param smooth_poly Savitzky-Golay polynomial order. Default: 2.
#' @param verbose Print progress messages. Default: TRUE.
#'
#' @return A data frame with columns:
#'   \item{gene}{Gene name}
#'   \item{I_ND_r1}{I_ND at first (nearest) distance}
#'   \item{I_ND_max}{Maximum I_ND across all distances}
#'   \item{delta_I}{Signed delta I (max - min, signed by trend)}
#'   \item{I_ND_20, I_ND_40, ...}{I_ND at each specified radius}
#'
#' @details
#' This function implements the (factor, sender_type, receiver_type) triplet structure
#' for analyzing cell-cell communication in single-cell spatial data:
#'
#' \strong{Mathematical Formulation:}
#' \deqn{I_{ND}(f,g,r)_{U \rightarrow V} = \frac{z_U^f \cdot \tilde{W}_r \cdot z_V^g}{\|z_U^f\| \times \|\tilde{W}_r \cdot z_V^g\|}}
#'
#' Where:
#' \itemize{
#'   \item \eqn{z_U^f} = standardized factor expression in sender cells (type U)
#'   \item \eqn{z_V^g} = standardized gene expression in receiver cells (type V)
#'   \item \eqn{\tilde{W}_r} = row-normalized weight matrix at distance r
#' }
#'
#' \strong{Interpretation:}
#' \itemize{
#'   \item delta_I > 0: Decay pattern (spatial responder)
#'   \item delta_I < 0: Increase pattern (avoidance)
#'   \item delta_I ~ 0: No distance dependence (constitutive)
#' }
#'
#' @examples
#' \dontrun{
#' # Load CosMx data
#' expr <- readRDS("cosmx_expression.rds")
#' meta <- read.csv("cosmx_metadata.csv")
#'
#' # Compute I_ND for IL1B: Macrophage -> Fibroblast
#' signatures <- compute_IND_sc(
#'     expr_matrix = expr,
#'     cell_meta = meta,
#'     factor_gene = "IL1B",
#'     sender_type = "Macrophage",
#'     receiver_type = "Fibroblast",
#'     radii = seq(20, 200, 20)
#' )
#'
#' # Top responder genes
#' head(signatures[order(-signatures$I_ND_r1), ], 20)
#' }
#'
#' @export
compute_IND_sc <- function(expr_matrix,
                           cell_meta,
                           factor_gene,
                           sender_type,
                           receiver_type,
                           radii = seq(20, 200, 20),
                           min_cells = 30,
                           smooth_window = 5,
                           smooth_poly = 2,
                           verbose = TRUE) {

  # Convert sparse matrix if needed
  if (inherits(expr_matrix, "dgCMatrix") || inherits(expr_matrix, "dgTMatrix")) {
    expr_matrix <- as.matrix(expr_matrix)
  }
  if (!is.matrix(expr_matrix)) {
    expr_matrix <- as.matrix(expr_matrix)
  }

  # Check required columns in metadata
  required_cols <- c("cell_id", "x", "y", "cell_type")
  missing_cols <- setdiff(required_cols, colnames(cell_meta))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in cell_meta: ", paste(missing_cols, collapse = ", "))
  }

  # Match cell order between expression matrix and metadata
  if (!all(colnames(expr_matrix) == cell_meta$cell_id)) {
    if (!all(colnames(expr_matrix) %in% cell_meta$cell_id)) {
      stop("Column names of expr_matrix must match cell_id in cell_meta")
    }
    cell_meta <- cell_meta[match(colnames(expr_matrix), cell_meta$cell_id), ]
  }

  # Extract coordinates and cell types
  coords <- as.matrix(cell_meta[, c("x", "y")])
  cell_types <- as.character(cell_meta$cell_type)

  # Validate cell types
  available_types <- unique(cell_types)
  if (!sender_type %in% available_types) {
    stop("Sender type '", sender_type, "' not found. Available types: ",
         paste(available_types, collapse = ", "))
  }
  if (!receiver_type %in% available_types) {
    stop("Receiver type '", receiver_type, "' not found. Available types: ",
         paste(available_types, collapse = ", "))
  }

  # Get sender and receiver indices (0-based for C++)
  sender_idx <- which(cell_types == sender_type) - 1L
  receiver_idx <- which(cell_types == receiver_type) - 1L

  n_senders <- length(sender_idx)
  n_receivers <- length(receiver_idx)

  if (n_senders < min_cells) {
    stop("Insufficient sender cells: ", n_senders, " < ", min_cells)
  }
  if (n_receivers < min_cells) {
    stop("Insufficient receiver cells: ", n_receivers, " < ", min_cells)
  }

  # Get factor gene index
  gene_names <- rownames(expr_matrix)
  if (is.null(gene_names)) {
    gene_names <- paste0("Gene_", seq_len(nrow(expr_matrix)))
  }

  if (is.character(factor_gene)) {
    factor_idx <- which(gene_names == factor_gene)
    if (length(factor_idx) == 0) {
      stop("Factor gene '", factor_gene, "' not found in expression matrix")
    }
    factor_idx <- factor_idx[1] - 1L
    factor_name <- factor_gene
  } else {
    factor_idx <- as.integer(factor_gene) - 1L
    factor_name <- gene_names[factor_idx + 1]
  }

  # Ensure smooth_window is odd
  if (smooth_window %% 2 == 0) {
    smooth_window <- smooth_window + 1L
  }

  if (verbose) {
    message("Computing I_ND signatures for single-cell ST:")
    message("  Factor: ", factor_name)
    message("  Triplet: ", sender_type, " -> ", receiver_type)
    message("  Senders: ", n_senders, " cells")
    message("  Receivers: ", n_receivers, " cells")
    message("  Radii: ", paste(radii, collapse = ", "), " um")
    message("  ", nrow(expr_matrix), " genes x ", ncol(expr_matrix), " cells")
  }

  # Call C++ implementation
  result <- cpp_compute_IND_signatures_sc(
    expr_matrix = expr_matrix,
    gene_names = gene_names,
    factor_idx = factor_idx,
    coords = coords,
    sender_idx = as.integer(sender_idx),
    receiver_idx = as.integer(receiver_idx),
    radii = as.numeric(radii),
    smooth_window = as.integer(smooth_window),
    smooth_poly = as.integer(smooth_poly),
    verbose = verbose
  )

  # Add metadata attributes
  attr(result, "factor_gene") <- factor_name
  attr(result, "sender_type") <- sender_type
  attr(result, "receiver_type") <- receiver_type
  attr(result, "n_senders") <- n_senders
  attr(result, "n_receivers") <- n_receivers
  attr(result, "radii") <- radii

  return(result)
}


#' Get I_ND Curve for a Single Gene Pair (Single-Cell)
#'
#' Computes the full I_ND(r) curve for a specific (factor, target, sender_type, receiver_type)
#' combination.
#'
#' @param expr_matrix Gene expression matrix (genes x cells).
#' @param cell_meta Data frame with cell metadata.
#' @param factor_gene Name or row index of the factor gene.
#' @param target_gene Name or row index of the target gene.
#' @param sender_type Cell type for sender cells.
#' @param receiver_type Cell type for receiver cells.
#' @param radii Distance bins in micrometers. Default: seq(20, 200, 20).
#' @param smooth_window Smoothing window size. Default: 5.
#' @param smooth_poly Polynomial order. Default: 2.
#'
#' @return A list with:
#'   \item{radii}{Distance values}
#'   \item{I_raw}{Raw I_ND values at each distance}
#'   \item{I_smooth}{Smoothed I_ND curve}
#'   \item{delta_I}{Signed delta I}
#'   \item{I_max, I_min}{Curve extrema}
#'   \item{factor_gene, target_gene}{Gene names}
#'   \item{sender_type, receiver_type}{Cell types}
#'
#' @export
get_IND_curve_sc <- function(expr_matrix,
                             cell_meta,
                             factor_gene,
                             target_gene,
                             sender_type,
                             receiver_type,
                             radii = seq(20, 200, 20),
                             smooth_window = 5,
                             smooth_poly = 2) {

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

  # Get indices
  factor_idx <- .get_gene_index(factor_gene, gene_names)
  target_idx <- .get_gene_index(target_gene, gene_names)

  factor_name <- gene_names[factor_idx]
  target_name <- gene_names[target_idx]

  # Get cell type indices (0-based)
  sender_idx <- which(cell_types == sender_type) - 1L
  receiver_idx <- which(cell_types == receiver_type) - 1L

  if (length(sender_idx) == 0) stop("No cells of sender type: ", sender_type)
  if (length(receiver_idx) == 0) stop("No cells of receiver type: ", receiver_type)

  # Extract expression vectors
  factor_expr <- as.numeric(expr_matrix[factor_idx, sender_idx + 1])
  target_expr <- as.numeric(expr_matrix[target_idx, receiver_idx + 1])

  if (smooth_window %% 2 == 0) {
    smooth_window <- smooth_window + 1L
  }

  # Call C++ to compute curve
  result <- cpp_compute_IND_curve_sc(
    factor_expr = factor_expr,
    gene_expr = target_expr,
    coords = coords,
    sender_idx = as.integer(sender_idx),
    receiver_idx = as.integer(receiver_idx),
    radii = as.numeric(radii),
    smooth_window = as.integer(smooth_window),
    smooth_poly = as.integer(smooth_poly)
  )

  result$factor_gene <- factor_name
  result$target_gene <- target_name
  result$sender_type <- sender_type
  result$receiver_type <- receiver_type

  return(result)
}


#' Plot I_ND Curve for Single-Cell ST
#'
#' Visualize the distance-dependent I_ND curve for a cell-type-specific triplet.
#'
#' @param curve_result Result from \code{\link{get_IND_curve_sc}}.
#' @param show_raw Show raw curve points. Default: TRUE.
#' @param show_smooth Show smoothed curve line. Default: TRUE.
#' @param title Custom plot title. If NULL, auto-generated.
#' @param ... Additional arguments passed to plot().
#'
#' @return Invisibly returns the curve_result.
#'
#' @export
plot_IND_curve_sc <- function(curve_result,
                              show_raw = TRUE,
                              show_smooth = TRUE,
                              title = NULL,
                              ...) {

  radii <- curve_result$radii
  I_raw <- curve_result$I_raw
  I_smooth <- curve_result$I_smooth

  if (is.null(title)) {
    title <- sprintf("%s: %s -> %s\nTarget: %s | delta_I = %.4f",
                     curve_result$factor_gene,
                     curve_result$sender_type,
                     curve_result$receiver_type,
                     curve_result$target_gene,
                     curve_result$delta_I_signed)
  }

  all_vals <- c(I_raw, I_smooth)
  all_vals <- all_vals[is.finite(all_vals)]
  ylim <- if (length(all_vals) == 0) c(-1, 1) else range(all_vals) * 1.1
  if (ylim[1] > 0) ylim[1] <- 0
  if (ylim[2] < 0) ylim[2] <- 0

  plot(radii, I_raw, type = "n",
       xlab = "Distance (um)", ylab = "I_ND",
       main = title, ylim = ylim, ...)

  abline(h = 0, lty = 2, col = "gray50")

  if (show_raw) {
    points(radii, I_raw, pch = 19, col = "gray40", cex = 1.2)
  }

  if (show_smooth) {
    valid_idx <- which(is.finite(I_smooth))
    if (length(valid_idx) > 1) {
      lines(radii[valid_idx], I_smooth[valid_idx], col = "blue", lwd = 2)
    }
  }

  invisible(curve_result)
}


#' Compute I_ND for Multiple Factors (Optimized)
#'
#' Efficiently computes I_ND signatures for multiple factor genes against a single
#' (sender_type, receiver_type) pair. This is much faster than calling compute_IND_sc
#' multiple times because weight matrices are computed only once.
#'
#' @param expr_matrix Gene expression matrix (genes x cells).
#' @param cell_meta Data frame with cell metadata.
#' @param factor_genes Character vector of factor gene names.
#' @param sender_type Cell type for sender cells.
#' @param receiver_type Cell type for receiver cells.
#' @param radii Distance bins in micrometers. Default: seq(20, 200, 20).
#' @param min_cells Minimum cells per type required. Default: 30.
#' @param smooth_window Smoothing window size. Default: 5.
#' @param smooth_poly Polynomial order. Default: 2.
#' @param verbose Print progress messages. Default: TRUE.
#'
#' @return A named list of data frames, one per factor gene.
#'
#' @details
#' \strong{Optimization:} This function precomputes:
#' \enumerate{
#'   \item Weight matrices for all radii (once per sender-receiver pair)
#'   \item Spatial lags for all genes (W * z_g for each radius)
#' }
#' Then it only needs to compute factor-specific operations for each factor,
#' which is O(n_factors * n_genes * n_radii) dot products instead of
#' O(n_factors * n_senders * n_receivers) distance calculations.
#'
#' @examples
#' \dontrun{
#' # Compute signatures for multiple cytokines efficiently
#' results <- compute_IND_multi_factor(
#'     expr_matrix = expr,
#'     cell_meta = meta,
#'     factor_genes = c("IL1B", "TGFB1", "IFNG", "TNF"),
#'     sender_type = "Macrophage",
#'     receiver_type = "Fibroblast"
#' )
#'
#' # Access results for each factor
#' il1b_sig <- results$IL1B
#' }
#'
#' @export
compute_IND_multi_factor <- function(expr_matrix,
                                     cell_meta,
                                     factor_genes,
                                     sender_type,
                                     receiver_type,
                                     radii = seq(20, 200, 20),
                                     min_cells = 30,
                                     smooth_window = 5,
                                     smooth_poly = 2,
                                     verbose = TRUE) {

  # Convert sparse matrix if needed
  if (inherits(expr_matrix, "dgCMatrix") || inherits(expr_matrix, "dgTMatrix")) {
    expr_matrix <- as.matrix(expr_matrix)
  }
  if (!is.matrix(expr_matrix)) {
    expr_matrix <- as.matrix(expr_matrix)
  }

  # Check metadata
  required_cols <- c("cell_id", "x", "y", "cell_type")
  missing_cols <- setdiff(required_cols, colnames(cell_meta))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # Match cell order
  if (!all(colnames(expr_matrix) == cell_meta$cell_id)) {
    if (!all(colnames(expr_matrix) %in% cell_meta$cell_id)) {
      stop("Column names of expr_matrix must match cell_id in cell_meta")
    }
    cell_meta <- cell_meta[match(colnames(expr_matrix), cell_meta$cell_id), ]
  }

  coords <- as.matrix(cell_meta[, c("x", "y")])
  cell_types <- as.character(cell_meta$cell_type)

  # Get sender/receiver indices
  sender_idx <- which(cell_types == sender_type) - 1L
  receiver_idx <- which(cell_types == receiver_type) - 1L

  if (length(sender_idx) < min_cells) {
    stop("Insufficient sender cells: ", length(sender_idx), " < ", min_cells)
  }
  if (length(receiver_idx) < min_cells) {
    stop("Insufficient receiver cells: ", length(receiver_idx), " < ", min_cells)
  }

  # Get gene names and validate factors
  gene_names <- rownames(expr_matrix)
  if (is.null(gene_names)) {
    gene_names <- paste0("Gene_", seq_len(nrow(expr_matrix)))
  }

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

  if (smooth_window %% 2 == 0) {
    smooth_window <- smooth_window + 1L
  }

  if (verbose) {
    message("Computing I_ND for ", length(factor_genes), " factors:")
    message("  Factors: ", paste(factor_genes, collapse = ", "))
    message("  Triplet: ", sender_type, " -> ", receiver_type)
    message("  Senders: ", length(sender_idx), " | Receivers: ", length(receiver_idx))
    message("  ", nrow(expr_matrix), " genes x ", ncol(expr_matrix), " cells")
  }

  result <- cpp_compute_IND_multi_factor(
    expr_matrix = expr_matrix,
    gene_names = gene_names,
    factor_indices = as.integer(factor_indices),
    factor_names = factor_genes,
    coords = coords,
    sender_idx = as.integer(sender_idx),
    receiver_idx = as.integer(receiver_idx),
    radii = as.numeric(radii),
    smooth_window = as.integer(smooth_window),
    smooth_poly = as.integer(smooth_poly),
    verbose = verbose
  )

  # Add metadata
  attr(result, "sender_type") <- sender_type
  attr(result, "receiver_type") <- receiver_type
  attr(result, "radii") <- radii

  return(result)
}


#' Compute I_ND for All Cell Type Pairs and Multiple Factors (Optimized)
#'
#' Efficiently computes I_ND signatures for multiple factor genes across all
#' sender-receiver cell type combinations. Optimized to precompute weight matrices
#' once per cell type pair.
#'
#' @param expr_matrix Gene expression matrix (genes x cells).
#' @param cell_meta Data frame with cell metadata.
#' @param factor_genes Character vector of factor gene names. If NULL, uses a single factor.
#' @param cell_types_include Cell types to include. Default: all types with sufficient cells.
#' @param radii Distance bins in micrometers. Default: seq(20, 200, 20).
#' @param min_cells Minimum cells per type required. Default: 30.
#' @param verbose Print progress messages. Default: TRUE.
#'
#' @return A nested named list: result[[pair_name]][[factor_name]] gives the data.frame
#'   for that combination. If only one factor, result[[pair_name]] gives the data.frame directly.
#'
#' @details
#' \strong{Optimization Strategy:}
#' For each (sender_type, receiver_type) pair:
#' \enumerate{
#'   \item Compute weight matrices once (expensive: O(n_senders * n_receivers))
#'   \item Compute spatial lags for all genes once
#'   \item Process all factors using precomputed data (cheap: just dot products)
#' }
#'
#' This is much faster than the naive approach of calling compute_IND_sc for each
#' (factor, sender, receiver) combination separately.
#'
#' @examples
#' \dontrun{
#' # Compute for multiple cytokines across all cell type pairs
#' all_results <- compute_IND_all_pairs(
#'     expr_matrix = expr,
#'     cell_meta = meta,
#'     factor_genes = c("IL1B", "TGFB1", "IFNG"),
#'     cell_types_include = c("Macrophage", "Fibroblast", "Epithelial")
#' )
#'
#' # Access results
#' il1b_macro_fibro <- all_results$Macrophage_to_Fibroblast$IL1B
#' }
#'
#' @export
compute_IND_all_pairs <- function(expr_matrix,
                                  cell_meta,
                                  factor_genes,
                                  cell_types_include = NULL,
                                  radii = seq(20, 200, 20),
                                  min_cells = 30,
                                  verbose = TRUE) {

  # Handle single factor case (backward compatibility)
  if (length(factor_genes) == 1) {
    single_factor <- TRUE
  } else {
    single_factor <- FALSE
  }

  # Get available cell types
  all_types <- unique(as.character(cell_meta$cell_type))

  # Filter by type counts
  type_counts <- table(cell_meta$cell_type)
  valid_types <- names(type_counts[type_counts >= min_cells])

  if (!is.null(cell_types_include)) {
    valid_types <- intersect(valid_types, cell_types_include)
  }

  if (length(valid_types) == 0) {
    stop("No cell types have >= ", min_cells, " cells")
  }

  if (verbose) {
    message("Computing I_ND for ", length(factor_genes), " factor(s) x ",
            length(valid_types), " cell types")
    message("  Factors: ", paste(factor_genes, collapse = ", "))
    message("  Cell types: ", paste(valid_types, collapse = ", "))
    message("  Total pairs: ", length(valid_types)^2)
  }

  # Generate all pairs
  pairs <- expand.grid(sender = valid_types, receiver = valid_types,
                       stringsAsFactors = FALSE)

  results <- list()

  for (i in seq_len(nrow(pairs))) {
    sender <- pairs$sender[i]
    receiver <- pairs$receiver[i]
    pair_name <- paste0(sender, "_to_", receiver)

    if (verbose) {
      message("\nProcessing pair ", i, "/", nrow(pairs), ": ", pair_name)
    }

    tryCatch({
      pair_result <- compute_IND_multi_factor(
        expr_matrix = expr_matrix,
        cell_meta = cell_meta,
        factor_genes = factor_genes,
        sender_type = sender,
        receiver_type = receiver,
        radii = radii,
        min_cells = min_cells,
        verbose = FALSE
      )

      # For single factor, return the data.frame directly (backward compatible)
      if (single_factor) {
        results[[pair_name]] <- pair_result[[1]]
      } else {
        results[[pair_name]] <- pair_result
      }

    }, error = function(e) {
      if (verbose) message("  Skipped: ", e$message)
    })
  }

  attr(results, "factor_genes") <- factor_genes
  attr(results, "cell_types") <- valid_types
  attr(results, "radii") <- radii

  return(results)
}


#' Summarize Single-Cell ST Cell Type Distribution
#'
#' Provides a summary of cell types and their counts in the dataset.
#'
#' @param cell_meta Data frame with cell metadata.
#' @param min_cells Minimum cells to highlight. Default: 30.
#'
#' @return A data frame with cell type counts and whether they meet the threshold.
#'
#' @export
summarize_cell_types <- function(cell_meta, min_cells = 30) {
  type_counts <- as.data.frame(table(cell_meta$cell_type))
  colnames(type_counts) <- c("cell_type", "count")
  type_counts <- type_counts[order(-type_counts$count), ]
  type_counts$sufficient <- type_counts$count >= min_cells
  rownames(type_counts) <- NULL
  return(type_counts)
}
