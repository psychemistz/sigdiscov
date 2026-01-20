#' @title Signed Delta I Signatures for Spatial Transcriptomics
#' @description Functions for computing distance-dependent spatial correlation
#' signatures (signed delta I) between genes in spatial transcriptomics data.
#' @name signed_delta_I
NULL


#' Compute Signed Delta I for a Single Factor
#'
#' Computes distance-dependent spatial correlation between a factor gene and all
#' other genes, then derives signed delta I signatures.
#'
#' @param expr_matrix A matrix of gene expression values (genes x spots). VST-normalized recommended.
#' @param spot_coords Spot coordinates as a data frame with row/col columns, or a 2-column matrix.
#' @param factor_gene Name or row index of the factor gene.
#' @param radii Numeric vector of distance bin edges. Default: seq(100, 600, 100).
#' @param mode Computation mode: "bivariate" (Moran's I) or "directional" (I_ND). Default: "bivariate".
#' @param sender_percentile For directional mode: percentile threshold for senders. Default: 75.
#' @param smooth_window Savitzky-Golay smoothing window size (must be odd). Default: 5.
#' @param smooth_poly Savitzky-Golay polynomial order. Default: 2.
#' @param coord_scale Scale factor for coordinates. Default: 1.
#' @param verbose Print progress messages. Default: TRUE.
#'
#' @return A data frame with columns: gene, delta_I, delta_I_signed, sign, I_short, I_long, I_max, I_min, auc, peak_radius_idx.
#'
#' @details
#' Interpretation:
#' \itemize{
#'   \item delta_I_signed > 0: Decay pattern (TRUE RESPONDER)
#'   \item delta_I_signed < 0: Increase pattern (AVOIDANCE)
#'   \item delta_I_signed ~ 0: Flat pattern (CONSTITUTIVE)
#' }
#'
#' @examples
#' \dontrun{
#' sig <- compute_signed_delta_I(expr_vst, coords, factor_gene = "IL1B")
#' head(sig[order(-sig$delta_I_signed), ])
#' }
#'
#' @export
compute_signed_delta_I <- function(expr_matrix,
                                    spot_coords,
                                    factor_gene,
                                    radii = seq(100, 600, 100),
                                    mode = c("bivariate", "directional"),
                                    sender_percentile = 75,
                                    smooth_window = 5,
                                    smooth_poly = 2,
                                    coord_scale = 1,
                                    verbose = TRUE) {

  mode <- match.arg(mode)
  mode_int <- if (mode == "bivariate") 0L else 1L

  if (!is.matrix(expr_matrix)) {
    expr_matrix <- as.matrix(expr_matrix)
  }

  # Get factor index
  if (is.character(factor_gene)) {
    factor_idx <- which(rownames(expr_matrix) == factor_gene)
    if (length(factor_idx) == 0) {
      stop("Factor gene '", factor_gene, "' not found in expression matrix")
    }
    factor_idx <- factor_idx[1] - 1L
  } else {
    factor_idx <- as.integer(factor_gene) - 1L
  }

  # Prepare coordinates
  coords <- .prepare_coords(spot_coords, coord_scale)

  if (ncol(expr_matrix) != nrow(coords)) {
    stop("Number of spots in expr_matrix and spot_coords must match")
  }

  gene_names <- rownames(expr_matrix)
  if (is.null(gene_names)) {
    gene_names <- paste0("Gene_", seq_len(nrow(expr_matrix)))
  }

  if (smooth_window %% 2 == 0) {
    smooth_window <- smooth_window + 1L
  }

  if (verbose) {
    message("Computing signed delta I signatures...")
    message("  Mode: ", mode)
    message("  Factor: ", if (is.character(factor_gene)) factor_gene else gene_names[factor_idx + 1])
    message("  ", nrow(expr_matrix), " genes x ", ncol(expr_matrix), " spots")
  }

  result <- cpp_compute_all_signatures(
    expr_matrix = expr_matrix,
    gene_names = gene_names,
    factor_idx = factor_idx,
    coords = coords,
    radii = as.numeric(radii),
    mode = mode_int,
    sender_percentile = as.numeric(sender_percentile),
    smooth_window = as.integer(smooth_window),
    smooth_poly = as.integer(smooth_poly),
    verbose = verbose
  )

  attr(result, "mode") <- mode
  attr(result, "factor_gene") <- if (is.character(factor_gene)) factor_gene else gene_names[factor_idx + 1]
  attr(result, "radii") <- radii

  return(result)
}


#' Get Moran's I Curve for a Gene Pair
#'
#' Computes the full I(r) curve for a single factor-gene pair.
#'
#' @param expr_matrix A matrix of gene expression values (genes x spots).
#' @param spot_coords Spot coordinates as a data frame or matrix.
#' @param factor_gene Name or row index of the factor gene.
#' @param target_gene Name or row index of the target gene.
#' @param radii Numeric vector of distance bin edges. Default: seq(100, 600, 100).
#' @param mode Computation mode: "bivariate" or "directional". Default: "bivariate".
#' @param sender_percentile For directional mode: percentile threshold. Default: 75.
#' @param smooth_window Savitzky-Golay window size. Default: 5.
#' @param smooth_poly Polynomial order. Default: 2.
#' @param coord_scale Scale factor for coordinates. Default: 1.
#'
#' @return A list with radii, I_raw, I_smooth, delta_I_signed, and other metrics.
#'
#' @examples
#' \dontrun{
#' curve <- get_moran_curve(expr_vst, coords, "IL1B", "COL1A1")
#' plot(curve$radii, curve$I_smooth, type = "l")
#' }
#'
#' @export
get_moran_curve <- function(expr_matrix,
                            spot_coords,
                            factor_gene,
                            target_gene,
                            radii = seq(100, 600, 100),
                            mode = c("bivariate", "directional"),
                            sender_percentile = 75,
                            smooth_window = 5,
                            smooth_poly = 2,
                            coord_scale = 1) {

  mode <- match.arg(mode)
  mode_int <- if (mode == "bivariate") 0L else 1L

  if (!is.matrix(expr_matrix)) {
    expr_matrix <- as.matrix(expr_matrix)
  }

  gene_names <- rownames(expr_matrix)
  if (is.null(gene_names)) {
    gene_names <- paste0("Gene_", seq_len(nrow(expr_matrix)))
  }

  # Get indices
  factor_idx <- .get_gene_index(factor_gene, gene_names)
  target_idx <- .get_gene_index(target_gene, gene_names)

  factor_name <- gene_names[factor_idx]
  target_name <- gene_names[target_idx]

  coords <- .prepare_coords(spot_coords, coord_scale)

  if (smooth_window %% 2 == 0) {
    smooth_window <- smooth_window + 1L
  }

  factor_expr <- as.numeric(expr_matrix[factor_idx, ])
  gene_expr <- as.numeric(expr_matrix[target_idx, ])

  result <- cpp_compute_single_signature(
    factor_expr = factor_expr,
    gene_expr = gene_expr,
    coords = coords,
    radii = as.numeric(radii),
    mode = mode_int,
    sender_percentile = as.numeric(sender_percentile),
    smooth_window = as.integer(smooth_window),
    smooth_poly = as.integer(smooth_poly)
  )

  result$factor_gene <- factor_name
  result$target_gene <- target_name
  result$mode <- mode

  return(result)
}


#' Plot Moran's I Curve
#'
#' Visualize the distance-dependent correlation curve for a gene pair.
#'
#' @param curve_result Result from \code{\link{get_moran_curve}}.
#' @param show_raw Show raw curve points. Default: TRUE.
#' @param show_smooth Show smoothed curve line. Default: TRUE.
#' @param title Custom plot title. If NULL, auto-generated.
#' @param ... Additional arguments passed to plot().
#'
#' @return Invisibly returns the curve_result.
#'
#' @export
plot_moran_curve <- function(curve_result,
                              show_raw = TRUE,
                              show_smooth = TRUE,
                              title = NULL,
                              ...) {

  radii <- curve_result$radii
  I_raw <- curve_result$I_raw
  I_smooth <- curve_result$I_smooth

  if (is.null(title)) {
    sign_text <- if (curve_result$sign > 0) "decay" else "increase"
    title <- sprintf("%s -> %s\ndelta_I_signed = %.4f (%s)",
                     curve_result$factor_gene,
                     curve_result$target_gene,
                     curve_result$delta_I_signed,
                     sign_text)
  }

  all_vals <- c(I_raw, I_smooth)
  all_vals <- all_vals[is.finite(all_vals)]
  ylim <- if (length(all_vals) == 0) c(-1, 1) else range(all_vals) * 1.1
  if (ylim[1] > 0) ylim[1] <- 0
  if (ylim[2] < 0) ylim[2] <- 0

  plot(radii, I_raw, type = "n",
       xlab = "Distance", ylab = "Spatial Correlation (I)",
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


#' Compute Delta I Matrix
#'
#' Compute delta I matrix for ALL gene pairs. Supports 4 methods:
#' \itemize{
#'   \item ring + moran: Ring-based Moran's I (default)
#'   \item ring + ind: Ring-based I_ND
#'   \item circular + moran: Circular Moran's I
#'   \item circular + ind: Circular I_ND
#' }
#'
#' @param expr_matrix A matrix of gene expression values (genes x spots). VST-normalized recommended.
#' @param spot_coords Spot coordinates as a data frame with row/col columns, or a 2-column matrix.
#' @param radii Numeric vector of distance bin edges. Default: seq(150, 650, 100).
#' @param weight_type Neighbor definition: "ring" or "circular". Default: "ring".
#' @param correlation_type Correlation measure: "moran" or "ind". Default: "moran".
#' @param chunk_size Number of genes per chunk for memory efficiency. Default: 1000.
#' @param smooth_window Savitzky-Golay smoothing window size (must be odd). Default: 5.
#' @param smooth_poly Savitzky-Golay polynomial order. Default: 2.
#' @param coord_scale Scale factor for coordinates. Default: 1.
#' @param verbose Print progress messages. Default: TRUE.
#'
#' @return A list with:
#' \itemize{
#'   \item delta_I_signed: n_genes x n_genes matrix
#'   \item n_genes, n_spots, n_radii: dimensions
#'   \item weight_type, correlation_type: method parameters
#' }
#'
#' @details
#' \strong{Weight Types:}
#' \itemize{
#'   \item Ring: Neighbors in distance band [r_inner, r_outer)
#'   \item Circular: All neighbors in cumulative disk [0, r_outer)
#' }
#'
#' \strong{Correlation Types:}
#' \itemize{
#'   \item Moran's I: I = z_f' * W * z_g / n
#'   \item I_ND: I_ND = dot(z_f, W*z_g) / (||z_f|| * ||W*z_g||)
#' }
#'
#' @examples
#' \dontrun{
#' # Ring-based Moran's I (default)
#' result <- compute_delta_I_matrix(expr_vst, coords)
#'
#' # Circular I_ND
#' result <- compute_delta_I_matrix(expr_vst, coords,
#'                                   weight_type = "circular",
#'                                   correlation_type = "ind")
#' }
#'
#' @export
compute_delta_I_matrix <- function(expr_matrix,
                                    spot_coords,
                                    radii = seq(150, 650, 100),
                                    weight_type = c("ring", "circular"),
                                    correlation_type = c("moran", "ind"),
                                    chunk_size = 1000,
                                    smooth_window = 5,
                                    smooth_poly = 2,
                                    coord_scale = 1,
                                    verbose = TRUE) {

  weight_type <- match.arg(weight_type)
  correlation_type <- match.arg(correlation_type)

  wt_int <- if (weight_type == "ring") 0L else 1L
  ct_int <- if (correlation_type == "moran") 0L else 1L

  if (!is.matrix(expr_matrix)) {
    expr_matrix <- as.matrix(expr_matrix)
  }

  coords <- .prepare_coords(spot_coords, coord_scale)

  if (ncol(expr_matrix) != nrow(coords)) {
    stop("Number of spots in expr_matrix and spot_coords must match")
  }

  if (smooth_window %% 2 == 0) {
    smooth_window <- smooth_window + 1L
    if (verbose) message("Adjusted smooth_window to ", smooth_window)
  }

  if (verbose) {
    message("Computing delta I matrix...")
    message("  Weight type: ", weight_type)
    message("  Correlation type: ", correlation_type)
    message("  ", nrow(expr_matrix), " genes x ", ncol(expr_matrix), " spots")
    message("  ", length(radii), " distance bins")
  }

  result <- cpp_compute_delta_I_matrix_unified(
    expr_matrix = expr_matrix,
    coords = coords,
    radii = as.numeric(radii),
    weight_type = wt_int,
    correlation_type = ct_int,
    chunk_size = as.integer(chunk_size),
    smooth_window = as.integer(smooth_window),
    smooth_poly = as.integer(smooth_poly),
    verbose = verbose
  )

  gene_names <- rownames(expr_matrix)
  if (!is.null(gene_names)) {
    rownames(result$delta_I_signed) <- gene_names
    colnames(result$delta_I_signed) <- gene_names
  }

  return(result)
}


#' Compute All 4 Delta Types for Factor Genes
#'
#' Computes all 4 delta types (ring/circular x moran/ind) for specified factor genes.
#'
#' @param expr_matrix A matrix of gene expression values (genes x spots).
#' @param spot_coords Spot coordinates as a data frame or matrix.
#' @param factor_genes Character vector of factor gene names.
#' @param radii Numeric vector of distance bin edges. Default: seq(150, 650, 100).
#' @param smooth_window Savitzky-Golay window size. Default: 5.
#' @param smooth_poly Polynomial order. Default: 2.
#' @param coord_scale Scale factor for coordinates. Default: 1.
#' @param verbose Print progress messages. Default: TRUE.
#'
#' @return If single factor: DataFrame with columns gene, delta_i_ring_moran, delta_i_cir_moran,
#'   delta_i_ring_ind, delta_i_cir_ind. If multiple factors: named list of DataFrames.
#'
#' @examples
#' \dontrun{
#' results <- compute_four_deltas(expr_vst, coords, c("IFNG", "TGFB1"))
#' cor(results$IFNG[, 2:5])
#' }
#'
#' @export
compute_four_deltas <- function(expr_matrix,
                                 spot_coords,
                                 factor_genes,
                                 radii = seq(150, 650, 100),
                                 smooth_window = 5,
                                 smooth_poly = 2,
                                 coord_scale = 1,
                                 verbose = TRUE) {

  if (!is.matrix(expr_matrix)) {
    expr_matrix <- as.matrix(expr_matrix)
  }

  coords <- .prepare_coords(spot_coords, coord_scale)

  if (ncol(expr_matrix) != nrow(coords)) {
    stop("Number of spots in expr_matrix and spot_coords must match")
  }

  gene_names <- rownames(expr_matrix)
  if (is.null(gene_names)) {
    gene_names <- paste0("Gene_", seq_len(nrow(expr_matrix)))
  }

  missing_genes <- setdiff(factor_genes, gene_names)
  if (length(missing_genes) > 0) {
    warning("Factor genes not found: ", paste(missing_genes, collapse = ", "))
    factor_genes <- intersect(factor_genes, gene_names)
  }

  if (length(factor_genes) == 0) {
    stop("No valid factor genes found")
  }

  if (smooth_window %% 2 == 0) {
    smooth_window <- smooth_window + 1L
  }

  if (verbose) {
    message("Computing all 4 delta types...")
    message("  Factors: ", paste(factor_genes, collapse = ", "))
    message("  ", nrow(expr_matrix), " genes x ", ncol(expr_matrix), " spots")
  }

  result <- cpp_compute_four_deltas_multi(
    expr_matrix = expr_matrix,
    gene_names = gene_names,
    factor_names = factor_genes,
    coords = coords,
    radii = as.numeric(radii),
    smooth_window = as.integer(smooth_window),
    smooth_poly = as.integer(smooth_poly),
    verbose = verbose
  )

  if (length(factor_genes) == 1) {
    return(result[[1]])
  }

  return(result)
}


# =============================================================================
# Internal helper functions
# =============================================================================

#' @keywords internal
.prepare_coords <- function(spot_coords, coord_scale = 1) {
  if (is.data.frame(spot_coords)) {
    if ("row" %in% names(spot_coords) && "col" %in% names(spot_coords)) {
      coords <- as.matrix(spot_coords[, c("row", "col")])
    } else {
      coords <- as.matrix(spot_coords[, 1:2])
    }
  } else if (is.matrix(spot_coords)) {
    coords <- spot_coords[, 1:2, drop = FALSE]
  } else {
    stop("spot_coords must be a data frame or matrix")
  }
  coords * coord_scale
}

#' @keywords internal
.get_gene_index <- function(gene, gene_names) {
  if (is.character(gene)) {
    idx <- which(gene_names == gene)
    if (length(idx) == 0) {
      stop("Gene '", gene, "' not found")
    }
    idx[1]
  } else {
    as.integer(gene)
  }
}
