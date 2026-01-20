#' Compute Signed Delta I Signatures for Spatial Transcriptomics
#'
#' Computes distance-dependent spatial correlation between a factor gene and all
#' other genes, then derives signed delta I signatures. This helps distinguish
#' true responsive genes from constitutively expressed genes.
#'
#' @param expr_matrix A matrix of gene expression values (genes x spots).
#'   VST-normalized values are recommended.
#' @param spot_coords Spot coordinates as a data frame with 'row' and 'col' columns,
#'   or a 2-column matrix. Coordinates should be in physical units (e.g., micrometers)
#'   or use \code{coord_scale} to convert.
#' @param factor_gene Name or row index of the factor gene.
#' @param radii Numeric vector of distance bin edges in same units as coordinates.
#'   Default: \code{seq(100, 600, 100)} for Visium in micrometers.
#' @param mode Computation mode: "bivariate" (symmetric, all spots) or "directional"
#'   (I_ND with sender/receiver split). Default: "bivariate".
#' @param sender_percentile For directional mode: percentile threshold for defining
#'   sender spots. Spots with factor expression >= this percentile are senders.
#'   Default: 75 (top 25% expressors are senders).
#' @param smooth_window Savitzky-Golay smoothing window size (must be odd).
#'   Default: 5.
#' @param smooth_poly Savitzky-Golay polynomial order. Default: 2.
#' @param coord_scale Scale factor for coordinates. For Visium grid coordinates,
#'   use ~100 to convert to micrometers. Default: 1.
#' @param verbose Print progress messages. Default: TRUE.
#'
#' @return A data frame with one row per gene containing: gene (name),
#'   delta_I (unsigned magnitude), delta_I_signed (positive=decay/responder,
#'   negative=increase/avoidance), sign (+1 or -1), I_short (correlation at
#'   shortest distance), I_long (correlation at longest distance), I_max,
#'   I_min, auc (area under curve), peak_radius_idx (radius of max |I|).
#'
#' @details
#' The function computes bivariate Moran's I (or directional I_ND) at multiple
#' distance bins, creating an I(r) curve. This curve is smoothed using a
#' Savitzky-Golay filter, then signed delta I is computed based on the trend.
#'
#' Interpretation: delta_I_signed > 0 indicates decay pattern (TRUE RESPONDER),
#' delta_I_signed < 0 indicates increase pattern (AVOIDANCE), and values near
#' zero indicate flat pattern (CONSTITUTIVE expression).
#'
#' Mode "bivariate" uses all spots symmetrically for general spatial pattern
#' discovery. Mode "directional" uses explicit sender-receiver directionality
#' (I_ND cosine similarity) for ligand-receptor signaling analysis.
#'
#' @examples
#' \dontrun{
#' # Load Visium data
#' visium <- load_visium_data("path/to/spaceranger_output")
#' visium <- filter_in_tissue(visium)
#' expr_vst <- vst_transform(visium$counts)
#' coords <- get_spot_coords(visium)
#'
#' # Bivariate mode (symmetric)
#' sig <- compute_signed_delta_I(
#'   expr_matrix = expr_vst,
#'   spot_coords = coords,
#'   factor_gene = "IL1B",
#'   mode = "bivariate",
#'   coord_scale = 100  # Convert grid to micrometers
#' )
#'
#' # View top responders (decay pattern)
#' head(sig[order(-sig$delta_I_signed), ])
#'
#' # Directional mode (sender -> receiver)
#' sig_dir <- compute_signed_delta_I(
#'   expr_matrix = expr_vst,
#'   spot_coords = coords,
#'   factor_gene = "IL1B",
#'   mode = "directional",
#'   sender_percentile = 75
#' )
#' }
#'
#' @seealso \code{\link{get_moran_curve}} for visualizing individual gene curves,
#'   \code{\link{pairwise_moran}} for standard pairwise Moran's I.
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

  # Input validation
  if (!is.matrix(expr_matrix)) {
    expr_matrix <- as.matrix(expr_matrix)
  }

  # Get factor index
  if (is.character(factor_gene)) {
    factor_idx <- which(rownames(expr_matrix) == factor_gene)
    if (length(factor_idx) == 0) {
      stop("Factor gene '", factor_gene, "' not found in expression matrix")
    }
    if (length(factor_idx) > 1) {
      warning("Multiple matches for '", factor_gene, "', using first match")
      factor_idx <- factor_idx[1]
    }
    factor_idx <- factor_idx - 1L  # Convert to 0-based for C++
  } else {
    factor_idx <- as.integer(factor_gene) - 1L
    if (factor_idx < 0 || factor_idx >= nrow(expr_matrix)) {
      stop("factor_gene index out of range")
    }
  }

  # Prepare coordinates
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
  coords <- coords * coord_scale

  # Validate dimensions
  if (ncol(expr_matrix) != nrow(coords)) {
    stop("Number of columns in expr_matrix (", ncol(expr_matrix),
         ") must match number of spots in spot_coords (", nrow(coords), ")")
  }

  # Gene names
  gene_names <- rownames(expr_matrix)
  if (is.null(gene_names)) {
    gene_names <- paste0("Gene_", seq_len(nrow(expr_matrix)))
  }

  # Ensure window is odd
  if (smooth_window %% 2 == 0) {
    smooth_window <- smooth_window + 1L
    if (verbose) message("Adjusted smooth_window to ", smooth_window, " (must be odd)")
  }

  if (verbose) {
    message("Computing signed delta I signatures...")
    message("  Mode: ", mode)
    message("  Factor gene: ", if (is.character(factor_gene)) factor_gene else gene_names[factor_idx + 1])
    message("  ", nrow(expr_matrix), " genes x ", ncol(expr_matrix), " spots")
    message("  ", length(radii), " distance bins: ", paste(radii, collapse = ", "))
  }

  # Call C++ implementation
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

  # Add metadata as attributes
  attr(result, "mode") <- mode
  attr(result, "factor_gene") <- if (is.character(factor_gene)) factor_gene else gene_names[factor_idx + 1]
  attr(result, "radii") <- radii
  attr(result, "sender_percentile") <- sender_percentile

  return(result)
}


#' Get Moran's I Curve for a Gene Pair
#'
#' Computes the full I(r) curve for a single factor-gene pair, useful for
#' visualization and detailed analysis.
#'
#' @param expr_matrix A matrix of gene expression values (genes x spots).
#' @param spot_coords Spot coordinates as a data frame or matrix.
#' @param factor_gene Name or row index of the factor gene.
#' @param target_gene Name or row index of the target gene.
#' @param radii Numeric vector of distance bin edges. Default: \code{seq(100, 600, 100)}.
#' @param mode Computation mode: "bivariate" or "directional". Default: "bivariate".
#' @param sender_percentile For directional mode: percentile threshold for senders.
#'   Default: 75.
#' @param smooth_window Savitzky-Golay smoothing window size. Default: 5.
#' @param smooth_poly Savitzky-Golay polynomial order. Default: 2.
#' @param coord_scale Scale factor for coordinates. Default: 1.
#'
#' @return A list containing: radii (distance bin edges), I_raw (raw I(r) values),
#'   I_smooth (smoothed values), delta_I, delta_I_signed, sign (+1 or -1),
#'   I_max, I_min, I_short, I_long, auc, peak_idx, factor_gene, target_gene,
#'   mode, and for directional mode: n_senders, n_receivers.
#'
#' @examples
#' \dontrun{
#' # Get curve for a specific gene pair
#' curve <- get_moran_curve(
#'   expr_matrix = expr_vst,
#'   spot_coords = coords,
#'   factor_gene = "IL1B",
#'   target_gene = "COL1A1",
#'   mode = "bivariate",
#'   coord_scale = 100
#' )
#'
#' # Plot the curve
#' plot(curve$radii, curve$I_raw, type = "p", pch = 19,
#'      xlab = "Distance (um)", ylab = "Moran's I",
#'      main = paste(curve$factor_gene, "->", curve$target_gene))
#' lines(curve$radii, curve$I_smooth, col = "blue", lwd = 2)
#' abline(h = 0, lty = 2, col = "gray")
#' legend("topright",
#'        legend = paste("delta_I_signed =", round(curve$delta_I_signed, 4)),
#'        bty = "n")
#' }
#'
#' @seealso \code{\link{compute_signed_delta_I}} for batch processing,
#'   \code{\link{plot_moran_curve}} for visualization.
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

  # Input validation
  if (!is.matrix(expr_matrix)) {
    expr_matrix <- as.matrix(expr_matrix)
  }

  # Gene names
  gene_names <- rownames(expr_matrix)
  if (is.null(gene_names)) {
    gene_names <- paste0("Gene_", seq_len(nrow(expr_matrix)))
  }

  # Get factor index
  if (is.character(factor_gene)) {
    factor_idx <- which(gene_names == factor_gene)
    if (length(factor_idx) == 0) {
      stop("Factor gene '", factor_gene, "' not found in expression matrix")
    }
    factor_idx <- factor_idx[1]
    factor_name <- factor_gene
  } else {
    factor_idx <- as.integer(factor_gene)
    if (factor_idx < 1 || factor_idx > nrow(expr_matrix)) {
      stop("factor_gene index out of range")
    }
    factor_name <- gene_names[factor_idx]
  }

  # Get target index
  if (is.character(target_gene)) {
    target_idx <- which(gene_names == target_gene)
    if (length(target_idx) == 0) {
      stop("Target gene '", target_gene, "' not found in expression matrix")
    }
    target_idx <- target_idx[1]
    target_name <- target_gene
  } else {
    target_idx <- as.integer(target_gene)
    if (target_idx < 1 || target_idx > nrow(expr_matrix)) {
      stop("target_gene index out of range")
    }
    target_name <- gene_names[target_idx]
  }

  # Prepare coordinates
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
  coords <- coords * coord_scale

  # Ensure window is odd
  if (smooth_window %% 2 == 0) {
    smooth_window <- smooth_window + 1L
  }

  # Get expression vectors
  factor_expr <- as.numeric(expr_matrix[factor_idx, ])
  gene_expr <- as.numeric(expr_matrix[target_idx, ])

  # Call C++ implementation
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

  # Add gene names and mode
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
#' @param title Custom plot title. If NULL, auto-generated from gene names.
#' @param ... Additional arguments passed to \code{plot()}.
#'
#' @return Invisibly returns the curve_result.
#'
#' @examples
#' \dontrun{
#' curve <- get_moran_curve(expr_vst, coords, "IL1B", "COL1A1")
#' plot_moran_curve(curve)
#' }
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

  # Auto-generate title
  if (is.null(title)) {
    sign_text <- if (curve_result$sign > 0) "decay" else "increase"
    title <- sprintf("%s -> %s\ndelta_I_signed = %.4f (%s)",
                     curve_result$factor_gene,
                     curve_result$target_gene,
                     curve_result$delta_I_signed,
                     sign_text)
  }

  # Determine y-axis limits
  all_vals <- c(I_raw, I_smooth)
  all_vals <- all_vals[is.finite(all_vals)]
  if (length(all_vals) == 0) {
    ylim <- c(-1, 1)
  } else {
    ylim <- range(all_vals) * 1.1
    if (ylim[1] > 0) ylim[1] <- 0
    if (ylim[2] < 0) ylim[2] <- 0
  }

  # Initialize plot
  plot(radii, I_raw, type = "n",
       xlab = expression(paste("Distance (", mu, "m)")),
       ylab = "Spatial Correlation (I)",
       main = title,
       ylim = ylim,
       ...)

  # Add horizontal line at 0

abline(h = 0, lty = 2, col = "gray50")

  # Add reference lines for I_max and I_min
  abline(h = curve_result$I_max, lty = 3, col = "red", lwd = 0.5)
  abline(h = curve_result$I_min, lty = 3, col = "red", lwd = 0.5)

  # Plot raw points
  if (show_raw) {
    points(radii, I_raw, pch = 19, col = "gray40", cex = 1.2)
  }

  # Plot smoothed line
  if (show_smooth) {
    # Handle NAs for smooth line
    valid_idx <- which(is.finite(I_smooth))
    if (length(valid_idx) > 1) {
      lines(radii[valid_idx], I_smooth[valid_idx], col = "blue", lwd = 2)
    }
  }

  # Add legend
  legend_items <- c()
  legend_cols <- c()
  legend_pch <- c()
  legend_lty <- c()

  if (show_raw) {
    legend_items <- c(legend_items, "Raw I(r)")
    legend_cols <- c(legend_cols, "gray40")
    legend_pch <- c(legend_pch, 19)
    legend_lty <- c(legend_lty, NA)
  }

  if (show_smooth) {
    legend_items <- c(legend_items, "Smoothed")
    legend_cols <- c(legend_cols, "blue")
    legend_pch <- c(legend_pch, NA)
    legend_lty <- c(legend_lty, 1)
  }

  if (length(legend_items) > 0) {
    legend("topright",
           legend = legend_items,
           col = legend_cols,
           pch = legend_pch,
           lty = legend_lty,
           lwd = 2,
           bty = "n")
  }

  invisible(curve_result)
}
