#' Load Visium Data from Space Ranger Output
#'
#' Loads expression matrix and spatial coordinates from a Space Ranger
#' output directory.
#'
#' @param path Character. Path to Space Ranger output directory.
#' @param use_image_coords Logical. Use image (pixel) coordinates instead of
#'   array coordinates (default: FALSE).
#' @param scale_factor Numeric. Scale factor for array coordinates to convert
#'   to approximate micron scale (default: 100). Ignored if use_image_coords
#'   is TRUE.
#'
#' @return A VisiumData list with components:
#' \describe{
#'   \item{expr}{Numeric matrix. Expression values (genes x spots).}
#'   \item{coords}{Numeric matrix. Spot coordinates (n_spots x 2).}
#'   \item{spot_ids}{Character vector. Spot barcodes.}
#'   \item{gene_names}{Character vector. Gene names/IDs.}
#'   \item{n_spots}{Integer. Number of spots.}
#'   \item{n_genes}{Integer. Number of genes.}
#' }
#'
#' @details
#' Space Ranger outputs two types of coordinates:
#' \itemize{
#'   \item Array coordinates (row, col): Grid positions, typically with ~100um
#'     spacing. Default option, scaled by scale_factor.
#'   \item Image coordinates (pxl_row, pxl_col): Pixel positions in the tissue
#'     image. May need scaling depending on image resolution.
#' }
#'
#' Requires the Seurat package for reading 10X matrix format.
#'
#' @examples
#' \dontrun{
#' # Load from Space Ranger output
#' data <- load_data_visium("path/to/spaceranger/outs")
#'
#' # Use image coordinates
#' data <- load_data_visium("path/to/spaceranger/outs", use_image_coords = TRUE)
#' }
#'
#' @seealso \code{\link{as_data_visium}} for converting Seurat objects
#'
#' @export
load_data_visium <- function(path, use_image_coords = FALSE, scale_factor = 100) {
    # Check directory structure
    matrix_dir <- file.path(path, "filtered_feature_bc_matrix")
    spatial_dir <- file.path(path, "spatial")

    if (!dir.exists(matrix_dir)) {
        stop("Matrix directory not found: ", matrix_dir)
    }
    if (!dir.exists(spatial_dir)) {
        stop("Spatial directory not found: ", spatial_dir)
    }

    # Load expression matrix
    if (requireNamespace("Seurat", quietly = TRUE)) {
        expr <- Seurat::Read10X(matrix_dir)
    } else {
        stop("Seurat package required for loading 10X data. ",
             "Install with: install.packages('Seurat')")
    }

    # Load spatial coordinates
    tissue_file <- file.path(spatial_dir, "tissue_positions_list.csv")
    if (!file.exists(tissue_file)) {
        tissue_file <- file.path(spatial_dir, "tissue_positions.csv")
    }
    if (!file.exists(tissue_file)) {
        stop("Tissue positions file not found in: ", spatial_dir)
    }

    tissue_positions <- utils::read.csv(tissue_file, header = FALSE, row.names = 1)
    colnames(tissue_positions) <- c("in_tissue", "array_row", "array_col",
                                     "pxl_row", "pxl_col")

    # Filter to spots in tissue
    in_tissue <- tissue_positions$in_tissue == 1
    tissue_positions <- tissue_positions[in_tissue, ]

    # Match spots between expression and coordinates
    common_spots <- intersect(colnames(expr), rownames(tissue_positions))
    if (length(common_spots) == 0) {
        stop("No matching spots between expression matrix and tissue positions")
    }

    expr <- expr[, common_spots]
    tissue_positions <- tissue_positions[common_spots, ]

    # Select coordinate type
    if (use_image_coords) {
        coords <- as.matrix(tissue_positions[, c("pxl_row", "pxl_col")])
    } else {
        coords <- as.matrix(tissue_positions[, c("array_row", "array_col")])
        coords <- coords * scale_factor
    }

    structure(
        list(
            expr = as.matrix(expr),
            coords = coords,
            spot_ids = colnames(expr),
            gene_names = rownames(expr),
            n_spots = ncol(expr),
            n_genes = nrow(expr)
        ),
        class = c("VisiumData", "list")
    )
}

#' Convert Seurat Object to VisiumData
#'
#' Extracts expression matrix and spatial coordinates from a Seurat object
#' containing Visium data.
#'
#' @param seurat_obj A Seurat object with spatial data.
#' @param assay Character. Assay name (default: "Spatial").
#' @param slot Character. Data slot: "counts", "data", or "scale.data"
#'   (default: "data" for log-normalized data).
#'
#' @return A VisiumData list (see \code{\link{load_data_visium}}).
#'
#' @details
#' This function provides a convenient interface for users who have already
#' processed their Visium data in Seurat. The returned VisiumData object
#' can be used directly with \code{\link{compute_signature_visium}}.
#'
#' @examples
#' \dontrun{
#' library(Seurat)
#' seurat_obj <- Load10X_Spatial("path/to/spaceranger")
#' data <- as_data_visium(seurat_obj)
#' }
#'
#' @seealso \code{\link{load_data_visium}}
#'
#' @export
as_data_visium <- function(seurat_obj, assay = "Spatial", slot = "data") {
    if (!requireNamespace("Seurat", quietly = TRUE)) {
        stop("Seurat package required for this function")
    }

    # Extract expression data
    expr <- Seurat::GetAssayData(seurat_obj, assay = assay, layer = slot)

    # Extract coordinates
    coords <- Seurat::GetTissueCoordinates(seurat_obj)

    structure(
        list(
            expr = as.matrix(expr),
            coords = as.matrix(coords[, 1:2]),
            spot_ids = colnames(expr),
            gene_names = rownames(expr),
            n_spots = ncol(expr),
            n_genes = nrow(expr)
        ),
        class = c("VisiumData", "list")
    )
}

#' Print Method for VisiumData
#'
#' @param x A VisiumData object.
#' @param ... Additional arguments (ignored).
#'
#' @return Invisibly returns the input object.
#'
#' @export
print.VisiumData <- function(x, ...) {
    cat("VisiumData object\n")
    cat("  Spots:", x$n_spots, "\n")
    cat("  Genes:", x$n_genes, "\n")
    cat("  Coordinate range:\n")
    cat("    X:", round(range(x$coords[, 1]), 1), "\n")
    cat("    Y:", round(range(x$coords[, 2]), 1), "\n")
    invisible(x)
}
