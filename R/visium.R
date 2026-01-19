#' Load Visium Data from Space Ranger Output
#'
#' Load 10x Visium spatial transcriptomics data from Space Ranger output directory.
#'
#' @param spaceranger_dir Path to the Space Ranger output directory (e.g., "sample/outs"
#'   or the directory containing "filtered_feature_bc_matrix" and "spatial" folders).
#' @param use_filtered Whether to use filtered (TRUE) or raw (FALSE) feature barcode matrix.
#'   Default: TRUE.
#'
#' @return A list with components:
#'   \item{counts}{Sparse matrix (dgCMatrix) of gene counts. Rows are genes, columns are spots.}
#'   \item{features}{Data frame with gene information (id, name, type).}
#'   \item{barcodes}{Character vector of spot barcodes.}
#'   \item{spot_coords}{Data frame with spot coordinates (barcode, in_tissue, row, col, pixel_row, pixel_col).}
#'   \item{in_tissue_idx}{Integer vector of indices for in-tissue spots.}
#'
#' @examples
#' \dontrun{
#' # Load Visium data
#' visium <- load_visium_data("path/to/spaceranger_output")
#'
#' # Get in-tissue counts only
#' counts_tissue <- visium$counts[, visium$in_tissue_idx]
#' coords_tissue <- visium$spot_coords[visium$in_tissue_idx, ]
#' }
#'
#' @export
load_visium_data <- function(spaceranger_dir, use_filtered = TRUE) {

  # Determine matrix directory
  matrix_subdir <- if (use_filtered) "filtered_feature_bc_matrix" else "raw_feature_bc_matrix"

  # Check for different possible directory structures
  if (dir.exists(file.path(spaceranger_dir, matrix_subdir))) {
    matrix_dir <- file.path(spaceranger_dir, matrix_subdir)
    spatial_dir <- file.path(spaceranger_dir, "spatial")
  } else if (dir.exists(file.path(spaceranger_dir, "outs", matrix_subdir))) {
    matrix_dir <- file.path(spaceranger_dir, "outs", matrix_subdir)
    spatial_dir <- file.path(spaceranger_dir, "outs", "spatial")
  } else {
    stop("Cannot find ", matrix_subdir, " directory in ", spaceranger_dir)
  }

  if (!dir.exists(spatial_dir)) {
    stop("Cannot find spatial directory in ", spaceranger_dir)
  }

  # Load matrix files
  matrix_file <- file.path(matrix_dir, "matrix.mtx.gz")
  features_file <- file.path(matrix_dir, "features.tsv.gz")
  barcodes_file <- file.path(matrix_dir, "barcodes.tsv.gz")

  # Check if files exist (also check for uncompressed versions)
  if (!file.exists(matrix_file)) {
    matrix_file <- file.path(matrix_dir, "matrix.mtx")
    if (!file.exists(matrix_file)) {
      stop("Cannot find matrix.mtx or matrix.mtx.gz in ", matrix_dir)
    }
  }

  if (!file.exists(features_file)) {
    features_file <- file.path(matrix_dir, "features.tsv")
    if (!file.exists(features_file)) {
      # Try genes.tsv for older versions
      features_file <- file.path(matrix_dir, "genes.tsv.gz")
      if (!file.exists(features_file)) {
        features_file <- file.path(matrix_dir, "genes.tsv")
        if (!file.exists(features_file)) {
          stop("Cannot find features.tsv or genes.tsv in ", matrix_dir)
        }
      }
    }
  }

  if (!file.exists(barcodes_file)) {
    barcodes_file <- file.path(matrix_dir, "barcodes.tsv")
    if (!file.exists(barcodes_file)) {
      stop("Cannot find barcodes.tsv or barcodes.tsv.gz in ", matrix_dir)
    }
  }

  # Read matrix (Matrix Market format)
  counts <- Matrix::readMM(matrix_file)
  counts <- methods::as(counts, "dgCMatrix")

  # Read features/genes
  features <- read.table(features_file, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  if (ncol(features) >= 3) {
    colnames(features) <- c("id", "name", "type")
  } else if (ncol(features) == 2) {
    colnames(features) <- c("id", "name")
    features$type <- "Gene Expression"
  } else {
    colnames(features) <- "id"
    features$name <- features$id
    features$type <- "Gene Expression"
  }

  # Read barcodes
  barcodes <- read.table(barcodes_file, sep = "\t", header = FALSE, stringsAsFactors = FALSE)$V1

  # Set matrix dimnames
  rownames(counts) <- features$name
  colnames(counts) <- barcodes

  # Read spatial coordinates
  positions_file <- file.path(spatial_dir, "tissue_positions_list.csv")
  if (!file.exists(positions_file)) {
    positions_file <- file.path(spatial_dir, "tissue_positions.csv")
    if (!file.exists(positions_file)) {
      stop("Cannot find tissue_positions_list.csv or tissue_positions.csv in ", spatial_dir)
    }
  }

  # tissue_positions_list.csv format: barcode, in_tissue, array_row, array_col, pxl_row, pxl_col
  spot_coords <- read.csv(positions_file, header = FALSE, stringsAsFactors = FALSE)
  colnames(spot_coords) <- c("barcode", "in_tissue", "row", "col", "pixel_row", "pixel_col")

  # Match barcodes to get coordinates in same order as matrix columns
  barcode_order <- match(barcodes, spot_coords$barcode)
  if (any(is.na(barcode_order))) {
    warning("Some barcodes in matrix not found in tissue positions file")
    barcode_order[is.na(barcode_order)] <- 1  # Placeholder, will be filtered later
  }
  spot_coords <- spot_coords[barcode_order, ]
  rownames(spot_coords) <- NULL

  # Get indices of in-tissue spots
  in_tissue_idx <- which(spot_coords$in_tissue == 1)

  return(list(
    counts = counts,
    features = features,
    barcodes = barcodes,
    spot_coords = spot_coords,
    in_tissue_idx = in_tissue_idx
  ))
}


#' Filter Visium Data to In-Tissue Spots
#'
#' Subset Visium data to include only spots that are within the tissue.
#'
#' @param visium_data A list returned by \code{load_visium_data}.
#'
#' @return A list with the same structure as input, but filtered to in-tissue spots only.
#'
#' @examples
#' \dontrun{
#' visium <- load_visium_data("path/to/spaceranger_output")
#' visium_tissue <- filter_in_tissue(visium)
#' }
#'
#' @export
filter_in_tissue <- function(visium_data) {
  idx <- visium_data$in_tissue_idx

  list(
    counts = visium_data$counts[, idx, drop = FALSE],
    features = visium_data$features,
    barcodes = visium_data$barcodes[idx],
    spot_coords = visium_data$spot_coords[idx, , drop = FALSE],
    in_tissue_idx = seq_along(idx)  # All spots are now in-tissue
  )
}


#' Get Spot Coordinates Data Frame
#'
#' Extract spot coordinates in the format expected by pairwise_moran.
#'
#' @param visium_data A list returned by \code{load_visium_data} or \code{filter_in_tissue}.
#'
#' @return A data frame with columns 'row' and 'col' for spot array coordinates.
#'
#' @export
get_spot_coords <- function(visium_data) {
  data.frame(
    row = visium_data$spot_coords$row,
    col = visium_data$spot_coords$col
  )
}
