#' VST Normalization for Spatial Transcriptomics Data
#'
#' Apply variance stabilizing transformation (VST) using sctransform to normalize
#' gene expression data. This is the recommended normalization for Moran's I computation.
#'
#' @param counts A matrix or sparse matrix of raw counts. Rows are genes, columns are spots.
#' @param min_cells Minimum number of cells expressing a gene for it to be included.
#'   Default: 5.
#' @param verbose Print progress messages. Default: TRUE.
#' @param ... Additional arguments passed to \code{sctransform::vst}.
#'
#' @return A dense matrix of VST-normalized values. Rows are genes, columns are spots.
#'   Note: genes that don't meet the min_cells threshold are excluded.
#'
#' @details
#' This function requires the \code{sctransform} package to be installed.
#' Install it with: \code{install.packages("sctransform")}
#'
#' The VST normalization:
#' \enumerate{
#'   \item Models the relationship between mean expression and variance
#'   \item Applies a variance stabilizing transformation
#'   \item Returns Pearson residuals that are approximately normally distributed
#' }
#'
#' @examples
#' \dontrun{
#' # Load Visium data
#' visium <- load_visium_data("path/to/spaceranger_output")
#' visium <- filter_in_tissue(visium)
#'
#' # Apply VST normalization
#' data_vst <- vst_transform(visium$counts)
#'
#' # Compute pairwise Moran's I
#' coords <- get_spot_coords(visium)
#' result <- pairwise_moran(data_vst, coords)
#' }
#'
#' @export
vst_transform <- function(counts, min_cells = 5, verbose = TRUE, ...) {

  # Check if sctransform is available
if (!requireNamespace("sctransform", quietly = TRUE)) {
    stop("Package 'sctransform' is required for VST normalization.\n",
         "Install it with: install.packages('sctransform')")
  }

  if (verbose) {
    cat("Applying VST normalization...\n")
    cat("Input:", nrow(counts), "genes x", ncol(counts), "spots\n")
  }

  # Convert sparse to dense if needed (sctransform handles both)
  if (methods::is(counts, "sparseMatrix")) {
    counts_mat <- as.matrix(counts)
  } else {
    counts_mat <- as.matrix(counts)
  }

  # Apply VST
  vst_result <- sctransform::vst(counts_mat, min_cells = min_cells, verbosity = if(verbose) 1 else 0, ...)

  # Extract normalized values (Pearson residuals)
  data_vst <- vst_result$y

  if (verbose) {
    cat("Output:", nrow(data_vst), "genes x", ncol(data_vst), "spots\n")
    cat("VST normalization complete.\n")
  }

  return(data_vst)
}


#' Full Pipeline: Load, Preprocess, and Compute Moran's I
#'
#' Convenience function that runs the complete analysis pipeline from raw
#' Space Ranger output to pairwise Moran's I results.
#'
#' @param spaceranger_dir Path to the Space Ranger output directory.
#' @param min_cells Minimum number of cells for VST. Default: 5.
#' @param max_radius Maximum grid radius for Moran's I computation. Default: 5.
#' @param platform Platform type: "visium" (default) or "old".
#' @param same_spot Whether to consider the same spot in computation. Default: TRUE.
#' @param mode Computation mode: "paired" (all pairs), "first" (pairs with first gene),
#'   or "single" (diagonal only). Default: "paired".
#' @param verbose Print progress messages. Default: TRUE.
#'
#' @return A list with components:
#'   \item{moran}{Matrix of pairwise Moran's I values.}
#'   \item{data_vst}{VST-normalized expression matrix.}
#'   \item{spot_coords}{Data frame of spot coordinates.}
#'   \item{gene_names}{Character vector of gene names.}
#'   \item{n_genes}{Number of genes after filtering.}
#'   \item{n_spots}{Number of in-tissue spots.}
#'
#' @examples
#' \dontrun{
#' # Run complete pipeline
#' result <- run_pipeline("path/to/spaceranger_output")
#'
#' # Access results
#' moran_matrix <- result$moran
#' gene_names <- result$gene_names
#' }
#'
#' @export
run_pipeline <- function(spaceranger_dir,
                         min_cells = 5,
                         max_radius = 5,
                         platform = c("visium", "old"),
                         same_spot = TRUE,
                         mode = c("paired", "first", "single"),
                         verbose = TRUE) {

  platform <- match.arg(platform)
  mode <- match.arg(mode)

  # Step 1: Load data
  if (verbose) cat("=== Step 1: Loading Visium data ===\n")
  visium <- load_visium_data(spaceranger_dir)
  visium <- filter_in_tissue(visium)

  if (verbose) {
    cat("Loaded", nrow(visium$counts), "genes x", ncol(visium$counts), "spots\n\n")
  }

  # Step 2: VST normalization
  if (verbose) cat("=== Step 2: VST normalization ===\n")
  data_vst <- vst_transform(visium$counts, min_cells = min_cells, verbose = verbose)

  if (verbose) cat("\n")

  # Step 3: Get spot coordinates
  spot_coords <- get_spot_coords(visium)

  # Step 4: Compute pairwise Moran's I
  if (verbose) cat("=== Step 3: Computing pairwise Moran's I ===\n")
  moran_result <- pairwise_moran(
    data_vst,
    spot_coords,
    max_radius = max_radius,
    platform = platform,
    same_spot = same_spot,
    mode = mode,
    verbose = verbose
  )

  if (verbose) cat("\n=== Pipeline complete ===\n")

  return(list(
    moran = moran_result,
    data_vst = data_vst,
    spot_coords = spot_coords,
    gene_names = rownames(data_vst),
    n_genes = nrow(data_vst),
    n_spots = ncol(data_vst)
  ))
}


#' Save Moran's I Results
#'
#' Save pairwise Moran's I results to a file in lower triangular format.
#'
#' @param moran_matrix Matrix of pairwise Moran's I values.
#' @param output_file Path to output file.
#' @param verbose Print progress messages. Default: TRUE.
#'
#' @return Invisibly returns the output file path.
#'
#' @export
save_moran_result <- function(moran_matrix, output_file, verbose = TRUE) {
  if (verbose) cat("Saving result to", output_file, "...\n")

  n <- nrow(moran_matrix)
  con <- file(output_file, "w")

  for (i in 1:n) {
    row_values <- moran_matrix[i, 1:i]
    writeLines(paste(row_values, collapse = "\t"), con)
  }

  close(con)

  if (verbose) cat("Done.\n")

  invisible(output_file)
}


#' Load Moran's I Results
#'
#' Load pairwise Moran's I results from a file in lower triangular format.
#'
#' @param input_file Path to input file.
#'
#' @return A symmetric matrix of pairwise Moran's I values.
#'
#' @export
load_moran_result <- function(input_file) {
  lines <- readLines(input_file)
  n <- length(lines)
  mat <- matrix(0, nrow = n, ncol = n)

  for (i in 1:n) {
    vals <- as.numeric(strsplit(trimws(lines[i]), "\\s+")[[1]])
    mat[i, 1:i] <- vals
    if (i > 1) {
      mat[1:(i-1), i] <- vals[1:(i-1)]
    }
  }

  return(mat)
}
