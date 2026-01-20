#' VST Normalization for Spatial Transcriptomics Data
#'
#' Apply variance stabilizing transformation (VST) using sctransform to normalize
#' gene expression data. This is the recommended normalization for Moran's I computation.
#'
#' @param counts A matrix or sparse matrix of raw counts. Rows are genes, columns are spots.
#' @param min_cells Minimum number of cells expressing a gene for it to be included.
#'   Default: 5.
#' @param n_genes Number of genes to use for fitting the model. Default: 2000.
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
vst_transform <- function(counts, min_cells = 5, n_genes = 2000, verbose = TRUE, ...) {

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

  # Pre-filter genes: remove genes with zero variance or all zeros
  gene_sums <- rowSums(counts_mat)
  gene_nonzero <- rowSums(counts_mat > 0)

  # Keep genes that have at least min_cells non-zero values and non-zero total
  keep_genes <- (gene_nonzero >= min_cells) & (gene_sums > 0)

  if (verbose) {
    cat("Filtering genes: keeping", sum(keep_genes), "of", length(keep_genes), "genes\n")
  }

  counts_filtered <- counts_mat[keep_genes, , drop = FALSE]

  # Increase future.globals.maxSize for large datasets
  old_option <- getOption("future.globals.maxSize")
  options(future.globals.maxSize = 2 * 1024^3)  # 2 GB
  on.exit(options(future.globals.maxSize = old_option), add = TRUE)

  # Apply VST with error handling
  tryCatch({
    vst_result <- sctransform::vst(
      counts_filtered,
      min_cells = min_cells,
      n_genes = min(n_genes, nrow(counts_filtered)),
      verbosity = if(verbose) 1 else 0,
      ...
    )

    # Extract normalized values (Pearson residuals)
    data_vst <- vst_result$y

    if (verbose) {
      cat("Output:", nrow(data_vst), "genes x", ncol(data_vst), "spots\n")
      cat("VST normalization complete.\n")
    }

    return(data_vst)

  }, error = function(e) {
    stop("VST transformation failed: ", e$message, "\n",
         "Try adjusting min_cells or n_genes parameters.")
  })
}


#' Full Pipeline: Load, Preprocess, and Compute Moran's I
#'
#' Convenience function that runs the complete analysis pipeline from raw
#' Space Ranger output to pairwise Moran's I results.
#'
#' @param spaceranger_dir Path to the Space Ranger output directory.
#' @param min_cells Minimum number of cells for VST. Default: 5.
#' @param n_genes Number of genes to use for VST model fitting. Default: 2000.
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
                         n_genes = 2000,
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
  data_vst <- vst_transform(visium$counts, min_cells = min_cells, n_genes = n_genes, verbose = verbose)

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
#' @param result Matrix of pairwise Moran's I values, or a list containing a
#'   'moran' element (as returned by \code{run_pipeline}).
#' @param output_file Path to output file.
#' @param verbose Print progress messages. Default: TRUE.
#'
#' @return Invisibly returns the output file path.
#'
#' @export
save_moran_result <- function(result, output_file, verbose = TRUE) {
  if (verbose) cat("Saving result to", output_file, "...\n")

  # Handle list input (from run_pipeline)
  if (is.list(result) && "moran" %in% names(result)) {
    moran_matrix <- result$moran
  } else {
    moran_matrix <- result
  }

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
#' @param gene_names Optional character vector of gene names for row/column names.
#'   If NULL, tries to read from a file with same name but .genes extension.
#'
#' @return A symmetric matrix of pairwise Moran's I values with row/column names.
#'
#' @export
load_moran_result <- function(input_file, gene_names = NULL) {
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

  # Set row/column names
  if (is.null(gene_names)) {
    # Try to read from .genes file
    genes_file <- sub("\\.[^.]+$", ".genes", input_file)
    if (file.exists(genes_file)) {
      gene_names <- readLines(genes_file)
    }
  }

  if (!is.null(gene_names) && length(gene_names) == n) {
    rownames(mat) <- gene_names
    colnames(mat) <- gene_names
  }

  return(mat)
}
