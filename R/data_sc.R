#' Load CosMx Data
#'
#' Load single-cell spatial transcriptomics data from CosMx CSV files.
#'
#' @param expr_file Expression file (CSV or RDS). Rows = genes, cols = cells.
#' @param meta_file Metadata CSV with columns: cell_id, x, y, cell_type
#' @return SCData list with components: expr, coords, cell_types, cell_ids,
#'   gene_names, n_cells, n_genes
#' @export
load_data_cosmx <- function(expr_file, meta_file) {
    # Load expression
    if (endsWith(expr_file, ".rds") || endsWith(expr_file, ".RDS")) {
        expr <- readRDS(expr_file)
    } else if (endsWith(expr_file, ".csv") || endsWith(expr_file, ".CSV")) {
        expr <- as.matrix(read.csv(expr_file, row.names = 1, check.names = FALSE))
    } else {
        stop("Unsupported expression file format. Use .rds or .csv")
    }

    # Load metadata
    meta <- read.csv(meta_file, stringsAsFactors = FALSE)

    # Validate required columns
    required_cols <- c("cell_id", "x", "y", "cell_type")
    missing <- setdiff(required_cols, colnames(meta))
    if (length(missing) > 0) {
        stop("Metadata missing columns: ", paste(missing, collapse = ", "))
    }

    # Match order
    meta <- meta[match(colnames(expr), meta$cell_id), ]

    structure(
        list(
            expr = expr,
            coords = as.matrix(meta[, c("x", "y")]),
            cell_types = as.character(meta$cell_type),
            cell_ids = meta$cell_id,
            gene_names = rownames(expr),
            n_cells = ncol(expr),
            n_genes = nrow(expr)
        ),
        class = c("SCData", "list")
    )
}

#' Load from AnnData H5AD
#'
#' Load single-cell spatial data from Python AnnData H5AD file.
#'
#' @param h5ad_file Path to .h5ad file
#' @param cell_type_col Column name for cell type in obs (default: "cell_type")
#' @param spatial_key Key in obsm for coordinates (default: "spatial")
#' @return SCData list
#' @export
load_data_anndata <- function(h5ad_file, cell_type_col = "cell_type",
                               spatial_key = "spatial") {
    if (!requireNamespace("rhdf5", quietly = TRUE)) {
        stop("Package 'rhdf5' required. Install from Bioconductor: ",
             "BiocManager::install('rhdf5')")
    }

    # Read components
    X <- rhdf5::h5read(h5ad_file, "X")
    obs_names <- rhdf5::h5read(h5ad_file, "obs/_index")
    var_names <- rhdf5::h5read(h5ad_file, "var/_index")
    cell_types <- rhdf5::h5read(h5ad_file, paste0("obs/", cell_type_col))
    coords <- rhdf5::h5read(h5ad_file, paste0("obsm/", spatial_key))

    # Handle sparse matrix
    if (is.list(X)) {
        X <- Matrix::sparseMatrix(
            i = X$indices + 1L,
            p = X$indptr,
            x = X$data,
            dims = c(X$shape[1], X$shape[2])
        )
    }

    # Convert to genes x cells (transpose)
    expr <- t(as.matrix(X))
    rownames(expr) <- var_names
    colnames(expr) <- obs_names

    rhdf5::H5close()

    structure(
        list(
            expr = expr,
            coords = coords,
            cell_types = as.character(cell_types),
            cell_ids = obs_names,
            gene_names = var_names,
            n_cells = ncol(expr),
            n_genes = nrow(expr)
        ),
        class = c("SCData", "list")
    )
}

#' Convert Seurat Object to SCData
#'
#' Convert a Seurat object with spatial data to SCData format.
#'
#' @param seurat_obj Seurat object with spatial data
#' @param assay Assay name (default: "RNA")
#' @param cell_type_col Column in metadata for cell type (default: "cell_type")
#' @return SCData list
#' @export
as_data_sc <- function(seurat_obj, assay = "RNA", cell_type_col = "cell_type") {
    if (!requireNamespace("Seurat", quietly = TRUE)) {
        stop("Package 'Seurat' required for this function")
    }

    expr <- Seurat::GetAssayData(seurat_obj, assay = assay, layer = "data")
    coords <- Seurat::GetTissueCoordinates(seurat_obj)

    if (!cell_type_col %in% colnames(seurat_obj@meta.data)) {
        stop("Cell type column '", cell_type_col, "' not found in metadata")
    }

    structure(
        list(
            expr = as.matrix(expr),
            coords = as.matrix(coords[, 1:2]),
            cell_types = as.character(seurat_obj@meta.data[[cell_type_col]]),
            cell_ids = colnames(expr),
            gene_names = rownames(expr),
            n_cells = ncol(expr),
            n_genes = nrow(expr)
        ),
        class = c("SCData", "list")
    )
}

#' Print Method for SCData
#'
#' @param x An SCData object.
#' @param ... Additional arguments (ignored).
#' @return Invisibly returns the input object.
#' @export
print.SCData <- function(x, ...) {
    cat("SCData object (Single-Cell Spatial)\n")
    cat("  Cells:", x$n_cells, "\n")
    cat("  Genes:", x$n_genes, "\n")
    cat("  Cell types:", length(unique(x$cell_types)), "\n")

    # Show cell type counts
    type_counts <- sort(table(x$cell_types), decreasing = TRUE)
    cat("  Top cell types:\n")
    for (i in seq_len(min(5, length(type_counts)))) {
        cat("    ", names(type_counts)[i], ": ", type_counts[i], "\n", sep = "")
    }

    invisible(x)
}
