#' Compute Spatial Signature (Auto-Detect Platform)
#'
#' High-level wrapper that automatically detects the platform type and calls
#' the appropriate analysis function.
#'
#' @param data Data object (VisiumData or SCData)
#' @param factor_gene Character. Name of the factor gene (e.g., "IL1B")
#' @param ... Additional arguments passed to platform-specific function
#'
#' @return VisiumSignature or ScSignature data.frame
#'
#' @details
#' This function inspects the data object and calls either
#' \code{\link{compute_signature_visium}} for Visium data or
#' \code{\link{compute_signature_sc}} for single-cell data.
#'
#' @examples
#' \dontrun{
#' # Visium data - auto-detected
#' visium_data <- load_data_visium("path/to/spaceranger")
#' sig <- compute_signature(visium_data, "IL1B")
#'
#' # Single-cell data - auto-detected
#' sc_data <- load_data_cosmx("expr.csv", "meta.csv")
#' sig <- compute_signature(sc_data, "IL1B",
#'                          sender_type = "Macrophage",
#'                          receiver_type = "Fibroblast")
#' }
#'
#' @seealso \code{\link{compute_signature_visium}}, \code{\link{compute_signature_sc}}
#'
#' @export
compute_signature <- function(data, factor_gene, ...) {
    if (inherits(data, "VisiumData")) {
        compute_signature_visium(data, factor_gene, ...)
    } else if (inherits(data, "SCData")) {
        compute_signature_sc(data, factor_gene, ...)
    } else if (is.list(data)) {
        # Try to auto-detect based on list contents
        if ("cell_types" %in% names(data)) {
            compute_signature_sc(data, factor_gene, ...)
        } else if ("spot_ids" %in% names(data) || "n_spots" %in% names(data)) {
            compute_signature_visium(data, factor_gene, ...)
        } else {
            stop("Cannot auto-detect platform. Use compute_signature_visium() or ",
                 "compute_signature_sc() directly.")
        }
    } else {
        stop("Unknown data type. Expected VisiumData, SCData, or compatible list.")
    }
}

#' Load Spatial Data (Auto-Detect Format)
#'
#' Convenience function to load spatial transcriptomics data from various formats.
#'
#' @param path Character. Path to data file or directory.
#' @param format Character. Data format (auto-detected if NULL):
#'   "spaceranger", "h5ad", "cosmx", "seurat"
#' @param ... Additional arguments passed to format-specific loader
#'
#' @return VisiumData or SCData object
#'
#' @examples
#' \dontrun{
#' # Auto-detect Space Ranger output
#' data <- load_data("path/to/spaceranger/outs")
#'
#' # Load AnnData file
#' data <- load_data("data.h5ad", format = "h5ad")
#' }
#'
#' @export
load_data <- function(path, format = NULL, ...) {
    if (is.null(format)) {
        # Auto-detect format
        if (dir.exists(path)) {
            if (file.exists(file.path(path, "filtered_feature_bc_matrix"))) {
                format <- "spaceranger"
            } else {
                stop("Cannot auto-detect format for directory: ", path)
            }
        } else if (file.exists(path)) {
            ext <- tolower(tools::file_ext(path))
            if (ext == "h5ad") {
                format <- "h5ad"
            } else {
                stop("Cannot auto-detect format for file: ", path)
            }
        } else {
            stop("Path does not exist: ", path)
        }
    }

    switch(format,
        spaceranger = load_data_visium(path, ...),
        h5ad = load_data_anndata(path, ...),
        stop("Unknown format: ", format)
    )
}

#' Create Spatial Weight Matrix (Auto-Detect Platform)
#'
#' Creates an appropriate spatial weight matrix based on data type.
#'
#' @param data Data object or coordinates matrix
#' @param radius Numeric. Distance threshold
#' @param platform Character. "visium" for binary weights, "sc" for Gaussian
#'   (auto-detected if data is VisiumData/SCData)
#' @param ... Additional arguments passed to platform-specific function
#'
#' @return Sparse weight matrix
#'
#' @export
create_weights <- function(data, radius, platform = NULL, ...) {
    if (is.null(platform)) {
        if (inherits(data, "VisiumData")) {
            platform <- "visium"
        } else if (inherits(data, "SCData")) {
            platform <- "sc"
        } else {
            stop("Cannot auto-detect platform. Specify platform = 'visium' or 'sc'")
        }
    }

    if (platform == "visium") {
        coords <- if (is.matrix(data)) data else data$coords
        create_weights_visium(coords, radius, ...)
    } else if (platform == "sc") {
        coords <- if (is.matrix(data)) data else data$coords
        create_weights_sc(coords, coords, radius, ...)
    } else {
        stop("Unknown platform: ", platform)
    }
}
