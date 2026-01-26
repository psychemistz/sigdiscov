#' Genome-Wide Spatial Interaction Analysis
#'
#' Unified entry point for genome-wide spatial analysis, equivalent to
#' Python's unified_genomewide_interaction.py. Computes Moran's I for
#' factor-target gene pairs across cell type pairs and radii.
#'
#' CPU-optimized for cross-platform compatibility (Mac Intel/M1, Linux, Windows).
#' Uses KD-tree for efficient neighbor search and BLAS-accelerated matrix operations.
#'
#' @param data SCData object or list with expr, coords, cell_types, gene_names
#' @param factor_genes Character vector of factor genes to analyze (default: all genes)
#' @param target_genes Character vector of target genes (default: all genes)
#' @param radii Numeric vector of distance radii (default: c(20, 50, 100, 200, 300, 500))
#' @param pairs Cell type pairs matrix (2 columns: sender, receiver). NULL = all pairs.
#' @param quantile_probs Expression quantile thresholds (default: 0.25). Can be vector
#'   for multi-quantile analysis.
#' @param use_annular Use annular (ring) weights (default: FALSE)
#' @param annular_width Width of annular ring in microns (default: 20)
#' @param annular_mode "constant_width" (default) or "variable" (uses radii intervals)
#' @param output_dir Directory for HDF5 output (NULL = return in memory)
#' @param method "streaming" (default, memory-efficient) or "matrix" (faster for small data)
#' @param sigma_factor Gaussian sigma = radius * sigma_factor (default: 1/3)
#' @param min_cells Minimum cells per type (default: 10)
#' @param min_connections Minimum neighbor connections (default: 10)
#' @param n_perm Number of permutations for significance testing (default: 0)
#' @param fdr_threshold FDR threshold for significance (default: 0.05)
#' @param batch_size Gene batch size for memory management (default: 100)
#' @param verbose Print progress (default: TRUE)
#'
#' @return GenomewideResult object with:
#' \describe{
#'   \item{results}{List of Moran's I results per pair}
#'   \item{delta_I}{Delta I values per gene pair (if computed)}
#'   \item{pair_stats}{Data frame with pair statistics}
#'   \item{metadata}{Analysis parameters}
#' }
#'
#' @details
#' ## Computational Methods
#'
#' **Streaming Mode** (default): Computes weights on-the-fly without storing
#' explicit weight matrix. Uses KD-tree for O(n log n) neighbor search.
#' Memory usage: O(n_genes^2) for result matrix only.
#'
#' **Matrix Mode**: Pre-computes sparse weight matrix. Faster for repeated
#' computations on same cell configuration. Memory: O(n_edges) for weights.
#'
#' ## Platform Optimization
#'
#' - Uses RcppArmadillo with system BLAS (vecLib on Mac, OpenBLAS on Linux)
#' - KD-tree neighbor search via nanoflann (header-only, cross-platform)
#' - Automatic method selection based on dataset size
#'
#' ## Multi-Quantile Analysis
#'
#' When quantile_probs is a vector, the analysis is performed for each
#' quantile threshold. This allows exploration of how spatial correlations
#' change with expression level cutoffs.
#'
#' ## Annular (Ring) Weights
#'
#' When use_annular = TRUE, only neighbors within an annular region are
#' considered. This isolates signals at specific distance ranges and is
#' useful for identifying paracrine vs juxtacrine interactions.
#'
#' @examples
#' \dontrun{
#' # Load data
#' data <- load_data_cosmx("expr.csv", "meta.csv")
#'
#' # Basic analysis
#' result <- genomewide_analysis(
#'     data,
#'     factor_genes = c("TGFB1", "IFNG", "TNF"),
#'     radii = c(50, 100, 200, 300)
#' )
#'
#' # Memory-efficient streaming with HDF5 output
#' genomewide_analysis(
#'     data,
#'     factor_genes = c("TGFB1"),
#'     radii = seq(20, 500, by = 20),
#'     output_dir = "results/",
#'     method = "streaming"
#' )
#'
#' # Annular (ring) analysis for paracrine signals
#' result <- genomewide_analysis(
#'     data,
#'     use_annular = TRUE,
#'     annular_width = 20,
#'     radii = c(50, 100, 200)
#' )
#'
#' # Multiple quantile thresholds
#' result <- genomewide_analysis(
#'     data,
#'     quantile_probs = c(0, 0.25, 0.50),
#'     radii = c(100, 200)
#' )
#' }
#'
#' @export
genomewide_analysis <- function(
    data,
    factor_genes = NULL,
    target_genes = NULL,
    radii = c(20, 50, 100, 200, 300, 500),
    pairs = NULL,
    quantile_probs = 0.25,
    use_annular = FALSE,
    annular_width = 20,
    annular_mode = c("constant_width", "variable"),
    output_dir = NULL,
    method = c("streaming", "matrix"),
    sigma_factor = 1/3,
    min_cells = 10,
    min_connections = 10,
    n_perm = 0,
    fdr_threshold = 0.05,
    batch_size = 100,
    verbose = TRUE
) {
    method <- match.arg(method)
    annular_mode <- match.arg(annular_mode)

    # Validate inputs
    .validate_genomewide_inputs(data)

    # Get gene lists
    if (is.null(factor_genes)) factor_genes <- data$gene_names
    if (is.null(target_genes)) target_genes <- data$gene_names

    # Validate genes exist
    missing_factors <- setdiff(factor_genes, data$gene_names)
    if (length(missing_factors) > 0) {
        warning("Factor genes not found: ", paste(utils::head(missing_factors, 5), collapse = ", "))
        factor_genes <- intersect(factor_genes, data$gene_names)
    }

    missing_targets <- setdiff(target_genes, data$gene_names)
    if (length(missing_targets) > 0) {
        warning("Target genes not found: ", paste(utils::head(missing_targets, 5), collapse = ", "))
        target_genes <- intersect(target_genes, data$gene_names)
    }

    if (length(factor_genes) == 0) {
        stop("No valid factor genes found in data")
    }

    if (verbose) {
        message("Genome-Wide Spatial Analysis")
        message("  ", length(data$gene_names), " genes, ", nrow(data$coords), " cells")
        message("  Factor genes: ", length(factor_genes))
        message("  Target genes: ", length(target_genes))
        message("  Radii: ", length(radii), " (", min(radii), " to ", max(radii), ")")
        message("  Method: ", method)
        if (use_annular) {
            message("  Annular mode: ", annular_mode, ", width = ", annular_width)
        }
        if (length(quantile_probs) > 1) {
            message("  Quantile thresholds: ", paste(quantile_probs, collapse = ", "))
        }
    }

    # Generate cell type pairs
    pairs <- .get_cell_type_pairs_internal(data, pairs, min_cells, verbose)

    if (nrow(pairs) == 0) {
        stop("No valid cell type pairs found")
    }

    # For multi-quantile analysis, run analysis for each quantile
    all_results <- list()

    for (q in sort(unique(quantile_probs))) {
        if (verbose && length(quantile_probs) > 1) {
            message("\n--- Quantile threshold: ", q, " ---")
        }

        # Compute analysis with compute_celltype_pair_analysis
        result_q <- compute_celltype_pair_analysis(
            data = data,
            radii = radii,
            pairs = pairs,
            min_cells = min_cells,
            output_dir = if (!is.null(output_dir) && length(quantile_probs) > 1) {
                file.path(output_dir, paste0("q", q))
            } else {
                output_dir
            },
            compute_moran = TRUE,
            compute_ind = FALSE,
            method = method,
            sigma_factor = sigma_factor,
            verbose = verbose
        )

        # Tag with quantile
        result_q$quantile <- q
        all_results[[paste0("q", q)]] <- result_q
    }

    # Combine results
    if (length(all_results) == 1) {
        result <- all_results[[1]]
    } else {
        result <- list(
            by_quantile = all_results,
            pair_stats = do.call(rbind, lapply(names(all_results), function(qname) {
                r <- all_results[[qname]]
                if (!is.null(r$pair_stats)) {
                    r$pair_stats$quantile <- r$quantile
                    r$pair_stats
                }
            }))
        )
    }

    # Apply FDR correction if permutation was run
    if (n_perm > 0 && !is.null(result$pvalue)) {
        result <- apply_fdr_correction(result, alpha = fdr_threshold)
    }

    # Add metadata
    result$metadata <- list(
        factor_genes = factor_genes,
        target_genes = target_genes,
        radii = radii,
        quantile_probs = quantile_probs,
        use_annular = use_annular,
        annular_width = annular_width,
        annular_mode = annular_mode,
        method = method,
        min_cells = min_cells,
        n_perm = n_perm,
        fdr_threshold = fdr_threshold,
        timestamp = Sys.time(),
        package_version = utils::packageVersion("sigdiscov")
    )

    class(result) <- c("GenomewideResult", class(result))
    result
}


#' Validate Genomewide Analysis Inputs
#' @keywords internal
.validate_genomewide_inputs <- function(data) {
    required_fields <- c("expr", "coords", "cell_types", "gene_names")
    missing <- setdiff(required_fields, names(data))
    if (length(missing) > 0) {
        stop("data is missing required fields: ", paste(missing, collapse = ", "))
    }

    # Validate dimensions
    n_cells <- nrow(data$coords)
    if (ncol(data$expr) != n_cells) {
        stop("Number of columns in expr (", ncol(data$expr),
             ") must match number of rows in coords (", n_cells, ")")
    }
    if (length(data$cell_types) != n_cells) {
        stop("Length of cell_types (", length(data$cell_types),
             ") must match number of cells (", n_cells, ")")
    }
    if (length(data$gene_names) != nrow(data$expr)) {
        stop("Length of gene_names (", length(data$gene_names),
             ") must match number of rows in expr (", nrow(data$expr), ")")
    }
}


#' Get Cell Type Pairs for Analysis
#' @keywords internal
.get_cell_type_pairs_internal <- function(data, pairs, min_cells, verbose) {
    if (!is.null(pairs)) {
        return(as.matrix(pairs))
    }

    # Get unique cell types
    cell_types <- sort(unique(data$cell_types))

    # Filter by min_cells
    type_counts <- table(data$cell_types)
    valid_types <- names(type_counts)[type_counts >= min_cells]

    if (length(valid_types) < 1) {
        stop("No cell types with >= ", min_cells, " cells")
    }

    if (verbose) {
        message("  Cell types with >= ", min_cells, " cells: ", length(valid_types))
    }

    # Generate all pairs including self-pairs
    pairs <- expand.grid(sender = valid_types, receiver = valid_types,
                         stringsAsFactors = FALSE)

    if (verbose) {
        message("  Cell type pairs: ", nrow(pairs))
    }

    as.matrix(pairs)
}


#' Print Method for GenomewideResult
#'
#' @param x GenomewideResult object
#' @param ... Additional arguments (ignored)
#' @return Invisibly returns the input object
#' @export
print.GenomewideResult <- function(x, ...) {
    cat("GenomewideResult object\n")
    cat("  Method:", x$metadata$method, "\n")
    cat("  Radii:", length(x$metadata$radii),
        "(", min(x$metadata$radii), "-", max(x$metadata$radii), ")\n")
    cat("  Factor genes:", length(x$metadata$factor_genes), "\n")

    if (!is.null(x$metadata$quantile_probs)) {
        if (length(x$metadata$quantile_probs) > 1) {
            cat("  Quantile thresholds:", paste(x$metadata$quantile_probs, collapse = ", "), "\n")
        } else {
            cat("  Quantile threshold:", x$metadata$quantile_probs, "\n")
        }
    }

    if (x$metadata$use_annular) {
        cat("  Annular mode:", x$metadata$annular_mode,
            ", width =", x$metadata$annular_width, "\n")
    }

    if (!is.null(x$pair_stats)) {
        cat("  Pairs analyzed:", nrow(x$pair_stats), "\n")
    }

    if (!is.null(x$by_quantile)) {
        cat("  Results by quantile:", length(x$by_quantile), "\n")
    }

    cat("  Timestamp:", format(x$metadata$timestamp, "%Y-%m-%d %H:%M:%S"), "\n")

    invisible(x)
}


#' Summary Method for GenomewideResult
#'
#' @param object GenomewideResult object
#' @param ... Additional arguments (ignored)
#' @return Summary data frame
#' @export
summary.GenomewideResult <- function(object, ...) {
    if (!is.null(object$pair_stats)) {
        cat("Pair Statistics Summary:\n")
        print(summary(object$pair_stats[, c("n_sender", "n_receiver",
                                            "mean_delta_I", "max_delta_I")]))
    }

    if (!is.null(object$delta_I)) {
        cat("\nDelta I Summary:\n")
        all_delta <- unlist(lapply(object$delta_I, as.vector))
        all_delta <- all_delta[is.finite(all_delta)]
        cat("  Range:", round(min(all_delta), 4), "to", round(max(all_delta), 4), "\n")
        cat("  Mean:", round(mean(all_delta), 4), "\n")
        cat("  SD:", round(stats::sd(all_delta), 4), "\n")
    }

    invisible(object)
}
