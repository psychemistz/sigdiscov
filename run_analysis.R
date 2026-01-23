#!/usr/bin/env Rscript
#' =============================================================================
#' SIGDISCOV: Unified Spatial Signature Analysis Script
#' =============================================================================
#'
#' This script provides a complete pipeline for comparing spatial correlation
#' metrics (Moran's I, I_ND) with reference signatures (CytoSig, SecAct).
#'
#' Usage:
#'   Rscript run_analysis.R --vst <vst_file> --output <output_dir> [options]
#'
#' Required arguments:
#'   --vst       Path to VST-transformed expression matrix (genes x spots)
#'   --output    Output directory for results
#'
#' Optional arguments:
#'   --cytosig   Path to CytoSig signature file (default: dataset/db/signature.centroid)
#'   --secact    Path to SecAct signature file (default: dataset/db/AllSigFilteredBy_*.tsv.gz)
#'   --radius    Spatial radius for weight matrix (default: 200)
#'   --sigma     RBF kernel bandwidth (default: 100)
#'   --n_perm    Number of permutations for z-scores (default: 1000, 0 to skip)
#'   --factors   Comma-separated list of factors to analyze (default: all common)
#'   --seed      Random seed for reproducibility (default: 42)
#'
#' =============================================================================

suppressPackageStartupMessages({
    library(data.table)
    library(sigdiscov)
    library(Matrix)
})

# =============================================================================
# CONFIGURATION
# =============================================================================

#' Parse command line arguments or use defaults
parse_args <- function() {
    args <- commandArgs(trailingOnly = TRUE)

    # Default configuration
    config <- list(
        vst_file = NULL,
        output_dir = "output",
        cytosig_file = "dataset/db/signature.centroid",
        secact_file = "dataset/db/AllSigFilteredBy_MoranI_TCGA_ICGC_0.25_ds3.tsv.gz",
        radius = 200,
        sigma = 100,
        n_perm = 1000,
        factors = NULL,  # NULL means all common factors
        seed = 42
    )

    # Parse arguments
    i <- 1
    while (i <= length(args)) {
        if (args[i] == "--vst") {
            config$vst_file <- args[i + 1]
            i <- i + 2
        } else if (args[i] == "--output") {
            config$output_dir <- args[i + 1]
            i <- i + 2
        } else if (args[i] == "--cytosig") {
            config$cytosig_file <- args[i + 1]
            i <- i + 2
        } else if (args[i] == "--secact") {
            config$secact_file <- args[i + 1]
            i <- i + 2
        } else if (args[i] == "--radius") {
            config$radius <- as.numeric(args[i + 1])
            i <- i + 2
        } else if (args[i] == "--sigma") {
            config$sigma <- as.numeric(args[i + 1])
            i <- i + 2
        } else if (args[i] == "--n_perm") {
            config$n_perm <- as.integer(args[i + 1])
            i <- i + 2
        } else if (args[i] == "--factors") {
            config$factors <- strsplit(args[i + 1], ",")[[1]]
            i <- i + 2
        } else if (args[i] == "--seed") {
            config$seed <- as.integer(args[i + 1])
            i <- i + 2
        } else if (args[i] == "--help" || args[i] == "-h") {
            cat("Usage: Rscript run_analysis.R --vst <file> --output <dir> [options]\n")
            cat("\nRun 'Rscript run_analysis.R --help' for more information.\n")
            quit(status = 0)
        } else {
            i <- i + 1
        }
    }

    config
}

# =============================================================================
# DATA LOADING FUNCTIONS
# =============================================================================

#' Load CytoSig reference signatures
load_cytosig <- function(file) {
    cat("Loading CytoSig...\n")
    header <- strsplit(readLines(file, n = 1), "\t")[[1]]
    data <- fread(file, header = FALSE, skip = 1, sep = "\t")
    genes <- data[[1]]
    mat <- as.matrix(data[, -1, with = FALSE])
    rownames(mat) <- genes
    colnames(mat) <- header
    cat("  Genes:", nrow(mat), ", Factors:", ncol(mat), "\n")
    mat
}

#' Load SecAct reference signatures
load_secact <- function(file) {
    cat("Loading SecAct...\n")
    if (grepl("\\.gz$", file)) {
        header <- strsplit(readLines(gzfile(file), n = 1), "\t")[[1]]
    } else {
        header <- strsplit(readLines(file, n = 1), "\t")[[1]]
    }
    data <- fread(file, header = FALSE, skip = 1, sep = "\t")
    genes <- data[[1]]
    mat <- as.matrix(data[, -1, with = FALSE])
    rownames(mat) <- genes
    colnames(mat) <- header
    cat("  Genes:", nrow(mat), ", Factors:", ncol(mat), "\n")
    mat
}

#' Load VST expression data
load_vst <- function(file) {
    cat("Loading VST data...\n")
    header <- strsplit(readLines(file, n = 1), "\t")[[1]]
    data <- fread(file, header = FALSE, skip = 1, sep = "\t")
    genes <- data[[1]]
    mat <- as.matrix(data[, -1, with = FALSE])
    rownames(mat) <- genes
    colnames(mat) <- header
    cat("  Genes:", nrow(mat), ", Spots:", ncol(mat), "\n")
    mat
}

# =============================================================================
# ANALYSIS FUNCTIONS
# =============================================================================

#' Compute pairwise Moran's I matrix
compute_moran_matrix <- function(data, coords, radius, sigma, sparse = TRUE) {
    cat("\nComputing Moran's I matrix...\n")
    t1 <- Sys.time()
    result <- pairwise_moran(data, coords, weight_type = "circular",
                             radius = radius, sigma = sigma, sparse_W = sparse)
    t2 <- Sys.time()
    cat("  Completed in", round(difftime(t2, t1, units = "secs"), 1), "sec\n")
    result$moran
}

#' Compute I_ND matrix (cosine similarity)
compute_ind_matrix <- function(data, W) {
    cat("\nComputing I_ND matrix...\n")
    t1 <- Sys.time()

    # Z-normalize data (population SD)
    data_z <- t(apply(data, 1, function(x) {
        n <- length(x)
        m <- mean(x)
        s <- sqrt(sum((x - m)^2) / n)
        if (s < 1e-10) return(rep(0, n))
        (x - m) / s
    }))

    # Compute spatial lags
    lag_G <- as.matrix(data_z %*% Matrix::t(W))

    # Compute norms
    norm_z <- sqrt(rowSums(data_z^2))
    norm_lag <- sqrt(rowSums(lag_G^2))

    # Compute I_ND (cosine similarity)
    dot_products <- tcrossprod(data_z, lag_G)
    ind_mat <- dot_products / (outer(norm_z, norm_lag))
    ind_mat[is.nan(ind_mat)] <- 0

    rownames(ind_mat) <- colnames(ind_mat) <- rownames(data)

    t2 <- Sys.time()
    cat("  Completed in", round(difftime(t2, t1, units = "secs"), 1), "sec\n")

    list(ind = ind_mat, data_z = data_z)
}

#' Run permutation test for z-scores
run_permutation_test <- function(data_z, W, moran_obs, ind_obs,
                                  gene_subset, n_perm = 1000, seed = 42) {
    cat("\nRunning permutation test (", n_perm, " permutations)...\n", sep = "")

    set.seed(seed)
    t1 <- Sys.time()

    n_genes <- length(gene_subset)
    n_spots <- ncol(data_z)

    # Subset data
    data_sub <- data_z[gene_subset, ]
    moran_sub <- moran_obs[gene_subset, gene_subset]
    ind_sub <- ind_obs[gene_subset, gene_subset]

    # Initialize accumulators
    moran_sum <- matrix(0, n_genes, n_genes)
    moran_sum_sq <- matrix(0, n_genes, n_genes)
    ind_sum <- matrix(0, n_genes, n_genes)
    ind_sum_sq <- matrix(0, n_genes, n_genes)

    W_dense <- as.matrix(W)
    S0 <- sum(W_dense)
    norm_z_sub <- sqrt(rowSums(data_sub^2))

    for (p in seq_len(n_perm)) {
        if (p %% 100 == 0) cat("  Permutation", p, "/", n_perm, "\n")

        # Permute spots
        perm_idx <- sample.int(n_spots)
        data_perm <- data_sub[, perm_idx]

        # Compute permuted Moran's I
        moran_perm <- tcrossprod(data_perm %*% W_dense, data_perm) / S0

        # Compute permuted I_ND
        lag_perm <- data_perm %*% t(W_dense)
        norm_lag_perm <- sqrt(rowSums(lag_perm^2))
        dot_perm <- tcrossprod(data_perm, lag_perm)
        ind_perm <- dot_perm / (outer(norm_z_sub, norm_lag_perm))
        ind_perm[is.nan(ind_perm)] <- 0

        # Update accumulators
        moran_sum <- moran_sum + moran_perm
        moran_sum_sq <- moran_sum_sq + moran_perm^2
        ind_sum <- ind_sum + ind_perm
        ind_sum_sq <- ind_sum_sq + ind_perm^2
    }

    # Compute z-scores
    moran_mean <- moran_sum / n_perm
    moran_var <- moran_sum_sq / n_perm - moran_mean^2
    moran_sd <- sqrt(pmax(moran_var, .Machine$double.eps))
    moran_z <- (moran_sub - moran_mean) / moran_sd

    ind_mean <- ind_sum / n_perm
    ind_var <- ind_sum_sq / n_perm - ind_mean^2
    ind_sd <- sqrt(pmax(ind_var, .Machine$double.eps))
    ind_z <- (ind_sub - ind_mean) / ind_sd

    rownames(moran_z) <- colnames(moran_z) <- gene_subset
    rownames(ind_z) <- colnames(ind_z) <- gene_subset

    t2 <- Sys.time()
    cat("  Completed in", round(difftime(t2, t1, units = "mins"), 1), "min\n")

    list(moran_z = moran_z, ind_z = ind_z)
}

#' Build correlation matrix for a factor
build_correlation <- function(factor_name, cytosig, secact, moran, ind,
                               moran_z, ind_z, common_genes) {
    # Check availability
    if (!factor_name %in% colnames(cytosig)) return(NULL)
    if (!factor_name %in% colnames(secact)) return(NULL)
    if (!factor_name %in% colnames(moran_z)) return(NULL)

    df <- data.frame(
        CytoSig = cytosig[common_genes, factor_name],
        SecAct = secact[common_genes, factor_name],
        Moran = moran[common_genes, factor_name],
        I_ND = ind[common_genes, factor_name],
        Moran_z = moran_z[common_genes, factor_name],
        I_ND_z = ind_z[common_genes, factor_name]
    )
    cor(df, use = "complete.obs")
}

# =============================================================================
# MAIN ANALYSIS PIPELINE
# =============================================================================

run_analysis <- function(config) {
    cat(strrep("=", 70), "\n")
    cat("SIGDISCOV: Spatial Signature Analysis\n")
    cat(strrep("=", 70), "\n\n")

    # Create output directory
    dir.create(config$output_dir, showWarnings = FALSE, recursive = TRUE)

    # Load reference data
    cat("--- Loading Reference Data ---\n")
    cytosig <- load_cytosig(config$cytosig_file)
    secact <- load_secact(config$secact_file)

    # Load VST data
    cat("\n--- Loading Expression Data ---\n")
    vst_data <- load_vst(config$vst_file)
    vst_genes <- rownames(vst_data)

    # Parse coordinates
    spot_names <- colnames(vst_data)
    spot_coords <- parse_spot_names(spot_names)
    coords <- data.frame(
        x = spot_coords$col * 100,
        y = spot_coords$row * 100 * sqrt(3) / 2
    )

    # Find common genes and factors
    common_genes <- Reduce(intersect, list(rownames(cytosig), rownames(secact), vst_genes))
    common_factors <- intersect(colnames(cytosig), colnames(secact))

    # Determine factors to analyze
    if (is.null(config$factors)) {
        target_factors <- intersect(common_factors, vst_genes)
    } else {
        target_factors <- intersect(config$factors, intersect(common_factors, vst_genes))
    }

    cat("\nCommon genes:", length(common_genes), "\n")
    cat("Target factors:", length(target_factors), "\n")

    # Genes for permutation: common_genes + target_factors
    perm_genes <- unique(c(common_genes, target_factors))
    cat("Genes for analysis:", length(perm_genes), "\n")

    # Compute spatial metrics
    cat("\n--- Computing Spatial Metrics ---\n")

    # Create weight matrix
    W_result <- create_weights_visium(coords, radius = config$radius,
                                       sigma = config$sigma, sparse = TRUE)
    W <- W_result$W

    # Compute Moran's I
    moran_mat <- compute_moran_matrix(vst_data, coords, config$radius, config$sigma)

    # Compute I_ND
    ind_result <- compute_ind_matrix(vst_data, W)
    ind_mat <- ind_result$ind
    data_z <- ind_result$data_z

    # Permutation test for z-scores
    if (config$n_perm > 0) {
        perm_result <- run_permutation_test(
            data_z, W, moran_mat, ind_mat,
            perm_genes, config$n_perm, config$seed
        )
        moran_z <- perm_result$moran_z
        ind_z <- perm_result$ind_z
    } else {
        cat("\nSkipping permutation test (n_perm = 0)\n")
        moran_z <- moran_mat[perm_genes, perm_genes]
        ind_z <- ind_mat[perm_genes, perm_genes]
    }

    # Build correlation matrices
    cat("\n--- Building Correlation Matrices ---\n")

    results <- list()
    cor_list <- list()

    for (factor in target_factors) {
        cor_mat <- build_correlation(factor, cytosig, secact, moran_mat, ind_mat,
                                     moran_z, ind_z, common_genes)
        if (!is.null(cor_mat)) {
            cor_list[[factor]] <- cor_mat
        }
    }

    # Print key results
    if ("IFNG" %in% names(cor_list)) {
        cat("\n=== IFNG Correlation (6x6) ===\n")
        print(round(cor_list[["IFNG"]], 4))
        results$ifng <- cor_list[["IFNG"]]
    }

    if ("TGFB1" %in% names(cor_list)) {
        cat("\n=== TGFB1 Correlation (6x6) ===\n")
        print(round(cor_list[["TGFB1"]], 4))
        results$tgfb1 <- cor_list[["TGFB1"]]
    }

    # Mean correlation
    if (length(cor_list) > 0) {
        mean_cor <- Reduce(`+`, cor_list) / length(cor_list)
        cat("\n=== Mean Correlation (", length(cor_list), " factors) ===\n", sep = "")
        print(round(mean_cor, 4))
        results$mean <- mean_cor
    }

    # Save results
    results$per_factor <- cor_list
    results$factors <- names(cor_list)
    results$config <- config
    results$moran <- moran_mat
    results$ind <- ind_mat
    results$moran_z <- moran_z
    results$ind_z <- ind_z

    output_file <- file.path(config$output_dir, "analysis_results.rds")
    saveRDS(results, output_file)
    cat("\nResults saved to:", output_file, "\n")

    # Summary
    cat("\n", strrep("=", 70), "\n", sep = "")
    cat("SUMMARY\n")
    cat(strrep("=", 70), "\n")
    cat("\nKey correlations (Moran vs CytoSig):\n")
    if (!is.null(results$ifng)) cat("  IFNG:", round(results$ifng["Moran", "CytoSig"], 4), "\n")
    if (!is.null(results$tgfb1)) cat("  TGFB1:", round(results$tgfb1["Moran", "CytoSig"], 4), "\n")
    if (!is.null(results$mean)) cat("  Mean:", round(results$mean["Moran", "CytoSig"], 4), "\n")

    cat("\nKey correlations (Moran vs SecAct):\n")
    if (!is.null(results$ifng)) cat("  IFNG:", round(results$ifng["Moran", "SecAct"], 4), "\n")
    if (!is.null(results$tgfb1)) cat("  TGFB1:", round(results$tgfb1["Moran", "SecAct"], 4), "\n")
    if (!is.null(results$mean)) cat("  Mean:", round(results$mean["Moran", "SecAct"], 4), "\n")

    cat("\nDone!\n")

    invisible(results)
}

# =============================================================================
# ENTRY POINT
# =============================================================================

# Run if executed directly
if (!interactive()) {
    config <- parse_args()

    # Check required arguments
    if (is.null(config$vst_file)) {
        # Use default for testing
        config$vst_file <- "dataset/1_vst.tsv"
    }

    if (!file.exists(config$vst_file)) {
        stop("VST file not found: ", config$vst_file)
    }

    run_analysis(config)
}
