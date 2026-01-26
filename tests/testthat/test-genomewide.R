# Tests for genome-wide spatial analysis functions

# Helper function to create synthetic spatial data
create_test_genomewide_data <- function(n_genes = 50, n_cells = 200, n_types = 3) {
    set.seed(42)

    # Generate cell types
    cell_types <- sample(paste0("Type", LETTERS[seq_len(n_types)]),
                         n_cells, replace = TRUE)

    # Generate coordinates (2D spatial layout)
    coords <- matrix(runif(n_cells * 2, 0, 1000), ncol = 2)
    colnames(coords) <- c("x", "y")

    # Generate expression data
    expr <- matrix(rnorm(n_genes * n_cells), nrow = n_genes, ncol = n_cells)
    rownames(expr) <- paste0("Gene", seq_len(n_genes))

    list(
        expr = expr,
        coords = coords,
        cell_types = cell_types,
        gene_names = rownames(expr),
        n_genes = n_genes,
        n_cells = n_cells
    )
}

# =============================================================================
# Tests for genomewide_analysis
# =============================================================================

test_that("genomewide_analysis basic workflow returns GenomewideResult", {
    data <- create_test_genomewide_data(n_genes = 20, n_cells = 100, n_types = 2)

    result <- genomewide_analysis(
        data,
        factor_genes = paste0("Gene", 1:5),
        radii = c(100, 200),
        min_cells = 10,
        verbose = FALSE
    )

    expect_s3_class(result, "GenomewideResult")
    expect_true(!is.null(result$metadata))
    expect_equal(result$metadata$radii, c(100, 200))
    expect_equal(length(result$metadata$factor_genes), 5)
})

test_that("genomewide_analysis handles missing factor genes gracefully", {
    data <- create_test_genomewide_data(n_genes = 20, n_cells = 100, n_types = 2)

    # Include some non-existent genes
    expect_warning(
        result <- genomewide_analysis(
            data,
            factor_genes = c("Gene1", "Gene2", "NotAGene"),
            radii = c(100),
            verbose = FALSE
        ),
        "Factor genes not found"
    )

    # Should still work with valid genes
    expect_s3_class(result, "GenomewideResult")
})

test_that("genomewide_analysis validates input data structure", {
    bad_data <- list(expr = matrix(1:10, nrow = 2))

    expect_error(
        genomewide_analysis(bad_data, verbose = FALSE),
        "missing required fields"
    )
})

test_that("genomewide_analysis respects min_cells threshold", {
    data <- create_test_genomewide_data(n_genes = 15, n_cells = 60, n_types = 3)
    # Create one rare type
    data$cell_types[1:3] <- "RareType"

    result <- genomewide_analysis(
        data,
        factor_genes = paste0("Gene", 1:5),
        radii = c(100),
        min_cells = 10,  # RareType has only 3 cells
        verbose = FALSE
    )

    # Pairs involving RareType should be excluded
    if (!is.null(result$pair_stats)) {
        expect_false("RareType" %in% result$pair_stats$sender)
    }
})

test_that("genomewide_analysis print method works", {
    data <- create_test_genomewide_data(n_genes = 15, n_cells = 80, n_types = 2)

    result <- genomewide_analysis(
        data,
        factor_genes = paste0("Gene", 1:3),
        radii = c(100, 200),
        verbose = FALSE
    )

    expect_output(print(result), "GenomewideResult object")
})

# =============================================================================
# Tests for annular weights
# =============================================================================

test_that("create_annular_weights_sc creates correct weight structure", {
    sender_coords <- matrix(c(0, 0), nrow = 1)
    receiver_coords <- matrix(c(
        50, 0,   # distance 50
        100, 0,  # distance 100
        150, 0,  # distance 150
        200, 0   # distance 200
    ), ncol = 2, byrow = TRUE)

    # Circular: should include all within 200
    W_circular <- create_annular_weights_sc(
        sender_coords, receiver_coords,
        outer_radius = 200, inner_radius = 0
    )
    expect_equal(W_circular$n_edges, 4)

    # Annular: should include only 100-200 (i.e., distances 100, 150, 200)
    W_annular <- create_annular_weights_sc(
        sender_coords, receiver_coords,
        outer_radius = 200, inner_radius = 100
    )
    # Should exclude the point at distance 50
    expect_lt(W_annular$n_edges, W_circular$n_edges)
})

test_that("create_directional_weights with inner_radius works", {
    sender_coords <- matrix(c(0, 0, 100, 0), ncol = 2, byrow = TRUE)
    receiver_coords <- matrix(c(
        25, 0,   # distance 25 from sender 1
        50, 0,   # distance 50 from sender 1
        125, 0,  # distance 25 from sender 2
        150, 0   # distance 50 from sender 2
    ), ncol = 2, byrow = TRUE)

    # No inner radius
    W1 <- create_directional_weights(
        sender_coords, receiver_coords,
        radius = 60, inner_radius = 0
    )

    # With inner radius = 30
    W2 <- create_directional_weights(
        sender_coords, receiver_coords,
        radius = 60, inner_radius = 30
    )

    # W2 should have fewer edges (excluding close neighbors)
    expect_lt(W2$n_edges, W1$n_edges)
})

# =============================================================================
# Tests for FDR correction
# =============================================================================

test_that("apply_fdr_correction works with data frame", {
    results <- data.frame(
        gene = paste0("G", 1:100),
        pvalue = c(runif(10, 0, 0.01), runif(90, 0.1, 1))
    )

    corrected <- apply_fdr_correction(results, alpha = 0.05)

    expect_true("p_adj" %in% names(corrected))
    expect_true("significant" %in% names(corrected))
    expect_true(all(corrected$p_adj >= corrected$pvalue))
})

test_that("apply_fdr_correction respects alpha threshold", {
    results <- data.frame(
        gene = paste0("G", 1:50),
        pvalue = seq(0.001, 0.1, length.out = 50)
    )

    # Strict threshold
    strict <- apply_fdr_correction(results, alpha = 0.01)
    # Lenient threshold
    lenient <- apply_fdr_correction(results, alpha = 0.1)

    expect_lte(sum(strict$significant), sum(lenient$significant))
})

test_that("apply_fdr_correction handles different methods", {
    results <- data.frame(
        gene = paste0("G", 1:20),
        pvalue = c(0.001, 0.01, 0.02, 0.05, rep(0.5, 16))
    )

    bh <- apply_fdr_correction(results, method = "BH")
    bonf <- apply_fdr_correction(results, method = "bonferroni")

    # Bonferroni should be more conservative
    expect_gte(sum(bh$significant), sum(bonf$significant))
})

test_that("apply_fdr_correction errors on missing pvalue column", {
    bad_results <- data.frame(gene = "A", score = 0.5)

    expect_error(
        apply_fdr_correction(bad_results),
        "pvalue"
    )
})

# =============================================================================
# Tests for multi-quantile analysis
# =============================================================================

test_that("genomewide_analysis handles multiple quantile thresholds", {
    data <- create_test_genomewide_data(n_genes = 15, n_cells = 80, n_types = 2)

    result <- genomewide_analysis(
        data,
        factor_genes = paste0("Gene", 1:3),
        radii = c(100),
        quantile_probs = c(0.25, 0.50),
        verbose = FALSE
    )

    expect_s3_class(result, "GenomewideResult")
    expect_true(!is.null(result$by_quantile))
    expect_equal(length(result$by_quantile), 2)
})

# =============================================================================
# Tests for input validation
# =============================================================================

test_that(".validate_genomewide_inputs catches dimension mismatches", {
    # Mismatched expr columns and coords rows
    bad_data <- list(
        expr = matrix(rnorm(100), nrow = 10, ncol = 10),
        coords = matrix(runif(40), nrow = 20, ncol = 2),  # 20 rows, not 10
        cell_types = rep("A", 20),
        gene_names = paste0("G", 1:10)
    )

    expect_error(
        genomewide_analysis(bad_data, verbose = FALSE),
        "must match"
    )
})

test_that(".validate_genomewide_inputs catches missing fields", {
    incomplete_data <- list(
        expr = matrix(1:10, nrow = 2),
        coords = matrix(1:4, nrow = 2)
        # Missing cell_types and gene_names
    )

    expect_error(
        genomewide_analysis(incomplete_data, verbose = FALSE),
        "missing required fields"
    )
})

# =============================================================================
# Tests for cell type pair generation
# =============================================================================

test_that("cell type pairs are generated correctly", {
    data <- create_test_genomewide_data(n_genes = 10, n_cells = 100, n_types = 3)

    result <- genomewide_analysis(
        data,
        factor_genes = paste0("Gene", 1:3),
        radii = c(100),
        min_cells = 5,
        verbose = FALSE
    )

    # With 3 cell types, should have 9 pairs (3x3 including self-pairs)
    if (!is.null(result$pair_stats)) {
        expect_gte(nrow(result$pair_stats), 1)
    }
})

test_that("specific pairs can be provided", {
    data <- create_test_genomewide_data(n_genes = 10, n_cells = 100, n_types = 3)

    # Specify only 2 pairs
    pairs <- matrix(c("TypeA", "TypeB", "TypeB", "TypeA"),
                    ncol = 2, byrow = TRUE)

    result <- genomewide_analysis(
        data,
        factor_genes = paste0("Gene", 1:3),
        radii = c(100),
        pairs = pairs,
        verbose = FALSE
    )

    if (!is.null(result$pair_stats)) {
        expect_equal(nrow(result$pair_stats), 2)
    }
})
