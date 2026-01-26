#!/usr/bin/env Rscript
# COAD TGFB1 Macrophage -> CAF I_ND Analysis using h5ad VST data
# Compare with output/results_macrophage_to_CAF_TGFB1_annular_20250818_150423.csv

library(sigdiscov)
library(Matrix)
library(rhdf5)

cat("=== COAD TGFB1 Analysis (h5ad VST) ===\n\n")

# =============================================================================
# 1. Load COAD Data from h5ad
# =============================================================================

cat("--- Loading COAD h5ad Data ---\n")

h5_file <- "dataset/cosmx/coad_spatial.h5ad"

# Load cell IDs
cell_ids <- h5read(h5_file, "obs/cell_ID")
n_cells <- length(cell_ids)
cat(sprintf("Cells: %d\n", n_cells))

# Load gene names
gene_names <- h5read(h5_file, "var/_index")
n_genes <- length(gene_names)
cat(sprintf("Genes: %d\n", n_genes))

# Load cell types (categorical)
cell_type_categories <- h5read(h5_file, "obs/cell_type/categories")
cell_type_codes <- h5read(h5_file, "obs/cell_type/codes")
cell_types <- cell_type_categories[cell_type_codes + 1]  # 0-indexed
cat(sprintf("Cell types: %d unique\n", length(unique(cell_types))))

# Load spatial coordinates (already in μm)
spatial <- h5read(h5_file, "obsm/spatial")
# spatial is 2 x n_cells, transpose to n_cells x 2
coords_um <- t(spatial)
colnames(coords_um) <- c("x", "y")
cat(sprintf("Coordinates (μm): X [%.1f, %.1f], Y [%.1f, %.1f]\n",
    min(coords_um[,1]), max(coords_um[,1]), min(coords_um[,2]), max(coords_um[,2])))

# Load VST expression matrix (sparse CSR format)
cat("Loading VST expression matrix (sparse)...\n")
x_data <- h5read(h5_file, "X/data")
x_indices <- h5read(h5_file, "X/indices")
x_indptr <- h5read(h5_file, "X/indptr")

# Convert to dgCMatrix (CSC format) - need to transpose CSR to get cells x genes
# CSR: rows are cells, cols are genes
# indptr has length n_cells + 1
# Create as dgRMatrix first, then transpose
cat("Converting to sparse matrix...\n")

# CSR format: indptr is row pointers, indices are column indices
# dgRMatrix is R's CSR format
expr_csr <- new("dgRMatrix",
    p = as.integer(x_indptr),
    j = as.integer(x_indices),
    x = as.numeric(x_data),
    Dim = c(n_cells, n_genes))

# Convert to dgCMatrix (CSC) by transposing twice or direct conversion
# For our purposes, we want genes x cells, so transpose once
expr <- t(expr_csr)  # Now genes x cells (dgCMatrix)

cat(sprintf("Expression matrix: %d genes x %d cells\n", nrow(expr), ncol(expr)))
cat(sprintf("Non-zero entries: %d (%.1f%% dense)\n",
    length(expr@x), 100 * length(expr@x) / prod(dim(expr))))

H5close()

# =============================================================================
# 2. Define Sender and Receiver
# =============================================================================

cat("\n--- Defining Sender/Receiver ---\n")

cat("Cell type distribution:\n")
print(table(cell_types))

# Define sender (macrophage) and receiver (CAF)
sender_type <- "macrophage"
receiver_type <- "CAF"

sender_mask <- cell_types == sender_type
receiver_mask <- cell_types == receiver_type

n_senders <- sum(sender_mask)
n_receivers <- sum(receiver_mask)

cat(sprintf("\nSender (%s): %d cells\n", sender_type, n_senders))
cat(sprintf("Receiver (%s): %d cells\n", receiver_type, n_receivers))

sender_idx <- which(sender_mask)
receiver_idx <- which(receiver_mask)

# =============================================================================
# 3. Prepare Data for Analysis
# =============================================================================

cat("\n--- Preparing Data ---\n")

# Check if TGFB1 exists
factor_gene <- "TGFB1"
if (!factor_gene %in% gene_names) {
    stop("TGFB1 not found in genes")
}
tgfb1_idx <- which(gene_names == factor_gene)
cat(sprintf("TGFB1 index: %d\n", tgfb1_idx))

# Set row names
rownames(expr) <- gene_names

# Z-normalize expression (gene-wise) - convert to dense for specific subsets
cat("Extracting sender/receiver expression...\n")

sender_coords <- coords_um[sender_idx, ]
receiver_coords <- coords_um[receiver_idx, ]

# Factor expression in senders
factor_expr_raw <- as.vector(expr[factor_gene, sender_idx])
factor_expr <- (factor_expr_raw - mean(factor_expr_raw)) / sd(factor_expr_raw)

cat(sprintf("TGFB1 in senders: mean=%.3f, sd=%.3f\n",
    mean(factor_expr_raw), sd(factor_expr_raw)))

# Receiver expression - keep sparse, normalize later
receiver_expr <- expr[, receiver_idx]
cat(sprintf("Receiver expression: %d genes x %d cells\n",
    nrow(receiver_expr), ncol(receiver_expr)))

# Z-normalize receiver expression (gene-wise)
cat("Z-normalizing receiver expression...\n")
# For sparse matrix, compute mean and sd per gene
gene_means <- Matrix::rowMeans(receiver_expr)
gene_sds <- sqrt(Matrix::rowMeans(receiver_expr^2) - gene_means^2)
gene_sds[gene_sds < 1e-10] <- 1  # Avoid division by zero

# Normalize: (x - mean) / sd
# This is tricky for sparse matrices - convert subset at a time
# For efficiency, we'll normalize during the lag computation

# =============================================================================
# 4. Compute I_ND at Multiple Radii (Annular/Ring)
# =============================================================================

cat("\n--- Computing I_ND with Annular Neighborhoods ---\n")

# Radii matching previous analysis
radii <- c(10, 20, 30, 50, 100, 200, 300, 500)
n_radii <- length(radii)

# Results storage
I_ND_matrix <- matrix(NA, nrow = n_genes, ncol = n_radii)
rownames(I_ND_matrix) <- gene_names
colnames(I_ND_matrix) <- paste0("r", radii)

cat(sprintf("Radii: %s μm\n", paste(radii, collapse = ", ")))
cat("Computing annular I_ND...\n\n")

for (r in seq_along(radii)) {
    radius <- radii[r]
    inner_radius <- if (r == 1) 0 else radii[r - 1]
    sigma <- radius / 3  # Gaussian sigma

    cat(sprintf("  Radius %d μm (ring %d-%d μm)...", radius, inner_radius, radius))

    # Create weight matrix with annular ring
    W <- create_weights_sc(
        sender_coords, receiver_coords,
        radius = radius,
        inner_radius = inner_radius,
        sigma = sigma
    )

    # Count connections
    mean_conn <- length(W@x) / nrow(sender_coords)

    if (length(W@x) < 10) {
        cat(sprintf(" skipped (%.1f conn/sender)\n", mean_conn))
        next
    }

    # Compute spatial lag for all genes
    # lag_G[s, g] = sum_r W[s,r] * ((expr[g,r] - mean_g) / sd_g)
    # = (1/sd_g) * (sum_r W[s,r] * expr[g,r] - mean_g * sum_r W[s,r])
    # Since W is row-normalized, sum_r W[s,r] = 1
    # = (1/sd_g) * (W %*% expr^T - mean_g)

    # Compute raw spatial lag
    lag_G_raw <- as.matrix(W %*% t(receiver_expr))  # n_senders x n_genes

    # Z-normalize using gene means and sds
    lag_G <- sweep(lag_G_raw, 2, gene_means, "-")
    lag_G <- sweep(lag_G, 2, gene_sds, "/")

    # Compute I_ND for each gene
    for (g in seq_len(n_genes)) {
        lag_g <- lag_G[, g]
        norm_lag <- sqrt(sum(lag_g^2))
        if (norm_lag < 1e-10 || is.na(norm_lag)) {
            I_ND_matrix[g, r] <- NA
        } else {
            I_ND_matrix[g, r] <- sum(factor_expr * lag_g) /
                (sqrt(sum(factor_expr^2)) * norm_lag)
        }
    }

    cat(sprintf(" done (%.1f conn/sender)\n", mean_conn))
}

# =============================================================================
# 5. Load Previous Results for Comparison
# =============================================================================

cat("\n--- Loading Previous Results ---\n")

prev_results <- read.csv("output/results_macrophage_to_CAF_TGFB1_annular_20250818_150423.csv")
cat(sprintf("Previous results: %d rows\n", nrow(prev_results)))

# Filter to quantile 0 (all senders)
prev_q0 <- prev_results[prev_results$min_expr_quantile == 0, ]
cat(sprintf("Previous results (quantile=0): %d rows\n", nrow(prev_q0)))

# =============================================================================
# 6. Compare Results
# =============================================================================

cat("\n--- Comparing Results ---\n\n")

for (r in seq_along(radii)) {
    radius <- radii[r]
    inner_radius <- if (r == 1) 0 else radii[r - 1]

    prev_r <- prev_q0[prev_q0$radius == radius, ]

    if (nrow(prev_r) == 0) next

    common_genes <- intersect(gene_names, prev_r$target_gene)
    new_values <- I_ND_matrix[common_genes, r]
    prev_idx <- match(common_genes, prev_r$target_gene)
    prev_values <- prev_r$morans_i_normalized[prev_idx]

    valid <- !is.na(new_values) & !is.na(prev_values)

    if (sum(valid) < 50) {
        cat(sprintf("Radius %d: Only %d valid pairs\n", radius, sum(valid)))
        next
    }

    cor_pearson <- cor(new_values[valid], prev_values[valid], method = "pearson")
    cor_spearman <- cor(new_values[valid], prev_values[valid], method = "spearman")
    ratio <- mean(abs(prev_values[valid])) / mean(abs(new_values[valid]))

    cat(sprintf("Radius %d μm (ring %d-%d):\n", radius, inner_radius, radius))
    cat(sprintf("  Common genes: %d\n", sum(valid)))
    cat(sprintf("  Pearson r:    %.4f\n", cor_pearson))
    cat(sprintf("  Spearman ρ:   %.4f\n", cor_spearman))
    cat(sprintf("  Scale ratio:  %.4f\n", ratio))
    cat("\n")
}

# =============================================================================
# 7. Top Genes Comparison
# =============================================================================

cat("\n--- Top Genes Comparison (r=10 μm) ---\n\n")

new_r10 <- I_ND_matrix[, "r10"]
prev_r10 <- prev_q0[prev_q0$radius == 10, ]

common_genes <- intersect(names(new_r10[!is.na(new_r10)]), prev_r10$target_gene)
new_r10_common <- new_r10[common_genes]
prev_idx <- match(common_genes, prev_r10$target_gene)
prev_r10_common <- prev_r10$morans_i_normalized[prev_idx]
names(prev_r10_common) <- common_genes

cat("Top 20 by NEW I_ND (r=10):\n")
top_new <- head(sort(new_r10_common, decreasing = TRUE), 20)
for (i in seq_along(top_new)) {
    gene <- names(top_new)[i]
    cat(sprintf("  %2d. %-10s  new=%.4f  prev=%.4f\n",
        i, gene, top_new[i], prev_r10_common[gene]))
}

cat("\nTop 20 by PREVIOUS I_ND (r=10):\n")
top_prev <- head(sort(prev_r10_common, decreasing = TRUE), 20)
for (i in seq_along(top_prev)) {
    gene <- names(top_prev)[i]
    cat(sprintf("  %2d. %-10s  prev=%.4f  new=%.4f\n",
        i, gene, top_prev[i], new_r10_common[gene]))
}

# =============================================================================
# 8. Save Results
# =============================================================================

cat("\n--- Saving Results ---\n")

comparison_df <- data.frame(
    gene = common_genes,
    I_ND_new = new_r10_common,
    I_ND_prev = prev_r10_common,
    stringsAsFactors = FALSE
)
comparison_df$diff <- comparison_df$I_ND_new - comparison_df$I_ND_prev
comparison_df <- comparison_df[order(-comparison_df$I_ND_new), ]

write.csv(comparison_df, "output/coad_tgfb1_h5ad_comparison.csv", row.names = FALSE)
cat("Comparison saved to: output/coad_tgfb1_h5ad_comparison.csv\n")

write.csv(I_ND_matrix, "output/coad_tgfb1_h5ad_IND_matrix.csv")
cat("I_ND matrix saved to: output/coad_tgfb1_h5ad_IND_matrix.csv\n")

cat("\n=== Analysis Complete ===\n")
