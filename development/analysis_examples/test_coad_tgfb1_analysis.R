#!/usr/bin/env Rscript
# COAD TGFB1 Macrophage -> CAF I_ND Analysis
# Compare with output/results_macrophage_to_CAF_TGFB1_annular_20250818_150423.csv

library(sigdiscov)
library(Matrix)

cat("=== COAD TGFB1 Analysis: Macrophage -> CAF ===\n\n")

# =============================================================================
# 1. Load COAD Data from RData
# =============================================================================

cat("--- Loading COAD Data ---\n")

load("dataset/cosmx/COAD.RData")

cat(sprintf("Expression (raw): %d cells x %d genes\n", nrow(raw), ncol(raw)))
cat(sprintf("Coordinates (xy): %d cells, units: mm\n", nrow(xy)))
cat(sprintf("Cell types: %d unique types\n", length(unique(clust))))

# Convert coordinates from mm to μm (same as h5ad)
coords_um <- xy * 1000
cat(sprintf("Coordinates converted to μm: X [%.1f, %.1f], Y [%.1f, %.1f]\n",
    min(coords_um[,1]), max(coords_um[,1]), min(coords_um[,2]), max(coords_um[,2])))

# =============================================================================
# 2. Define Sender and Receiver
# =============================================================================

cat("\n--- Defining Sender/Receiver ---\n")

# Cell types in RData
cat("Cell type distribution:\n")
print(table(clust))

# Define sender (macrophage) and receiver (CAF)
sender_type <- "macrophage"
receiver_type <- "CAF"

sender_mask <- clust == sender_type
receiver_mask <- clust == receiver_type

n_senders <- sum(sender_mask)
n_receivers <- sum(receiver_mask)

cat(sprintf("\nSender (%s): %d cells\n", sender_type, n_senders))
cat(sprintf("Receiver (%s): %d cells\n", receiver_type, n_receivers))

# Get indices
sender_idx <- which(sender_mask)
receiver_idx <- which(receiver_mask)

# =============================================================================
# 3. Prepare Data for Analysis
# =============================================================================

cat("\n--- Preparing Data ---\n")

# Expression matrix: raw is cells x genes, need genes x cells
# The raw matrix in RData is cells x genes
# For sigdiscov, we need genes x cells
expr <- t(raw)  # Now genes x cells
cat(sprintf("Expression matrix: %d genes x %d cells\n", nrow(expr), ncol(expr)))

# Check if TGFB1 exists
factor_gene <- "TGFB1"
if (!factor_gene %in% rownames(expr)) {
    stop("TGFB1 not found in genes")
}
cat(sprintf("TGFB1 found in %d/%d cells\n",
    sum(expr[factor_gene,] > 0), ncol(expr)))

# Z-normalize expression (gene-wise)
cat("Z-normalizing expression...\n")
expr_norm <- t(scale(t(expr)))
expr_norm[is.na(expr_norm)] <- 0

# Extract sender and receiver data
sender_coords <- coords_um[sender_idx, ]
receiver_coords <- coords_um[receiver_idx, ]

factor_expr_raw <- expr[factor_gene, sender_idx]
receiver_expr_norm <- expr_norm[, receiver_idx]

# Z-normalize factor expression
factor_expr <- (factor_expr_raw - mean(factor_expr_raw)) / sd(factor_expr_raw)

cat(sprintf("Sender coords: %d cells\n", nrow(sender_coords)))
cat(sprintf("Receiver coords: %d cells\n", nrow(receiver_coords)))

# =============================================================================
# 4. Compute I_ND at Multiple Radii (Annular/Ring)
# =============================================================================

cat("\n--- Computing I_ND with Annular Neighborhoods ---\n")

# Radii matching previous analysis
radii <- c(10, 20, 30, 50, 100, 200, 300, 500)

# Results storage
n_genes <- nrow(expr)
gene_names <- rownames(expr)
n_radii <- length(radii)

# Storage for I_ND at each radius
I_ND_matrix <- matrix(NA, nrow = n_genes, ncol = n_radii)
rownames(I_ND_matrix) <- gene_names
colnames(I_ND_matrix) <- paste0("r", radii)

# Connection counts
conn_counts <- numeric(n_radii)

cat(sprintf("Radii: %s μm\n", paste(radii, collapse = ", ")))
cat("Computing annular I_ND...\n\n")

for (r in seq_along(radii)) {
    radius <- radii[r]
    inner_radius <- if (r == 1) 0 else radii[r - 1]
    sigma <- (radius - inner_radius) / 3  # Gaussian sigma for the ring

    cat(sprintf("  Radius %d μm (ring %d-%d μm)...", radius, inner_radius, radius))

    # Create weight matrix with annular ring
    W <- create_weights_sc(
        sender_coords, receiver_coords,
        radius = radius,
        inner_radius = inner_radius,
        sigma = sigma
    )

    # Count connections
    conn_counts[r] <- length(W@x)
    mean_conn <- length(W@x) / nrow(sender_coords)

    if (length(W@x) < 10) {
        cat(sprintf(" skipped (%.0f connections)\n", mean_conn))
        next
    }

    # Compute spatial lag for all genes
    lag_G <- as.matrix(W %*% t(receiver_expr_norm))  # n_senders x n_genes

    # Compute I_ND for each gene
    for (g in seq_len(n_genes)) {
        lag_g <- lag_G[, g]
        norm_lag <- sqrt(sum(lag_g^2))
        if (norm_lag < 1e-10) {
            I_ND_matrix[g, r] <- NA
        } else {
            I_ND_matrix[g, r] <- sum(factor_expr * lag_g) /
                (sqrt(sum(factor_expr^2)) * norm_lag)
        }
    }

    cat(sprintf(" done (%.0f conn/sender)\n", mean_conn))
}

# =============================================================================
# 5. Load Previous Results for Comparison
# =============================================================================

cat("\n--- Loading Previous Results ---\n")

prev_results <- read.csv("output/results_macrophage_to_CAF_TGFB1_annular_20250818_150423.csv")
cat(sprintf("Previous results: %d rows\n", nrow(prev_results)))
cat(sprintf("Columns: %s\n", paste(head(names(prev_results), 10), collapse = ", ")))

# Filter to quantile 0 (all senders)
prev_q0 <- prev_results[prev_results$min_expr_quantile == 0, ]
cat(sprintf("Previous results (quantile=0): %d rows\n", nrow(prev_q0)))

# =============================================================================
# 6. Compare Results
# =============================================================================

cat("\n--- Comparing Results ---\n\n")

# For each radius, compare I_ND values
for (r in seq_along(radii)) {
    radius <- radii[r]
    inner_radius <- if (r == 1) 0 else radii[r - 1]

    # Get previous results for this radius
    prev_r <- prev_q0[prev_q0$radius == radius, ]

    if (nrow(prev_r) == 0) {
        cat(sprintf("Radius %d: No previous results\n", radius))
        next
    }

    # Match genes
    common_genes <- intersect(gene_names, prev_r$target_gene)

    if (length(common_genes) < 100) {
        cat(sprintf("Radius %d: Only %d common genes\n", radius, length(common_genes)))
        next
    }

    # Get values
    new_values <- I_ND_matrix[common_genes, r]
    prev_idx <- match(common_genes, prev_r$target_gene)
    prev_values <- prev_r$morans_i_normalized[prev_idx]

    # Filter valid pairs
    valid <- !is.na(new_values) & !is.na(prev_values)

    if (sum(valid) < 50) {
        cat(sprintf("Radius %d: Only %d valid pairs\n", radius, sum(valid)))
        next
    }

    # Correlation
    cor_pearson <- cor(new_values[valid], prev_values[valid], method = "pearson")
    cor_spearman <- cor(new_values[valid], prev_values[valid], method = "spearman")

    # Scale comparison (previous may be scaled differently)
    ratio <- mean(abs(prev_values[valid])) / mean(abs(new_values[valid]))

    cat(sprintf("Radius %d μm (ring %d-%d):\n", radius, inner_radius, radius))
    cat(sprintf("  Common genes: %d\n", sum(valid)))
    cat(sprintf("  Pearson r:    %.4f\n", cor_pearson))
    cat(sprintf("  Spearman ρ:   %.4f\n", cor_spearman))
    cat(sprintf("  Scale ratio:  %.4f (prev/new)\n", ratio))
    cat("\n")
}

# =============================================================================
# 7. Top Genes Comparison
# =============================================================================

cat("\n--- Top Genes Comparison (r=10 μm) ---\n\n")

# Get r=10 results
new_r10 <- I_ND_matrix[, "r10"]
prev_r10 <- prev_q0[prev_q0$radius == 10, ]

# Match and sort
common_genes <- intersect(names(new_r10[!is.na(new_r10)]), prev_r10$target_gene)
new_r10_common <- new_r10[common_genes]
prev_idx <- match(common_genes, prev_r10$target_gene)
prev_r10_common <- prev_r10$morans_i_normalized[prev_idx]
names(prev_r10_common) <- common_genes

# Top 20 by new I_ND
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
# 8. Save Comparison Results
# =============================================================================

cat("\n--- Saving Results ---\n")

# Create comparison data frame
comparison_df <- data.frame(
    gene = common_genes,
    I_ND_new = new_r10_common,
    I_ND_prev = prev_r10_common,
    stringsAsFactors = FALSE
)
comparison_df$diff <- comparison_df$I_ND_new - comparison_df$I_ND_prev
comparison_df <- comparison_df[order(-comparison_df$I_ND_new), ]

write.csv(comparison_df, "output/coad_tgfb1_comparison.csv", row.names = FALSE)
cat("Comparison saved to: output/coad_tgfb1_comparison.csv\n")

# Save I_ND matrix
write.csv(I_ND_matrix, "output/coad_tgfb1_IND_matrix.csv")
cat("I_ND matrix saved to: output/coad_tgfb1_IND_matrix.csv\n")

cat("\n=== Analysis Complete ===\n")
