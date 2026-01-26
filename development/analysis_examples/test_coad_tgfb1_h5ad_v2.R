#!/usr/bin/env Rscript
# COAD TGFB1 Analysis with GLOBAL normalization (matching Python implementation)

library(sigdiscov)
library(Matrix)
library(rhdf5)

cat("=== COAD TGFB1 Analysis (h5ad VST, Global Norm) ===\n\n")

# =============================================================================
# 1. Load COAD Data from h5ad
# =============================================================================

cat("--- Loading COAD h5ad Data ---\n")

h5_file <- "dataset/cosmx/coad_spatial.h5ad"

cell_ids <- h5read(h5_file, "obs/cell_ID")
n_cells <- length(cell_ids)
gene_names <- h5read(h5_file, "var/_index")
n_genes <- length(gene_names)

cell_type_categories <- h5read(h5_file, "obs/cell_type/categories")
cell_type_codes <- h5read(h5_file, "obs/cell_type/codes")
cell_types <- cell_type_categories[cell_type_codes + 1]

spatial <- h5read(h5_file, "obsm/spatial")
coords_um <- t(spatial)

cat(sprintf("Cells: %d, Genes: %d\n", n_cells, n_genes))

# Load VST expression
cat("Loading VST expression...\n")
x_data <- h5read(h5_file, "X/data")
x_indices <- h5read(h5_file, "X/indices")
x_indptr <- h5read(h5_file, "X/indptr")

expr_csr <- new("dgRMatrix",
    p = as.integer(x_indptr),
    j = as.integer(x_indices),
    x = as.numeric(x_data),
    Dim = c(n_cells, n_genes))

expr <- t(expr_csr)  # genes x cells
rownames(expr) <- gene_names

H5close()

# =============================================================================
# 2. Compute GLOBAL statistics (across ALL cells)
# =============================================================================

cat("\n--- Computing Global Statistics ---\n")

# Global mean and std for each gene (across ALL cells)
global_means <- Matrix::rowMeans(expr)
global_vars <- Matrix::rowMeans(expr^2) - global_means^2
global_sds <- sqrt(pmax(global_vars, 0))
global_sds[global_sds < 1e-10] <- 1

cat(sprintf("Global stats computed for %d genes\n", length(global_means)))

# =============================================================================
# 3. Define Sender and Receiver
# =============================================================================

cat("\n--- Defining Sender/Receiver ---\n")

sender_type <- "macrophage"
receiver_type <- "CAF"

sender_mask <- cell_types == sender_type
receiver_mask <- cell_types == receiver_type

sender_idx <- which(sender_mask)
receiver_idx <- which(receiver_mask)

n_senders <- length(sender_idx)
n_receivers <- length(receiver_idx)

cat(sprintf("Sender (%s): %d cells\n", sender_type, n_senders))
cat(sprintf("Receiver (%s): %d cells\n", receiver_type, n_receivers))

sender_coords <- coords_um[sender_idx, ]
receiver_coords <- coords_um[receiver_idx, ]

# =============================================================================
# 4. Prepare Factor and Receiver Expression with GLOBAL normalization
# =============================================================================

cat("\n--- Preparing Expression (Global Normalization) ---\n")

factor_gene <- "TGFB1"
tgfb1_idx <- which(gene_names == factor_gene)

# Factor expression in senders - normalize with GLOBAL stats
factor_expr_raw <- as.vector(expr[factor_gene, sender_idx])
factor_global_mean <- global_means[tgfb1_idx]
factor_global_sd <- global_sds[tgfb1_idx]
factor_expr <- (factor_expr_raw - factor_global_mean) / factor_global_sd

cat(sprintf("TGFB1 global mean=%.4f, sd=%.4f\n", factor_global_mean, factor_global_sd))
cat(sprintf("TGFB1 in senders: raw mean=%.4f, normalized mean=%.4f\n",
    mean(factor_expr_raw), mean(factor_expr)))

# Receiver expression (will normalize with global stats in loop)
receiver_expr <- expr[, receiver_idx]

# =============================================================================
# 5. Compute I_ND at Multiple Radii
# =============================================================================

cat("\n--- Computing I_ND (Global Normalization) ---\n")

radii <- c(10, 20, 30, 50, 100, 200, 300, 500)
n_radii <- length(radii)

I_ND_matrix <- matrix(NA, nrow = n_genes, ncol = n_radii)
rownames(I_ND_matrix) <- gene_names
colnames(I_ND_matrix) <- paste0("r", radii)

cat(sprintf("Radii: %s μm\n\n", paste(radii, collapse = ", ")))

for (r in seq_along(radii)) {
    radius <- radii[r]
    inner_radius <- if (r == 1) 0 else radii[r - 1]
    sigma <- radius / 3

    cat(sprintf("  Radius %d μm (ring %d-%d)...", radius, inner_radius, radius))

    W <- create_weights_sc(
        sender_coords, receiver_coords,
        radius = radius,
        inner_radius = inner_radius,
        sigma = sigma
    )

    mean_conn <- length(W@x) / n_senders

    if (length(W@x) < 10) {
        cat(sprintf(" skipped (%.1f conn)\n", mean_conn))
        next
    }

    # Compute spatial lag (raw)
    lag_G_raw <- as.matrix(W %*% t(receiver_expr))  # n_senders x n_genes

    # Normalize spatial lag with GLOBAL mean and sd
    lag_G <- sweep(lag_G_raw, 2, global_means, "-")
    lag_G <- sweep(lag_G, 2, global_sds, "/")

    # Compute I_ND
    for (g in seq_len(n_genes)) {
        lag_g <- lag_G[, g]
        norm_f <- sqrt(sum(factor_expr^2))
        norm_lag <- sqrt(sum(lag_g^2))

        if (norm_f < 1e-10 || norm_lag < 1e-10 || is.na(norm_lag)) {
            I_ND_matrix[g, r] <- NA
        } else {
            I_ND_matrix[g, r] <- sum(factor_expr * lag_g) / (norm_f * norm_lag)
        }
    }

    cat(sprintf(" done (%.1f conn/sender)\n", mean_conn))
}

# =============================================================================
# 6. Compare with Previous Results
# =============================================================================

cat("\n--- Loading Previous Results ---\n")

prev_results <- read.csv("output/results_macrophage_to_CAF_TGFB1_annular_20250818_150423.csv")
prev_q0 <- prev_results[prev_results$min_expr_quantile == 0, ]

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
    if (sum(valid) < 50) next

    cor_pearson <- cor(new_values[valid], prev_values[valid])
    cor_spearman <- cor(new_values[valid], prev_values[valid], method = "spearman")
    ratio <- mean(abs(prev_values[valid])) / mean(abs(new_values[valid]))

    # Also compute RMSE
    rmse <- sqrt(mean((new_values[valid] - prev_values[valid])^2))

    cat(sprintf("Radius %d μm:\n", radius))
    cat(sprintf("  Pearson r:   %.4f\n", cor_pearson))
    cat(sprintf("  Spearman ρ:  %.4f\n", cor_spearman))
    cat(sprintf("  Scale ratio: %.4f\n", ratio))
    cat(sprintf("  RMSE:        %.4f\n\n", rmse))
}

# =============================================================================
# 7. Top Genes Comparison
# =============================================================================

cat("--- Top 20 Genes (r=10 μm) ---\n\n")

new_r10 <- I_ND_matrix[, "r10"]
prev_r10 <- prev_q0[prev_q0$radius == 10, ]

common_genes <- intersect(names(new_r10[!is.na(new_r10)]), prev_r10$target_gene)
new_r10_common <- new_r10[common_genes]
prev_idx <- match(common_genes, prev_r10$target_gene)
prev_r10_common <- prev_r10$morans_i_normalized[prev_idx]
names(prev_r10_common) <- common_genes

cat("Top 20 by NEW I_ND:\n")
top_new <- head(sort(new_r10_common, decreasing = TRUE), 20)
for (i in seq_along(top_new)) {
    gene <- names(top_new)[i]
    diff <- top_new[i] - prev_r10_common[gene]
    cat(sprintf("  %2d. %-10s new=%.4f prev=%.4f diff=%+.4f\n",
        i, gene, top_new[i], prev_r10_common[gene], diff))
}

cat("\nTop 20 by PREVIOUS I_ND:\n")
top_prev <- head(sort(prev_r10_common, decreasing = TRUE), 20)
for (i in seq_along(top_prev)) {
    gene <- names(top_prev)[i]
    diff <- new_r10_common[gene] - top_prev[i]
    cat(sprintf("  %2d. %-10s prev=%.4f new=%.4f diff=%+.4f\n",
        i, gene, top_prev[i], new_r10_common[gene], diff))
}

# Save
write.csv(I_ND_matrix, "output/coad_tgfb1_h5ad_global_IND.csv")
cat("\n\nResults saved to: output/coad_tgfb1_h5ad_global_IND.csv\n")

cat("\n=== Analysis Complete ===\n")
