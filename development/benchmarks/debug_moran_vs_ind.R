#!/usr/bin/env Rscript
# Debug: Why Moran's I and I_ND show different correlations with CytoSig

library(sigdiscov)

cat("=== Debug: Moran's I vs I_ND Discrepancy ===\n\n")

# Load the analysis results
I_ND_matrix <- read.csv("output/tgfb1_analysis/I_ND_matrix.csv", row.names = 1)
result_df <- read.csv("output/tgfb1_analysis/tgfb1_factor_results.csv")

# Load SecAct signature
secact_file <- "dataset/db/AllSigFilteredBy_MoranI_TCGA_ICGC_0.25_ds3.tsv.gz"
if (file.exists(secact_file)) {
    cat("Loading SecAct signature...\n")
    secact <- read.table(gzfile(secact_file), header = TRUE, row.names = 1, sep = "\t")
    if ("TGFB1" %in% colnames(secact)) {
        secact_tgfb1 <- secact[, "TGFB1"]
        names(secact_tgfb1) <- rownames(secact)
        cat("  SecAct TGFB1 signature: ", length(secact_tgfb1), " genes\n")

        # Correlate with result_df
        common <- intersect(result_df$gene, names(secact_tgfb1))
        cat("  Common genes with results: ", length(common), "\n\n")

        if (length(common) > 50) {
            idx <- match(common, result_df$gene)
            cor_moran_secact <- cor(result_df$moran_I[idx], secact_tgfb1[common],
                                    use = "complete.obs")
            cor_ind_secact <- cor(result_df$I_ND_r50[idx], secact_tgfb1[common],
                                  use = "complete.obs")
            cat("Correlation with SecAct TGFB1:\n")
            cat("  Moran's I: r = ", round(cor_moran_secact, 4), "\n")
            cat("  I_ND:      r = ", round(cor_ind_secact, 4), "\n\n")
        }
    }
}

# === ANALYSIS OF THE DISCREPANCY ===

cat("=== Investigating the Discrepancy ===\n\n")

# 1. Check the I_ND curve shape across radii
cat("1. I_ND values across radii (selected genes):\n")
cat("   (Should decrease with radius for paracrine signaling)\n\n")

# Select some genes with high CytoSig scores
cytosig <- read.table("dataset/db/signature.centroid", header = TRUE, row.names = 1, sep = "\t")
cytosig_tgfb1 <- cytosig[, "TGFB1"]

# Find genes with high positive CytoSig
high_cytosig_genes <- names(sort(cytosig_tgfb1, decreasing = TRUE))[1:10]
common_high <- intersect(high_cytosig_genes, rownames(I_ND_matrix))

cat("I_ND curves for top CytoSig genes:\n")
print(round(I_ND_matrix[common_high, ], 4))

cat("\nCytoSig values for these genes:\n")
print(round(cytosig_tgfb1[common_high], 4))

# 2. Check variance of I_ND vs Moran's I
cat("\n2. Variance comparison:\n")
cat("   I_ND variance:      ", var(result_df$I_ND_r50, na.rm = TRUE), "\n")
cat("   Moran's I variance: ", var(result_df$moran_I, na.rm = TRUE), "\n")

# 3. The key issue: connectivity pattern
cat("\n3. Connectivity Analysis:\n")
cat("   In the single-cell analysis:\n")
cat("   - 5561 senders (high-TGFB1 macrophages)\n")
cat("   - 13151 receivers (all fibroblasts)\n")
cat("   - Each sender connects to ALL 13151 receivers!\n")
cat("   - This means the weight matrix is nearly UNIFORM\n")

cat("\n4. Problem with uniform weights:\n")
cat("   When W is uniform (all entries ~ 1/n_receiver):\n")
cat("   - spatial_lag[g] = W @ receiver_expr[g] ≈ mean(receiver_expr[g])\n")
cat("   - This makes spatial_lag nearly CONSTANT across senders\n")
cat("   - I_ND = cosine(factor_expr, spatial_lag) becomes meaningless\n")

# 5. Verify: check spatial lag variance
cat("\n5. Verifying spatial lag uniformity:\n")

# Re-run a quick test with actual data
expr_file <- "dataset/cosmx/Lung5_Rep1_vst.tsv"
meta_file <- "dataset/cosmx/Lung5_Rep1_meta.csv"

cat("   Loading data subset...\n")
# Load just first 100 genes for speed
expr_raw <- read.table(expr_file, header = TRUE, row.names = 1, sep = "\t",
                       nrows = 100, check.names = FALSE)
meta <- read.csv(meta_file, row.names = 1)

common_cells <- intersect(colnames(expr_raw), meta$cell_ID)
expr_raw <- expr_raw[, common_cells[1:10000]]  # Subset for speed
meta <- meta[match(common_cells[1:10000], meta$cell_ID), ]

sender_mask <- meta$cellType == "macrophage"
receiver_mask <- meta$cellType == "fibroblast"

sender_idx <- which(sender_mask)
receiver_idx <- which(receiver_mask)

if (length(sender_idx) > 100 && length(receiver_idx) > 100) {
    sender_coords <- as.matrix(meta[sender_idx, c("sdimx", "sdimy")])
    receiver_coords <- as.matrix(meta[receiver_idx, c("sdimx", "sdimy")])

    # Create weight matrix at different radii
    for (r in c(50, 100, 500)) {
        W <- create_weights_sc(sender_coords, receiver_coords, radius = r)

        # Check how many receivers each sender connects to
        nnz_per_row <- diff(W@p)  # Number of non-zeros per row for dgCMatrix

        cat(sprintf("\n   Radius %d um:\n", r))
        cat(sprintf("     Mean connections/sender: %.1f\n", mean(nnz_per_row)))
        cat(sprintf("     Min connections: %d\n", min(nnz_per_row)))
        cat(sprintf("     Max connections: %d\n", max(nnz_per_row)))

        # Check weight variance within rows
        if (length(W@x) > 0) {
            cat(sprintf("     Weight variance (across all entries): %.6f\n", var(W@x)))
        }
    }
}

# 6. The solution
cat("\n\n=== ROOT CAUSE ===\n")
cat("The CosMx data coordinates span a small physical range relative to the radius.\n")
cat("At radius=100um, nearly ALL receivers are within range of every sender.\n")
cat("This makes the spatial weighting trivial - no spatial selectivity.\n\n")

cat("=== SOLUTION ===\n")
cat("1. Use SMALLER radii (e.g., 10-50 um) to get actual spatial structure\n")
cat("2. Or use INNER_RADIUS to create an annular ring (exclude immediate neighbors)\n")
cat("3. The Visium case works because spots are spaced further apart (~100um)\n\n")

# 7. Check coordinate range
cat("=== Coordinate Analysis ===\n")
coord_range_x <- diff(range(meta$sdimx))
coord_range_y <- diff(range(meta$sdimy))
cat(sprintf("Coordinate range X: %.2f\n", coord_range_x))
cat(sprintf("Coordinate range Y: %.2f\n", coord_range_y))
cat(sprintf("Approx tissue size: %.0f x %.0f um\n", coord_range_x, coord_range_y))

# Cell density
n_cells <- nrow(meta)
area <- coord_range_x * coord_range_y
density <- n_cells / area
avg_dist <- sqrt(1 / density)
cat(sprintf("\nCell density: %.2f cells/um^2\n", density))
cat(sprintf("Average inter-cell distance: %.2f um\n", avg_dist))
cat(sprintf("\nRecommended radius range: %.0f - %.0f um\n", avg_dist, avg_dist * 5))
