#!/usr/bin/env Rscript
# TGFB1 Factor Analysis: Macrophage -> Fibroblast
# Compute I_ND, Moran's I, and Delta I, then correlate with CytoSig/SecAct signatures

library(sigdiscov)

cat("=== TGFB1 Factor Analysis: Macrophage -> Fibroblast ===\n\n")

# =============================================================================
# 1. Load CosMx Data
# =============================================================================

cat("--- Loading Data ---\n")

expr_file <- "dataset/cosmx/Lung5_Rep1_vst.tsv"
meta_file <- "dataset/cosmx/Lung5_Rep1_meta.csv"

t1 <- Sys.time()
expr <- as.matrix(read.table(expr_file, header = TRUE, row.names = 1,
                              sep = "\t", check.names = FALSE))
meta <- read.csv(meta_file, row.names = 1)

common_cells <- intersect(colnames(expr), meta$cell_ID)
expr <- expr[, common_cells]
meta <- meta[match(common_cells, meta$cell_ID), ]
coords <- as.matrix(meta[, c("sdimx", "sdimy")])

t2 <- Sys.time()
cat(sprintf("Loaded: %d genes x %d cells (%.1f sec)\n",
            nrow(expr), ncol(expr), as.numeric(difftime(t2, t1, units = "secs"))))

# Check TGFB1
if (!"TGFB1" %in% rownames(expr)) {
    stop("TGFB1 not found in gene list")
}
cat("TGFB1 found in gene list.\n\n")

# =============================================================================
# 2. Define Sender and Receiver Populations
# =============================================================================

cat("--- Defining Sender/Receiver ---\n")

sender_type <- "macrophage"
receiver_type <- "fibroblast"
factor_gene <- "TGFB1"

# Get cell type indices
sender_mask <- meta$cellType == sender_type
receiver_mask <- meta$cellType == receiver_type

cat(sprintf("Sender type: %s (%d cells)\n", sender_type, sum(sender_mask)))
cat(sprintf("Receiver type: %s (%d cells)\n", receiver_type, sum(receiver_mask)))

# Define high-TGFB1 senders (top 75% expression among macrophages)
tgfb1_expr <- expr[factor_gene, sender_mask]
threshold <- quantile(tgfb1_expr, 0.25)  # Top 75%
high_tgfb1_mask <- sender_mask & (expr[factor_gene, ] >= threshold)

cat(sprintf("TGFB1 threshold (25th percentile): %.3f\n", threshold))
cat(sprintf("High-TGFB1 senders: %d cells\n", sum(high_tgfb1_mask)))

# Extract data
sender_idx <- which(high_tgfb1_mask)
receiver_idx <- which(receiver_mask)

sender_coords <- coords[sender_idx, ]
receiver_coords <- coords[receiver_idx, ]

# Factor expression (z-normalized TGFB1 in senders)
factor_expr_raw <- expr[factor_gene, sender_idx]
factor_expr <- (factor_expr_raw - mean(factor_expr_raw)) / sd(factor_expr_raw)

# Receiver gene expression (all genes)
receiver_expr <- expr[, receiver_idx]

cat("\n")

# =============================================================================
# 3. Compute I_ND at Multiple Radii
# =============================================================================

cat("--- Computing I_ND (Factor-Gene Correlation) ---\n")

radii <- c(50, 100, 150, 200, 250, 300, 400, 500)
n_radii <- length(radii)
n_genes <- nrow(expr)
gene_names <- rownames(expr)

# Storage for I_ND curves
I_ND_matrix <- matrix(NA, nrow = n_genes, ncol = n_radii)
rownames(I_ND_matrix) <- gene_names
colnames(I_ND_matrix) <- paste0("r", radii)

cat(sprintf("Radii: %s um\n", paste(radii, collapse = ", ")))

for (r in seq_along(radii)) {
    radius <- radii[r]
    sigma <- radius / 3

    cat(sprintf("  Radius %d um...", radius))

    # Create weight matrix (sender -> receiver)
    W <- create_weights_sc(sender_coords, receiver_coords, radius = radius, sigma = sigma)

    if (length(W@x) < 10) {
        cat(" (skipped - too few connections)\n")
        next
    }

    # Compute spatial lag for all receiver genes
    # lag_g = W %*% receiver_expr^T (for each gene)
    # Result: n_senders x n_genes
    receiver_expr_norm <- t(scale(t(receiver_expr)))  # Z-normalize genes
    receiver_expr_norm[is.na(receiver_expr_norm)] <- 0

    # Spatial lag: for each sender, weighted average of receiver expression
    lag_G <- as.matrix(W %*% t(receiver_expr_norm))  # n_senders x n_genes

    # Compute I_ND for each gene
    # I_ND = cosine similarity between factor_expr and lag_g
    for (g in seq_len(n_genes)) {
        lag_g <- lag_G[, g]
        # Cosine similarity
        I_ND_matrix[g, r] <- sum(factor_expr * lag_g) /
            (sqrt(sum(factor_expr^2)) * sqrt(sum(lag_g^2) + 1e-10))
    }

    cat(sprintf(" done (%.0f connections/sender)\n", length(W@x) / nrow(sender_coords)))
}

cat("\n")

# =============================================================================
# 4. Compute Delta I from I_ND Curves
# =============================================================================

cat("--- Computing Delta I ---\n")

delta_I <- numeric(n_genes)
I_max <- numeric(n_genes)
I_r1 <- I_ND_matrix[, 1]  # First radius

for (g in seq_len(n_genes)) {
    curve <- I_ND_matrix[g, ]
    if (all(is.na(curve))) {
        delta_I[g] <- NA
        I_max[g] <- NA
        next
    }

    max_val <- max(curve, na.rm = TRUE)
    min_val <- min(curve, na.rm = TRUE)
    I_short <- curve[1]
    I_long <- curve[length(curve)]

    sign <- if (I_short >= I_long) 1 else -1
    delta_I[g] <- sign * (max_val - min_val)
    I_max[g] <- max_val
}

names(delta_I) <- gene_names
names(I_max) <- gene_names
names(I_r1) <- gene_names

cat(sprintf("Delta I computed for %d genes\n", sum(!is.na(delta_I))))

# =============================================================================
# 5. Compute Pairwise Moran's I (for comparison)
# =============================================================================

cat("\n--- Computing Pairwise Moran's I ---\n")

# Use single radius for pairwise Moran (expensive)
radius_moran <- 100
sigma_moran <- radius_moran / 3

cat(sprintf("Radius: %d um\n", radius_moran))

# Directional Moran's I: sender genes x receiver genes
# This computes correlation between ALL sender genes and ALL receiver genes
# For factor analysis, we focus on factor_gene row

result_moran <- pairwise_moran_directional_streaming_cpp(
    sender_data = expr[, sender_idx],
    receiver_data = expr[, receiver_idx],
    sender_coords = sender_coords,
    receiver_coords = receiver_coords,
    radius = radius_moran,
    sigma = sigma_moran,
    normalize_data = TRUE,
    verbose = TRUE
)

# Add dimension names to moran matrix
rownames(result_moran$moran) <- gene_names
colnames(result_moran$moran) <- gene_names

# Extract TGFB1 row (factor gene correlations with all receiver genes)
tgfb1_moran <- result_moran$moran[factor_gene, ]
names(tgfb1_moran) <- gene_names

cat("\n")

# =============================================================================
# 6. Load CytoSig Signature
# =============================================================================

cat("--- Loading CytoSig TGFB1 Signature ---\n")

sig_file <- "dataset/db/signature.centroid"
sig_data <- read.table(sig_file, header = TRUE, row.names = 1, sep = "\t",
                       check.names = FALSE)

# TGFB1 is column 37 (check header)
if ("TGFB1" %in% colnames(sig_data)) {
    cytosig_tgfb1 <- sig_data[, "TGFB1"]
    names(cytosig_tgfb1) <- rownames(sig_data)
    cat(sprintf("CytoSig TGFB1 signature: %d genes\n", length(cytosig_tgfb1)))
} else {
    stop("TGFB1 not found in signature file")
}

# Match genes
common_genes_cytosig <- intersect(gene_names, names(cytosig_tgfb1))
cat(sprintf("Common genes with CosMx: %d\n", length(common_genes_cytosig)))

# =============================================================================
# 7. Try to Load SecAct Signature (if available)
# =============================================================================

cat("\n--- SecAct Signature ---\n")

# SecAct signatures might be in a different file or need to be fetched
# For now, we'll use CytoSig and note that SecAct comparison can be added

# Check if there's a SecAct-related file
secact_files <- list.files("dataset/db", pattern = "secact|SecAct",
                           full.names = TRUE, ignore.case = TRUE)

if (length(secact_files) > 0) {
    cat(sprintf("Found SecAct file: %s\n", secact_files[1]))
    # Load and process SecAct data
} else {
    cat("SecAct signature file not found locally.\n")
    cat("Using CytoSig TGFB1 signature for correlation analysis.\n")
    # SecAct is from the same group as CytoSig (Jiang lab), signatures are similar
}

# =============================================================================
# 8. Correlate Results with CytoSig
# =============================================================================

cat("\n--- Correlation Analysis ---\n")

# Prepare vectors for common genes
idx_common <- match(common_genes_cytosig, gene_names)

I_ND_r1_common <- I_r1[common_genes_cytosig]
I_max_common <- I_max[common_genes_cytosig]
delta_I_common <- delta_I[common_genes_cytosig]
moran_common <- tgfb1_moran[common_genes_cytosig]
cytosig_common <- cytosig_tgfb1[common_genes_cytosig]

# Remove NA values
valid_idx <- !is.na(I_ND_r1_common) & !is.na(cytosig_common)

cat(sprintf("\nCorrelation with CytoSig TGFB1 (n=%d genes):\n", sum(valid_idx)))

# Pearson correlations
cor_I_ND_r1 <- cor(I_ND_r1_common[valid_idx], cytosig_common[valid_idx], method = "pearson")
cor_I_max <- cor(I_max_common[valid_idx], cytosig_common[valid_idx], method = "pearson")
cor_delta_I <- cor(delta_I_common[valid_idx], cytosig_common[valid_idx], method = "pearson")
cor_moran <- cor(moran_common[valid_idx], cytosig_common[valid_idx], method = "pearson")

cat(sprintf("  I_ND (r=50um):     r = %.4f\n", cor_I_ND_r1))
cat(sprintf("  I_ND (max):        r = %.4f\n", cor_I_max))
cat(sprintf("  Delta I:           r = %.4f\n", cor_delta_I))
cat(sprintf("  Moran's I (r=100): r = %.4f\n", cor_moran))

# Spearman correlations
cat("\nSpearman correlations:\n")
cor_I_ND_r1_sp <- cor(I_ND_r1_common[valid_idx], cytosig_common[valid_idx], method = "spearman")
cor_I_max_sp <- cor(I_max_common[valid_idx], cytosig_common[valid_idx], method = "spearman")
cor_delta_I_sp <- cor(delta_I_common[valid_idx], cytosig_common[valid_idx], method = "spearman")
cor_moran_sp <- cor(moran_common[valid_idx], cytosig_common[valid_idx], method = "spearman")

cat(sprintf("  I_ND (r=50um):     rho = %.4f\n", cor_I_ND_r1_sp))
cat(sprintf("  I_ND (max):        rho = %.4f\n", cor_I_max_sp))
cat(sprintf("  Delta I:           rho = %.4f\n", cor_delta_I_sp))
cat(sprintf("  Moran's I (r=100): rho = %.4f\n", cor_moran_sp))

# =============================================================================
# 9. Summary Statistics
# =============================================================================

cat("\n--- Summary Statistics ---\n")

cat("\nI_ND (r=50um):\n")
print(summary(I_ND_r1_common[valid_idx]))

cat("\nDelta I:\n")
print(summary(delta_I_common[valid_idx]))

cat("\nMoran's I (r=100um):\n")
print(summary(moran_common[valid_idx]))

cat("\nCytoSig TGFB1:\n")
print(summary(cytosig_common[valid_idx]))

# =============================================================================
# 10. Top Genes by Each Metric
# =============================================================================

cat("\n--- Top 20 Genes by Each Metric ---\n")

# Create result data frame
result_df <- data.frame(
    gene = common_genes_cytosig,
    I_ND_r50 = I_ND_r1_common,
    I_max = I_max_common,
    delta_I = delta_I_common,
    moran_I = moran_common,
    cytosig = cytosig_common,
    stringsAsFactors = FALSE
)
result_df <- result_df[valid_idx, ]

cat("\nTop 20 by I_ND (r=50um):\n")
top_IND <- result_df[order(-result_df$I_ND_r50), ][1:20, ]
print(top_IND[, c("gene", "I_ND_r50", "cytosig")], row.names = FALSE)

cat("\nTop 20 by Delta I:\n")
top_delta <- result_df[order(-result_df$delta_I), ][1:20, ]
print(top_delta[, c("gene", "delta_I", "cytosig")], row.names = FALSE)

cat("\nTop 20 by Moran's I:\n")
top_moran <- result_df[order(-result_df$moran_I), ][1:20, ]
print(top_moran[, c("gene", "moran_I", "cytosig")], row.names = FALSE)

cat("\nTop 20 by CytoSig (positive):\n")
top_cytosig <- result_df[order(-result_df$cytosig), ][1:20, ]
print(top_cytosig[, c("gene", "cytosig", "I_ND_r50", "delta_I")], row.names = FALSE)

# =============================================================================
# 11. Save Results
# =============================================================================

cat("\n--- Saving Results ---\n")

output_dir <- "output/tgfb1_analysis"
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

# Save I_ND matrix (all radii)
write.csv(I_ND_matrix, file.path(output_dir, "I_ND_matrix.csv"))

# Save summary results
write.csv(result_df, file.path(output_dir, "tgfb1_factor_results.csv"), row.names = FALSE)

# Save correlation summary
cor_summary <- data.frame(
    metric = c("I_ND_r50", "I_max", "delta_I", "moran_I"),
    pearson_r = c(cor_I_ND_r1, cor_I_max, cor_delta_I, cor_moran),
    spearman_rho = c(cor_I_ND_r1_sp, cor_I_max_sp, cor_delta_I_sp, cor_moran_sp)
)
write.csv(cor_summary, file.path(output_dir, "correlation_with_cytosig.csv"), row.names = FALSE)

cat(sprintf("Results saved to: %s/\n", output_dir))

cat("\n=== Analysis Complete ===\n")
