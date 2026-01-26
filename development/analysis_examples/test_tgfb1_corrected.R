#!/usr/bin/env Rscript
# TGFB1 Factor Analysis with CORRECTED radii
# Coordinates are in ~mm, not μm!

library(sigdiscov)

cat("=== TGFB1 Factor Analysis (CORRECTED RADII) ===\n\n")

# Load data
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

# Coordinate analysis
cat("\n--- Coordinate Analysis ---\n")
cat(sprintf("X range: %.4f - %.4f (span: %.4f)\n", 
    min(coords[,1]), max(coords[,1]), diff(range(coords[,1]))))
cat(sprintf("Y range: %.4f - %.4f (span: %.4f)\n", 
    min(coords[,2]), max(coords[,2]), diff(range(coords[,2]))))

# Based on median NN distance ~0.027 and typical cell diameter ~15um,
# coordinates likely in mm. 0.027 units ≈ 27 μm
# So radius should be in same units:
# - Short range (50-100 μm): 0.05-0.10 units
# - Medium range (200-300 μm): 0.20-0.30 units

cat("\n--- Defining Sender/Receiver ---\n")
sender_type <- "macrophage"
receiver_type <- "fibroblast"
factor_gene <- "TGFB1"

sender_mask <- meta$cellType == sender_type
receiver_mask <- meta$cellType == receiver_type

cat(sprintf("Sender type: %s (%d cells)\n", sender_type, sum(sender_mask)))
cat(sprintf("Receiver type: %s (%d cells)\n", receiver_type, sum(receiver_mask)))

# High-TGFB1 senders
tgfb1_expr <- expr[factor_gene, sender_mask]
threshold <- quantile(tgfb1_expr, 0.25)
high_tgfb1_mask <- sender_mask & (expr[factor_gene, ] >= threshold)

sender_idx <- which(high_tgfb1_mask)
receiver_idx <- which(receiver_mask)

cat(sprintf("High-TGFB1 senders: %d cells\n", length(sender_idx)))

sender_coords <- coords[sender_idx, ]
receiver_coords <- coords[receiver_idx, ]

factor_expr_raw <- expr[factor_gene, sender_idx]
factor_expr <- (factor_expr_raw - mean(factor_expr_raw)) / sd(factor_expr_raw)

receiver_expr <- expr[, receiver_idx]

# CORRECTED radii (in coordinate units, likely mm)
# 0.05 ≈ 50 μm, 0.10 ≈ 100 μm, etc.
radii <- c(0.02, 0.05, 0.10, 0.15, 0.20, 0.30, 0.50, 1.00)
n_radii <- length(radii)
n_genes <- nrow(expr)
gene_names <- rownames(expr)

cat("\n--- Computing I_ND (CORRECTED radii) ---\n")
cat(sprintf("Radii: %s units\n", paste(radii, collapse = ", ")))

I_ND_matrix <- matrix(NA, nrow = n_genes, ncol = n_radii)
rownames(I_ND_matrix) <- gene_names
colnames(I_ND_matrix) <- paste0("r", radii)

for (r in seq_along(radii)) {
    radius <- radii[r]
    sigma <- radius / 3
    
    cat(sprintf("  Radius %.2f units...", radius))
    
    W <- create_weights_sc(sender_coords, receiver_coords, radius = radius, sigma = sigma)
    
    # Count connections
    nnz_per_row <- diff(W@p)
    mean_conn <- mean(nnz_per_row)
    
    if (length(W@x) < 10) {
        cat(sprintf(" skipped (%.0f conn/sender)\n", mean_conn))
        next
    }
    
    receiver_expr_norm <- t(scale(t(receiver_expr)))
    receiver_expr_norm[is.na(receiver_expr_norm)] <- 0
    
    lag_G <- as.matrix(W %*% t(receiver_expr_norm))
    
    for (g in seq_len(n_genes)) {
        lag_g <- lag_G[, g]
        I_ND_matrix[g, r] <- sum(factor_expr * lag_g) /
            (sqrt(sum(factor_expr^2)) * sqrt(sum(lag_g^2) + 1e-10))
    }
    
    cat(sprintf(" done (%.0f connections/sender)\n", mean_conn))
}

# Compute Delta I
cat("\n--- Computing Delta I ---\n")
delta_I <- numeric(n_genes)
I_max <- numeric(n_genes)
I_r1 <- I_ND_matrix[, 1]

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
    
    sign <- if (!is.na(I_short) && !is.na(I_long) && I_short >= I_long) 1 else -1
    delta_I[g] <- sign * (max_val - min_val)
    I_max[g] <- max_val
}

names(delta_I) <- gene_names
names(I_max) <- gene_names
names(I_r1) <- gene_names

# Compute Pairwise Moran's I at r=0.10 (~100 μm)
cat("\n--- Computing Pairwise Moran's I (r=0.10) ---\n")
radius_moran <- 0.10
sigma_moran <- radius_moran / 3

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

rownames(result_moran$moran) <- gene_names
colnames(result_moran$moran) <- gene_names
tgfb1_moran <- result_moran$moran[factor_gene, ]
names(tgfb1_moran) <- gene_names

# Load signatures
cat("\n--- Loading Signatures ---\n")

# CytoSig
sig_data <- read.table("dataset/db/signature.centroid", header = TRUE, 
                       row.names = 1, sep = "\t", check.names = FALSE)
cytosig_tgfb1 <- sig_data[, "TGFB1"]
names(cytosig_tgfb1) <- rownames(sig_data)

# SecAct
secact_file <- "dataset/db/AllSigFilteredBy_MoranI_TCGA_ICGC_0.25_ds3.tsv.gz"
if (file.exists(secact_file)) {
    secact <- read.table(gzfile(secact_file), header = TRUE, row.names = 1, sep = "\t")
    secact_tgfb1 <- secact[, "TGFB1"]
    names(secact_tgfb1) <- rownames(secact)
    cat(sprintf("SecAct TGFB1: %d genes\n", length(secact_tgfb1)))
}

cat(sprintf("CytoSig TGFB1: %d genes\n", length(cytosig_tgfb1)))

# Correlation with CytoSig
common_genes_cytosig <- intersect(gene_names, names(cytosig_tgfb1))
cat(sprintf("Common genes (CytoSig): %d\n", length(common_genes_cytosig)))

I_ND_r1_common <- I_r1[common_genes_cytosig]
I_max_common <- I_max[common_genes_cytosig]
delta_I_common <- delta_I[common_genes_cytosig]
moran_common <- tgfb1_moran[common_genes_cytosig]
cytosig_common <- cytosig_tgfb1[common_genes_cytosig]

valid_idx <- !is.na(I_ND_r1_common) & !is.na(cytosig_common)

cat(sprintf("\n=== Correlation with CytoSig TGFB1 (n=%d genes) ===\n", sum(valid_idx)))

cor_I_ND_r1 <- cor(I_ND_r1_common[valid_idx], cytosig_common[valid_idx], method = "pearson")
cor_I_max <- cor(I_max_common[valid_idx], cytosig_common[valid_idx], method = "pearson")
cor_delta_I <- cor(delta_I_common[valid_idx], cytosig_common[valid_idx], method = "pearson")
cor_moran <- cor(moran_common[valid_idx], cytosig_common[valid_idx], method = "pearson")

cat(sprintf("  I_ND (r=0.02):  r = %.4f\n", cor_I_ND_r1))
cat(sprintf("  I_ND (max):     r = %.4f\n", cor_I_max))
cat(sprintf("  Delta I:        r = %.4f\n", cor_delta_I))
cat(sprintf("  Moran's I:      r = %.4f\n", cor_moran))

# Spearman
cat("\nSpearman correlations:\n")
cat(sprintf("  I_ND (r=0.02):  rho = %.4f\n", 
    cor(I_ND_r1_common[valid_idx], cytosig_common[valid_idx], method = "spearman")))
cat(sprintf("  Moran's I:      rho = %.4f\n", 
    cor(moran_common[valid_idx], cytosig_common[valid_idx], method = "spearman")))

# Correlation with SecAct
if (exists("secact_tgfb1")) {
    common_genes_secact <- intersect(gene_names, names(secact_tgfb1))
    cat(sprintf("\n=== Correlation with SecAct TGFB1 (n=%d genes) ===\n", length(common_genes_secact)))
    
    valid_secact <- common_genes_secact[!is.na(I_r1[common_genes_secact])]
    
    cor_I_ND_secact <- cor(I_r1[valid_secact], secact_tgfb1[valid_secact], method = "pearson")
    cor_moran_secact <- cor(tgfb1_moran[valid_secact], secact_tgfb1[valid_secact], method = "pearson")
    
    cat(sprintf("  I_ND (r=0.02):  r = %.4f\n", cor_I_ND_secact))
    cat(sprintf("  Moran's I:      r = %.4f\n", cor_moran_secact))
}

# Variance comparison
cat("\n=== Metric Variance ===\n")
cat(sprintf("I_ND variance:     %.6e\n", var(I_ND_r1_common[valid_idx], na.rm = TRUE)))
cat(sprintf("Moran's I variance: %.6e\n", var(moran_common[valid_idx], na.rm = TRUE)))
cat(sprintf("Delta I variance:   %.6e\n", var(delta_I_common[valid_idx], na.rm = TRUE)))

# Top genes
cat("\n=== Top 15 Genes by Each Metric ===\n")

result_df <- data.frame(
    gene = common_genes_cytosig,
    I_ND_r0.02 = I_ND_r1_common,
    delta_I = delta_I_common,
    moran_I = moran_common,
    cytosig = cytosig_common,
    stringsAsFactors = FALSE
)
result_df <- result_df[valid_idx, ]

cat("\nTop by Moran's I:\n")
top_moran <- result_df[order(-result_df$moran_I), ][1:15, ]
print(top_moran[, c("gene", "moran_I", "I_ND_r0.02", "cytosig")], row.names = FALSE)

cat("\nTop by I_ND:\n")
top_IND <- result_df[order(-result_df$I_ND_r0.02), ][1:15, ]
print(top_IND[, c("gene", "I_ND_r0.02", "moran_I", "cytosig")], row.names = FALSE)

cat("\nTop by CytoSig (expected targets):\n")
top_cytosig <- result_df[order(-result_df$cytosig), ][1:15, ]
print(top_cytosig[, c("gene", "cytosig", "moran_I", "I_ND_r0.02")], row.names = FALSE)

cat("\n=== Analysis Complete ===\n")
