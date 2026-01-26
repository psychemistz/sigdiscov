#!/usr/bin/env Rscript
# Test cell type pair analysis on CosMx Lung5_Rep1 data

library(sigdiscov)

# =============================================================================
# Load CosMx Data
# =============================================================================

cat("=== Cell Type Pair Analysis Test on CosMx Data ===\n\n")

expr_file <- "dataset/cosmx/Lung5_Rep1_vst.tsv"
meta_file <- "dataset/cosmx/Lung5_Rep1_meta.csv"

cat("Loading CosMx data...\n")
cat("  Expression:", expr_file, "\n")
cat("  Metadata:", meta_file, "\n")

t1 <- Sys.time()

# Load expression data
if (grepl("\\.tsv$", expr_file)) {
    expr <- as.matrix(read.table(expr_file, header = TRUE, row.names = 1,
                                  sep = "\t", check.names = FALSE))
} else {
    expr <- as.matrix(read.csv(expr_file, row.names = 1, check.names = FALSE))
}
cat(sprintf("  Expression matrix: %d genes x %d cells\n", nrow(expr), ncol(expr)))

# Load metadata
meta <- read.csv(meta_file, row.names = 1)
cat(sprintf("  Metadata: %d cells\n", nrow(meta)))
cat(sprintf("  Columns: %s\n", paste(colnames(meta), collapse = ", ")))

# Align expression and metadata
common_cells <- intersect(colnames(expr), meta$cell_ID)
cat(sprintf("  Common cells: %d\n", length(common_cells)))

expr <- expr[, common_cells]
meta <- meta[match(common_cells, meta$cell_ID), ]

# Create SCData object
data <- list(
    expr = expr,
    coords = as.matrix(meta[, c("sdimx", "sdimy")]),
    cell_types = meta$cellType,
    gene_names = rownames(expr),
    n_genes = nrow(expr),
    n_cells = ncol(expr)
)
class(data) <- c("SCData", "list")

t2 <- Sys.time()
cat(sprintf("Loaded: %d genes x %d cells (%.1f sec)\n\n",
            data$n_genes, data$n_cells, as.numeric(difftime(t2, t1, units = "secs"))))

# Show cell type distribution
cat("Cell type distribution:\n")
ct_table <- sort(table(data$cell_types), decreasing = TRUE)
print(ct_table)
cat("\n")

# =============================================================================
# Test 1: Single Cell Type Pair (small test)
# =============================================================================

cat("--- Test 1: Single Cell Type Pair ---\n")

# Pick two cell types with reasonable counts
sender_type <- names(ct_table)[1]  # Most common
receiver_type <- names(ct_table)[2]  # Second most common

cat(sprintf("Testing: %s -> %s\n", sender_type, receiver_type))

# Extract subsets
sender_idx <- which(data$cell_types == sender_type)
receiver_idx <- which(data$cell_types == receiver_type)

cat(sprintf("  Senders: %d cells\n", length(sender_idx)))
cat(sprintf("  Receivers: %d cells\n", length(receiver_idx)))

# Use first 50 genes for quick test
n_test_genes <- 50
test_genes <- seq_len(min(n_test_genes, data$n_genes))

sender_data <- data$expr[test_genes, sender_idx, drop = FALSE]
receiver_data <- data$expr[test_genes, receiver_idx, drop = FALSE]
sender_coords <- data$coords[sender_idx, , drop = FALSE]
receiver_coords <- data$coords[receiver_idx, , drop = FALSE]

# Run single radius test
cat("\nRunning pairwise_moran_celltype_pair (radius=100)...\n")
t1 <- Sys.time()
result <- pairwise_moran_celltype_pair(
    sender_data = sender_data,
    receiver_data = receiver_data,
    sender_coords = sender_coords,
    receiver_coords = receiver_coords,
    radius = 100,  # 100 um
    method = "streaming",
    verbose = TRUE
)
t2 <- Sys.time()

cat(sprintf("\nCompleted in %.1f sec\n", as.numeric(difftime(t2, t1, units = "secs"))))
cat(sprintf("Moran matrix: %d x %d\n", nrow(result$moran), ncol(result$moran)))
cat(sprintf("Weight sum (S0): %.2f\n", result$weight_sum))
cat(sprintf("Total edges: %d\n", result$n_edges))

# Show top Moran values
moran_vals <- result$moran[upper.tri(result$moran)]
cat(sprintf("\nMoran's I summary:\n"))
cat(sprintf("  Mean: %.4f\n", mean(moran_vals)))
cat(sprintf("  SD: %.4f\n", sd(moran_vals)))
cat(sprintf("  Range: [%.4f, %.4f]\n", min(moran_vals), max(moran_vals)))

# =============================================================================
# Test 2: Multi-radius Analysis
# =============================================================================

cat("\n--- Test 2: Multi-radius Analysis (C++ function) ---\n")

radii <- c(50, 100, 150, 200)
cat(sprintf("Radii: %s um\n", paste(radii, collapse = ", ")))

t1 <- Sys.time()
result_multi <- compute_celltype_pair_moran_cpp(
    sender_data = sender_data,
    receiver_data = receiver_data,
    sender_coords = as.matrix(sender_coords),
    receiver_coords = as.matrix(receiver_coords),
    radii = radii,
    sigma_factor = 1/3,
    verbose = TRUE
)
t2 <- Sys.time()

cat(sprintf("\nCompleted in %.1f sec\n", as.numeric(difftime(t2, t1, units = "secs"))))
cat(sprintf("Delta I matrix: %d x %d\n", nrow(result_multi$delta_I), ncol(result_multi$delta_I)))

# Show delta I summary
delta_vals <- result_multi$delta_I[upper.tri(result_multi$delta_I)]
cat(sprintf("\nDelta I summary:\n"))
cat(sprintf("  Mean: %.4f\n", mean(delta_vals)))
cat(sprintf("  SD: %.4f\n", sd(delta_vals)))
cat(sprintf("  Range: [%.4f, %.4f]\n", min(delta_vals), max(delta_vals)))

# Top gene pairs by delta I
delta_df <- data.frame(
    gene1 = rep(rownames(sender_data), each = length(test_genes)),
    gene2 = rep(rownames(sender_data), times = length(test_genes)),
    delta_I = as.vector(result_multi$delta_I)
)
delta_df <- delta_df[order(-abs(delta_df$delta_I)), ]
cat("\nTop 10 gene pairs by |delta I|:\n")
print(head(delta_df, 10), row.names = FALSE)

# =============================================================================
# Test 3: Full compute_celltype_pair_analysis (subset)
# =============================================================================

cat("\n--- Test 3: Full Analysis (3 cell types, 3 radii) ---\n")

# Create subset data with 3 most common cell types
top3_types <- names(ct_table)[1:3]
subset_idx <- which(data$cell_types %in% top3_types)

# Limit cells per type for faster testing
max_cells_per_type <- 2000
subset_idx_limited <- c()
for (ct in top3_types) {
    ct_idx <- which(data$cell_types == ct)
    if (length(ct_idx) > max_cells_per_type) {
        ct_idx <- sample(ct_idx, max_cells_per_type)
    }
    subset_idx_limited <- c(subset_idx_limited, ct_idx)
}

# Create subset SCData object
data_subset <- list(
    expr = data$expr[test_genes, subset_idx_limited, drop = FALSE],
    coords = data$coords[subset_idx_limited, , drop = FALSE],
    cell_types = data$cell_types[subset_idx_limited],
    gene_names = data$gene_names[test_genes],
    n_genes = length(test_genes),
    n_cells = length(subset_idx_limited)
)
class(data_subset) <- c("SCData", "list")

cat(sprintf("Subset: %d genes x %d cells\n", data_subset$n_genes, data_subset$n_cells))
cat("Cell types:", paste(top3_types, collapse = ", "), "\n")
cat("Pairs to analyze: 9 (3x3)\n")

# Run full analysis
t1 <- Sys.time()
result_full <- compute_celltype_pair_analysis(
    data = data_subset,
    radii = c(50, 100, 200),
    min_cells = 100,
    output_dir = NULL,  # Keep in memory
    verbose = TRUE
)
t2 <- Sys.time()

cat(sprintf("\nTotal time: %.1f sec\n", as.numeric(difftime(t2, t1, units = "secs"))))

# Print result
print(result_full)

# Show pair stats
cat("\nPair statistics:\n")
print(result_full$pair_stats, row.names = FALSE)

# Extract top delta I
cat("\nTop 20 gene pairs by |delta I| (across all cell type pairs):\n")
top_delta <- extract_top_delta_i(result_full, n_top = 20, min_delta_I = 0)
print(top_delta, row.names = FALSE)

cat("\n=== Test Complete ===\n")
