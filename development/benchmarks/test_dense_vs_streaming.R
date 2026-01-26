#!/usr/bin/env Rscript
# Compare dense BLAS vs streaming for large radii

library(sigdiscov)

cat("=== Dense BLAS vs Streaming Comparison ===\n\n")

# Load CosMx data
cat("Loading CosMx data...\n")
t0 <- Sys.time()

meta <- read.csv("dataset/cosmx/Lung5_Rep1_meta.csv")
vst <- read.csv("dataset/cosmx/Lung5_Rep1_vst_v2_rpy2.csv", row.names = 1)

cat(sprintf("Loaded: %d genes x %d cells (%.1f sec)\n\n",
            nrow(vst), ncol(vst), difftime(Sys.time(), t0, units = "secs")))

# Use subset for initial comparison
set.seed(42)
n_cells <- 20000
n_genes <- 50
cell_idx <- sample(ncol(vst), n_cells)
gene_idx <- sample(nrow(vst), n_genes)

data <- as.matrix(vst[gene_idx, cell_idx])
coords <- as.matrix(meta[cell_idx, c("sdimx", "sdimy")])

cat(sprintf("Test data: %d genes x %d cells\n\n", n_genes, n_cells))

# Test at different radii
radii_um <- c(100, 500, 1000, 2000)
radii_mm <- radii_um / 1000

for (i in seq_along(radii_mm)) {
    r <- radii_mm[i]
    cat(sprintf("--- Radius = %d um (%.3f mm) ---\n", radii_um[i], r))

    # Streaming
    t1 <- Sys.time()
    result_stream <- pairwise_moran_streaming_cpp(data, coords, r, r/3, verbose = FALSE)
    time_stream <- difftime(Sys.time(), t1, units = "secs")

    # Dense BLAS
    t1 <- Sys.time()
    result_dense <- pairwise_moran_dense_cpp(data, coords, r, r/3, chunk_size = 500, verbose = FALSE)
    time_dense <- difftime(Sys.time(), t1, units = "secs")

    # Compare
    avg_neighbors <- result_stream$n_edges / n_cells
    correlation <- cor(as.vector(result_stream$moran), as.vector(result_dense$moran))

    cat(sprintf("  Avg neighbors/cell: %.0f\n", avg_neighbors))
    cat(sprintf("  Streaming: %.2f sec\n", time_stream))
    cat(sprintf("  Dense BLAS: %.2f sec\n", time_dense))
    cat(sprintf("  Speedup: %.1fx\n", as.numeric(time_stream) / as.numeric(time_dense)))
    cat(sprintf("  Correlation: %.6f\n\n", correlation))
}

cat("=== Full Scale Test (all genes, all cells) ===\n\n")

# Full data
data_full <- as.matrix(vst)
coords_full <- as.matrix(meta[, c("sdimx", "sdimy")])
n_cells_full <- ncol(data_full)
n_genes_full <- nrow(data_full)

cat(sprintf("Full data: %d genes x %d cells\n\n", n_genes_full, n_cells_full))

# Test at 600 um (where streaming was slow)
r_test <- 0.6  # 600 um in mm

cat(sprintf("--- Radius = 600 um ---\n"))

t1 <- Sys.time()
result <- pairwise_moran_dense_cpp(data_full, coords_full, r_test, r_test/3,
                                    chunk_size = 500, verbose = TRUE)
time_dense <- difftime(Sys.time(), t1, units = "secs")

cat(sprintf("\nDense BLAS time: %.1f sec (%.1f min)\n", time_dense, time_dense/60))
cat(sprintf("Avg neighbors/cell: %.0f\n", result$n_edges / n_cells_full))

cat("\nDone.\n")
