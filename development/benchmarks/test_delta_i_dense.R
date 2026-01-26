#!/usr/bin/env Rscript
# Delta I computation with Dense BLAS (100 um spacing)
# 960 genes x 98k cells, 50 radii (100-5000 um)

library(sigdiscov)
library(rhdf5)

cat("=== Delta I Computation (Dense BLAS, 100 um spacing) ===\n")
cat("960 genes x 98k cells, 50 radii (100-5000 um)\n\n")

# Load CosMx data
cat("Loading CosMx data...\n")
t0 <- Sys.time()

meta <- read.csv("dataset/cosmx/Lung5_Rep1_meta.csv")
vst <- read.csv("dataset/cosmx/Lung5_Rep1_vst_v2_rpy2.csv", row.names = 1)

load_time <- difftime(Sys.time(), t0, units = "secs")
cat(sprintf("Loaded: %d genes x %d cells (%.1f sec)\n",
            nrow(vst), ncol(vst), load_time))

data <- as.matrix(vst)
coords <- as.matrix(meta[, c("sdimx", "sdimy")])
n_genes <- nrow(data)
n_cells <- ncol(data)
gene_names <- rownames(data)

cat(sprintf("Data: %d genes x %d cells\n", n_genes, n_cells))

# Define radii: 100 um to 5000 um in 100 um steps (50 radii)
radii_um <- seq(100, 5000, by = 100)
radii_mm <- radii_um / 1000
n_radii <- length(radii_um)

cat(sprintf("Radii: %d values from %d to %d um (100 um steps)\n",
            n_radii, min(radii_um), max(radii_um)))
cat(sprintf("Estimated time: ~%.1f hours (50 radii × ~260 sec)\n\n", 50 * 260 / 3600))

# Create HDF5 file
h5_file <- "output/moran_curves_dense.h5"
if (file.exists(h5_file)) file.remove(h5_file)
h5createFile(h5_file)

h5createDataset(h5_file, "moran", dims = c(n_genes, n_genes, n_radii),
                storage.mode = "double", chunk = c(n_genes, n_genes, 1),
                level = 4)

h5write(radii_um, h5_file, "radii_um")
h5write(gene_names, h5_file, "gene_names")
h5write(n_cells, h5_file, "n_cells")

cat(sprintf("HDF5 file: %s\n\n", h5_file))

# Compute Moran matrices using Dense BLAS
cat("--- Computing I(r) curves (Dense BLAS) ---\n")

total_time <- 0

for (i in seq_along(radii_mm)) {
    r <- radii_mm[i]
    t1 <- Sys.time()

    result <- pairwise_moran_dense_cpp(
        data, coords, radius = r, sigma = r / 3,
        chunk_size = 500, verbose = FALSE
    )

    elapsed <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
    total_time <- total_time + elapsed

    # Write to HDF5
    h5write(result$moran, h5_file, "moran", index = list(NULL, NULL, i))

    # Progress
    avg_neighbors <- result$n_edges / n_cells
    eta <- (total_time / i) * (n_radii - i) / 60
    cat(sprintf("  Radius %2d/%d (%4d um): %5.1f sec, %5.0f neighbors/cell, ETA: %5.1f min\n",
                i, n_radii, radii_um[i], elapsed, avg_neighbors, eta))
    flush(stdout())
}

H5close()
cat(sprintf("\nTotal Moran computation: %.1f min (%.2f hours)\n",
            total_time/60, total_time/3600))
cat(sprintf("HDF5 file size: %.1f MB\n", file.size(h5_file) / 1e6))

# Free memory
rm(data, coords, result)
gc()

cat("\n--- Computing Delta I from HDF5 ---\n")
t1 <- Sys.time()

radii_um <- h5read(h5_file, "radii_um")
gene_names <- h5read(h5_file, "gene_names")
n_genes <- length(gene_names)
n_radii <- length(radii_um)

delta_I_matrix <- matrix(NA, n_genes, n_genes)
rownames(delta_I_matrix) <- gene_names
colnames(delta_I_matrix) <- gene_names

argmax_matrix <- matrix(NA, n_genes, n_genes)
I_max_matrix <- matrix(NA, n_genes, n_genes)

for (i in 1:n_genes) {
    row_curves <- matrix(NA, n_genes, n_radii)
    for (r_idx in 1:n_radii) {
        row_curves[, r_idx] <- h5read(h5_file, "moran", index = list(i, NULL, r_idx))
    }
    for (j in 1:n_genes) {
        curve <- row_curves[j, ]
        delta_result <- compute_delta_i(curve, radii_um)
        delta_I_matrix[i, j] <- delta_result$delta_I
        argmax_matrix[i, j] <- delta_result$argmax
        I_max_matrix[i, j] <- delta_result$I_max
    }
    if (i %% 100 == 0) {
        cat(sprintf("  Processed %d/%d genes\n", i, n_genes))
        flush(stdout())
    }
}

H5close()
delta_time <- difftime(Sys.time(), t1, units = "secs")
cat(sprintf("Delta I computation: %.1f sec\n", delta_time))

# Save results
saveRDS(list(
    delta_I = delta_I_matrix,
    argmax = argmax_matrix,
    I_max = I_max_matrix,
    radii_um = radii_um,
    gene_names = gene_names,
    n_cells = n_cells
), "output/delta_I_dense_results.rds")

# Summary
cat("\n=== Results ===\n")
cat(sprintf("Delta I range: [%.4f, %.4f]\n",
            min(delta_I_matrix, na.rm = TRUE),
            max(delta_I_matrix, na.rm = TRUE)))
cat(sprintf("Mean |Delta I|: %.4f\n", mean(abs(delta_I_matrix), na.rm = TRUE)))

cat("\n--- Top 20 gene pairs by |Delta I| ---\n")
delta_vec <- as.vector(delta_I_matrix)
idx_order <- order(abs(delta_vec), decreasing = TRUE)

count <- 0
for (k in idx_order) {
    if (count >= 20) break
    i <- ((k - 1) %% n_genes) + 1
    j <- ((k - 1) %/% n_genes) + 1
    if (i != j) {
        count <- count + 1
        cat(sprintf("%2d. %s - %s: Delta I = %.4f, I_max = %.4f at %d um\n",
                    count, gene_names[i], gene_names[j],
                    delta_I_matrix[i, j], I_max_matrix[i, j],
                    radii_um[argmax_matrix[i, j]]))
    }
}

total_elapsed <- difftime(Sys.time(), t0, units = "mins")
cat(sprintf("\nTotal time: %.1f min (%.2f hours)\n", total_elapsed, as.numeric(total_elapsed)/60))
cat(sprintf("HDF5 file: %s\n", h5_file))
cat("Results saved to: output/delta_I_dense_results.rds\n")
