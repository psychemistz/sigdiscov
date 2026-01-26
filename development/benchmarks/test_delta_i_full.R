#!/usr/bin/env Rscript
# Full Delta I computation using Dense BLAS
# 960 genes x 98k cells, 250 radii (20-5000 um)

library(sigdiscov)

cat("=== Full Delta I Computation (Dense BLAS) ===\n")
cat("960 genes x 98k cells, 250 radii (20-5000 um)\n\n")

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

cat(sprintf("Data: %d genes x %d cells\n", n_genes, n_cells))

# Define radii: 20 um to 5000 um in 20 um steps
# CosMx coords are in mm, so convert
radii_um <- seq(20, 5000, by = 20)
radii_mm <- radii_um / 1000

cat(sprintf("Radii: %d values from %d to %d um\n\n",
            length(radii_um), min(radii_um), max(radii_um)))

# Compute pairwise Moran's I at each radius using Dense BLAS
cat("--- Computing I(r) curves (Dense BLAS) ---\n")

I_curves <- vector("list", length(radii_mm))
total_time <- 0

for (i in seq_along(radii_mm)) {
    r <- radii_mm[i]

    t1 <- Sys.time()
    result <- pairwise_moran_dense_cpp(
        data,
        coords,
        radius = r,
        sigma = r / 3,
        chunk_size = 500,
        verbose = FALSE
    )
    elapsed <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
    total_time <- total_time + elapsed

    I_curves[[i]] <- result$moran

    # Progress every 10 radii
    if (i %% 10 == 0 || i == 1 || i == length(radii_mm)) {
        avg_neighbors <- result$n_edges / n_cells
        eta <- (total_time / i) * (length(radii_mm) - i) / 60
        cat(sprintf("  Radius %3d/%d (%4d um): %5.1f sec, %.0f neighbors/cell, ETA: %.1f min\n",
                    i, length(radii_mm), radii_um[i], elapsed, avg_neighbors, eta))
    }
}

cat(sprintf("\nTotal Moran computation: %.1f min\n", total_time/60))

# Build I(r) array
cat("\n--- Building I(r) array ---\n")
I_array <- array(NA, dim = c(n_genes, n_genes, length(radii_mm)),
                  dimnames = list(rownames(data), rownames(data),
                                  paste0("r", radii_um)))

for (i in seq_along(I_curves)) {
    I_array[,,i] <- I_curves[[i]]
}

# Free memory
rm(I_curves)
gc()

# Compute delta I for each gene pair
cat("--- Computing Delta I matrix ---\n")
t1 <- Sys.time()

delta_I_matrix <- matrix(NA, n_genes, n_genes)
rownames(delta_I_matrix) <- rownames(data)
colnames(delta_I_matrix) <- rownames(data)

argmax_matrix <- matrix(NA, n_genes, n_genes)
I_max_matrix <- matrix(NA, n_genes, n_genes)

for (i in 1:n_genes) {
    for (j in 1:n_genes) {
        curve <- I_array[i, j, ]
        delta_result <- compute_delta_i(curve, radii_um)
        delta_I_matrix[i, j] <- delta_result$delta_I
        argmax_matrix[i, j] <- delta_result$argmax
        I_max_matrix[i, j] <- delta_result$I_max
    }
    if (i %% 100 == 0) {
        cat(sprintf("  Processed %d/%d genes\n", i, n_genes))
    }
}

delta_time <- difftime(Sys.time(), t1, units = "secs")
cat(sprintf("Delta I computation: %.1f sec\n", delta_time))

# Save results
cat("\n--- Saving results ---\n")
saveRDS(list(
    delta_I = delta_I_matrix,
    argmax = argmax_matrix,
    I_max = I_max_matrix,
    radii_um = radii_um,
    gene_names = rownames(data),
    n_cells = n_cells
), "output/delta_I_full_results.rds")

# Summary
cat("\n=== Results ===\n")
cat(sprintf("Delta I range: [%.4f, %.4f]\n",
            min(delta_I_matrix, na.rm = TRUE),
            max(delta_I_matrix, na.rm = TRUE)))
cat(sprintf("Mean |Delta I|: %.4f\n", mean(abs(delta_I_matrix), na.rm = TRUE)))

# Top gene pairs by |Delta I|
cat("\n--- Top 20 gene pairs by |Delta I| ---\n")
delta_vec <- as.vector(delta_I_matrix)
idx_order <- order(abs(delta_vec), decreasing = TRUE)

count <- 0
for (k in idx_order) {
    if (count >= 20) break
    i <- ((k - 1) %% n_genes) + 1
    j <- ((k - 1) %/% n_genes) + 1
    if (i != j) {  # Skip diagonal
        count <- count + 1
        cat(sprintf("%2d. %s - %s: Delta I = %.4f, I_max = %.4f at %d um\n",
                    count,
                    rownames(data)[i], rownames(data)[j],
                    delta_I_matrix[i, j],
                    I_max_matrix[i, j],
                    radii_um[argmax_matrix[i, j]]))
    }
}

total_elapsed <- difftime(Sys.time(), t0, units = "mins")
cat(sprintf("\nTotal time: %.1f min\n", total_elapsed))
cat("Results saved to: output/delta_I_full_results.rds\n")
