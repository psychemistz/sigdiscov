# sigdiscov

**Spatial Correlation Analysis for Spatial Transcriptomics**

[![Version](https://img.shields.io/badge/version-1.2.0-blue.svg)](https://github.com/psychemistz/sigdiscov/releases)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

Compute pairwise Moran's I statistics, signed delta I signatures, and directional spatial correlations (I_ND) for spatial transcriptomics data. Supports both spot-based platforms (10x Visium) and single-cell resolution platforms (CosMx, Xenium, MERFISH).

## Version History

| Version | Features |
|---------|----------|
| **v1.2.0** | Single-cell ST support (CosMx, Xenium, MERFISH), cell-type-aware I_ND, permutation testing |
| **v1.1.0** | Signed delta I signatures, distance-dependent spatial correlation curves |
| **v1.0.0** | Basic pairwise Moran's I for Visium |

## R Package Installation

```r
# Install latest version from GitHub
devtools::install_github("psychemistz/sigdiscov", subdir = "sigdiscov")

# Install specific version
devtools::install_github("psychemistz/sigdiscov@v1.2.0", subdir = "sigdiscov")
```

## Key Features

### 1. Pairwise Moran's I (All versions)
Compute spatial autocorrelation between all gene pairs.

```r
library(sigdiscov)

# Load and process Visium data
visium <- load_visium_data("path/to/spaceranger_output")
visium <- filter_in_tissue(visium)
data_vst <- vst_transform(visium$counts)
coords <- get_spot_coords(visium)

# Compute pairwise Moran's I
moran <- pairwise_moran(data_vst, coords)
```

### 2. Signed Delta I Signatures (v1.1.0+)
Distance-dependent spatial correlation curves to distinguish responsive genes.

```r
# Compute signed delta I for a factor gene
result <- compute_signed_delta_I(
    expr_matrix = data_vst,
    coords = coords,
    factor_gene = "IL1B",
    radii = seq(50, 300, 50)
)

# Plot the Moran curve
plot_moran_curve(result, gene = "CXCL8")
```

### 3. Single-Cell ST Analysis (v1.2.0)
Cell-type-aware directional I_ND for single-cell resolution platforms.

```r
# Load single-cell ST data (CosMx, Xenium, MERFISH)
cell_meta <- data.frame(
    cell_id = colnames(expr_matrix),
    x = spatial_coords$x,
    y = spatial_coords$y,
    cell_type = cell_annotations
)

# Compute I_ND for a specific cell-type triplet
result <- compute_IND_sc(
    expr_matrix = expr_matrix,
    cell_meta = cell_meta,
    factor_gene = "IL1B",
    sender_type = "Macrophage",
    receiver_type = "Fibroblast",
    radii = seq(20, 100, 20)
)

# Multi-factor analysis (efficient)
results <- compute_IND_multi_factor(
    expr_matrix = expr_matrix,
    cell_meta = cell_meta,
    factor_genes = c("IL1B", "TGFB1", "IFNG"),
    sender_type = "Macrophage",
    receiver_type = "Fibroblast"
)
```

### 4. Permutation Testing (v1.2.0)
Statistical significance testing for spatial correlations.

```r
# Permutation test for single factor
perm_result <- permutation_test_spatial(
    expr_matrix = expr_matrix,
    cell_meta = cell_meta,
    factor_gene = "IL1B",
    sender_type = "Macrophage",
    receiver_type = "Fibroblast",
    n_perm = 999
)

# Vectorized multi-factor test (9x faster)
results <- permutation_test_multi_factor(
    expr_matrix = expr_matrix,
    cell_meta = cell_meta,
    factor_genes = c("IL1B", "TGFB1", "IFNG", "TNF", "IL6"),
    sender_type = "Macrophage",
    receiver_type = "Fibroblast",
    n_perm = 999
)

# Significant genes (FDR < 0.05)
sig_genes <- results[["IL1B"]][results[["IL1B"]]$p_adj < 0.05, ]
```

## Standalone C++ Tool

For command-line usage without R:

```bash
# Build
cd accelerate/Release
make

# Run
./pairwise_moran_I -i input.tsv -o output.txt
```

### Options

```
-i    Input file (tab-separated, genes x spots)
-o    Output file
-r    Maximum grid radius (default: 5)
-p    Platform: 0=Visium, 1=Old ST (default: 0)
```

## Performance

- **BLAS backend**: Apple Accelerate (Mac), OpenBLAS/MKL (Linux/Windows)
- **Pairwise Moran's I**: ~2,000x speedup vs naive implementation
- **Single-cell I_ND**: Optimized for datasets with 100k+ cells
- **Permutation testing**: Vectorized BLAS operations for multi-factor tests

## Repository Structure

```
sigdiscov/
├── sigdiscov/           # R package
│   ├── R/
│   │   ├── moran.R           # Core Moran's I functions
│   │   ├── signed_delta_I.R  # Signed delta I signatures
│   │   ├── singlecell.R      # Single-cell ST functions
│   │   ├── permutation.R     # Permutation testing
│   │   └── visium.R          # Visium data loading
│   └── src/             # C++ source (RcppArmadillo)
├── accelerate/          # Standalone C++ tool
└── dataset/             # Example data
```

## Authors

- Peng Jiang (peng.jiang@nih.gov)
- Seongyong Park (seongyong.park@nih.gov)

## License

MIT License
