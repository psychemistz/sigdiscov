# sigdiscov

**Spatial Signature Discovery for Spatial Transcriptomics**

R package for computing spatial correlation metrics in spatial transcriptomics data. Supports both Visium (10x Genomics) and single-cell resolution platforms (CosMx, Xenium, MERFISH).

## Features

- **Visium support**: Binary weights for grid-based spots (~100um spacing)
- **Single-cell support**: Gaussian weights for continuous coordinates
- **Metrics**: Bivariate Moran's I and I_ND (cosine similarity)
- **Delta I**: Signed distance-dependent correlation signatures
- **Batch computation**: All cell type pairs with HDF5 output
- **High performance**: RcppArmadillo matrix operations

## Installation

```r
# Install from GitHub
devtools::install_github("psychemistz/sigdiscov")
```

### Dependencies

**Required:** R (>= 4.0.0), Rcpp, RcppArmadillo, Matrix

**Optional:** Seurat (data loading), rhdf5 (HDF5 I/O)

## Quick Start

### Visium: From VST File

```r
library(sigdiscov)

# Load VST-transformed data
data <- read.table("vst.tsv", header=TRUE, row.names=1)
colnames(data) <- gsub("X", "", colnames(data))
data <- as.matrix(data)

# Parse spot coordinates from column names (format: "ROW_COL")
spot_coords <- parse_spot_names(colnames(data))

# Compute pairwise Moran's I matrix
result <- pairwise_moran(data, spot_coords, max_radius = 3)
result[1:5, 1:5]
```

### Visium: Signature Analysis

```r
library(sigdiscov)

# Load from Space Ranger output
data <- load_data_visium("path/to/spaceranger/outs")

# Compute spatial signature for a factor gene
sig <- compute_signature_visium(
    data,
    factor_gene = "IL1B",
    radii = seq(100, 500, 100),
    metric = "ind",
    mode = "directional"
)

# Top correlated genes
head(sig[order(-sig$ind_r1_val), ], 20)
```

### Single-Cell: CosMx/Xenium/MERFISH

```r
library(sigdiscov)

# Load single-cell data
data <- load_data_cosmx("path/to/cosmx")
# or: data <- load_data_anndata("data.h5ad")

# Compute signature for cell type pair
sig <- compute_signature_sc(
    data,
    factor_gene = "TGFB1",
    sender_celltype = "Macrophage",
    receiver_celltype = "Fibroblast",
    radii = seq(10, 100, 10)
)
```

### Batch Computation (All Cell Type Pairs)

```r
# Compute for all cell type pairs at multiple radii
result <- compute_batch_sc(
    data,
    factor_genes = c("IL1B", "TGFB1", "TNF"),
    radii = seq(10, 100, 10)
)

# Save to HDF5 for Python
save_hdf5_sc(result, "signatures.h5")
```

## Main Functions

### Data Loading

| Function | Description |
|----------|-------------|
| `load_data_visium()` | Load Visium from Space Ranger output |
| `load_data_cosmx()` | Load CosMx data |
| `load_data_anndata()` | Load from AnnData H5AD file |
| `as_data_visium()` | Convert Seurat object |
| `parse_spot_names()` | Parse "ROW_COL" spot names to coordinates |

### Visium Analysis

| Function | Description |
|----------|-------------|
| `compute_signature_visium()` | Main Visium analysis |
| `pairwise_moran()` | Pairwise Moran's I matrix |
| `create_weights_visium()` | Binary weight matrix |
| `create_ring_weights_visium()` | Ring (annular) weights |

### Single-Cell Analysis

| Function | Description |
|----------|-------------|
| `compute_signature_sc()` | Main single-cell analysis |
| `compute_batch_sc()` | Batch all cell type pairs |
| `create_weights_sc()` | Gaussian weight matrix (sigma = radius/3) |
| `save_hdf5_sc()` / `load_hdf5_sc()` | HDF5 I/O |

### Core Metrics

| Function | Description |
|----------|-------------|
| `compute_moran_from_lag()` | Bivariate Moran's I |
| `compute_ind_from_lag()` | I_ND (cosine similarity) |
| `compute_delta_i()` | Signed delta I from curve |
| `batch_permutation_test()` | Significance testing |

## Metrics

### Moran's I
Classic spatial autocorrelation: `I = z_f' * W * z_g / n`

### I_ND (Normalized Directional)
Cosine similarity between factor and spatial lag: `I_ND = z_f' * lag_g / (||z_f|| * ||lag_g||)`

Bounded [-1, 1], interpretable as correlation.

### Delta I
Distance-dependent signature: `delta_I = sign * (I_max - I_min)`

| Pattern | delta_I | Interpretation |
|---------|---------|----------------|
| Decay | > 0 | Paracrine signaling (high near, low far) |
| Increase | < 0 | Avoidance pattern |
| Flat | ~ 0 | Constitutive expression |

## Weight Matrices

### Visium (Binary)
All spots within radius get equal weight (1/n_neighbors).

### Single-Cell (Gaussian)
Gaussian kernel with `sigma = radius / 3`:
```
w_ij = exp(-d_ij^2 / (2 * sigma^2))
```

## License

MIT

## Citation

Beibei Ru, Lanqi Gong, Emily Yang, Seongyong Park, George Zaki, Kenneth Aldape, Lalage Wakefield, Peng Jiang. Inference of secreted protein activities in intercellular communication. [Link](https://github.com/data2intelligence/SecAct)
