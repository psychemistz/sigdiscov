# sigdiscov

**Spatial Signature Discovery for Spatial Transcriptomics**

R package for computing spatial correlation metrics in spatial transcriptomics data. Supports both Visium (10x Genomics) and single-cell resolution platforms (CosMx, Xenium, MERFISH).

## Features

- **Visium support**: Circular RBF weights (SpaCET-compatible) as default
- **Single-cell support**: Gaussian weights for continuous coordinates
- **Metrics**: Bivariate Moran's I and I_ND (cosine similarity)
- **Delta I**: Signed distance-dependent correlation signatures
- **Batch computation**: All cell type pairs with HDF5 output
- **High performance**: RcppArmadillo matrix operations
- **SpaCET-compatible**: Uses population SD (N) for z-normalization

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

# Parse spot coordinates from column names (format: "ROWxCOL")
spot_coords <- parse_spot_names(colnames(data))

# Convert to physical coordinates (um)
coords <- data.frame(
  x = spot_coords$col * 100,
  y = spot_coords$row * 100 * sqrt(3) / 2
)

# Compute pairwise Moran's I matrix (circular weights, default)
result <- pairwise_moran(data, coords, radius = 200, sigma = 100)
result$moran[1:5, 1:5]

# Or use legacy grid-based weights
result_grid <- pairwise_moran(data, spot_coords, weight_type = "grid", max_radius = 3)
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

## Command-Line Analysis

A unified analysis script is provided for comparing spatial metrics with reference signatures:

```bash
# Basic usage (computes Moran's I, I_ND, compares with CytoSig/SecAct)
Rscript run_analysis.R --vst data.tsv --output results

# With permutation testing for z-scores (1000 permutations)
Rscript run_analysis.R --vst data.tsv --output results --n_perm 1000

# Analyze specific factors only
Rscript run_analysis.R --vst data.tsv --output results --factors IFNG,TGFB1

# Custom parameters
Rscript run_analysis.R --vst data.tsv --output results \
    --radius 200 --sigma 100 --n_perm 1000 --seed 42
```

**Options:**
| Argument | Description | Default |
|----------|-------------|---------|
| `--vst` | VST expression matrix (genes × spots) | Required |
| `--output` | Output directory | `output` |
| `--radius` | Spatial radius for weights | 200 |
| `--sigma` | RBF kernel bandwidth | 100 |
| `--n_perm` | Permutations for z-scores (0 to skip) | 1000 |
| `--factors` | Comma-separated factor list | All common |

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
| `pairwise_moran()` | Pairwise Moran's I matrix (circular weights default) |
| `create_weights_visium()` | Weight matrix (circular default, grid optional) |
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

### Visium (Circular RBF - Default)
Circular Euclidean distance with RBF kernel (SpaCET-compatible):
```
w_ij = exp(-d_ij^2 / (2 * sigma^2))  if d_ij <= radius, else 0
```
Default: `radius = 200`, `sigma = 100`

### Visium (Grid-based - Legacy)
Hexagonal grid-based Gaussian distance decay:
```
w_ij = exp(-d_ij^2 / (2 * sigma^2))    # sigma = 100
```
Use with `weight_type = "grid"` and `max_radius` parameter.

### Single-Cell (Gaussian)
Gaussian kernel with `sigma = radius / 3`:
```
w_ij = exp(-d_ij^2 / (2 * sigma^2))
```

## Benchmark: Grid-based Dense vs Sparse Weight Matrix

Tested on real Visium datasets with `weight_type = "grid"` and `max_radius = 3`:

| Dataset | Genes | Spots | Dense W Time | Sparse W Time | Dense W Memory | Sparse W Memory |
|---------|-------|-------|--------------|---------------|----------------|-----------------|
| 1_vst.tsv | 19,729 | 3,813 | 11.76s | 14.43s | 111 MB | 1.5 MB |
| 2_vst.tsv | 19,419 | 4,727 | 15.23s | 15.89s | 171 MB | 1.8 MB |
| 3_vst.tsv | 15,931 | 2,518 | 6.09s | 5.67s | 48 MB | 0.9 MB |

**Key findings:**
- **Identical outputs** (max difference < 1e-16)
- **~99% memory reduction** for weight matrix storage
- **Similar computation speed** (within 5-20%)

## Benchmark: Circular RBF vs SpaCET

Tested on Visium dataset (5,000 genes × 3,813 spots, radius=200, sigma=100):

| Implementation | Time | Speedup | W Memory |
|----------------|------|---------|----------|
| SpaCET R | 536.9s | 1x | 110.9 MB |
| sigdiscov Dense C++ | 25.7s | **20.9x** | 110.9 MB |
| sigdiscov Sparse C++ | 20.9s | **25.7x** | **0.35 MB** |

**Key findings:**
- **Identical outputs** (max difference < 1e-12, correlation = 1.0)
- **21-26x faster** than SpaCET R implementation
- **314x memory reduction** for weight matrix with sparse mode

## SpaCET Compatibility

The default circular weights implementation produces **identical results** to SpaCET:
- Weight matrix difference: < 1e-16
- Moran's I difference: < 1e-12
- Uses population SD (N) for z-normalization

## License

MIT

## Citation

Beibei Ru, Lanqi Gong, Emily Yang, Seongyong Park, George Zaki, Kenneth Aldape, Lalage Wakefield, Peng Jiang. Inference of secreted protein activities in intercellular communication. [Link](https://github.com/data2intelligence/SecAct)
