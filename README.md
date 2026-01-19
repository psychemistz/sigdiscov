# sigdiscov

**Pairwise Moran's I for Spatial Transcriptomics**

An R package for computing pairwise Moran's I statistics between genes in spatial transcriptomics data. Uses optimized BLAS matrix operations for high-performance computation.

## Features

- Fast pairwise Moran's I computation using matrix multiplication (BLAS)
- Support for 10x Visium spatial transcriptomics data
- Built-in VST normalization via sctransform
- Sparse matrix support
- Complete pipeline from raw Space Ranger output to Moran's I results

## Installation

### From GitHub

```r
# Install devtools if not already installed
install.packages("devtools")

# Install sigdiscov
devtools::install_github("psychemistz/sigdiscov")
```

### Dependencies

Required:
- R (>= 3.5.0)
- Rcpp (>= 1.0.0)
- RcppArmadillo
- Matrix

Optional (for VST normalization):
- sctransform

```r
# Install optional dependency for VST
install.packages("sctransform")
```

## Quick Start

### Full Pipeline (Recommended)

Run the complete analysis from raw Space Ranger output:

```r
library(sigdiscov)

# Run full pipeline
result <- run_pipeline("path/to/spaceranger_output")

# Access results
moran_matrix <- result$moran
gene_names <- result$gene_names
```

### Step-by-Step Analysis

For more control over each step:

```r
library(sigdiscov)

# Step 1: Load Visium data
visium <- load_visium_data("path/to/spaceranger_output")
visium <- filter_in_tissue(visium)

# Step 2: VST normalization (requires sctransform)
data_vst <- vst_transform(visium$counts, min_cells = 5)

# Step 3: Get spot coordinates
spot_coords <- get_spot_coords(visium)

# Step 4: Compute pairwise Moran's I
result <- pairwise_moran(
  data_vst,
  spot_coords,
  max_radius = 5,
  platform = "visium"
)
```

### Using Pre-processed Data

If you already have normalized expression data:

```r
library(sigdiscov)

# Load your data (genes x spots matrix)
data <- read.table("expression.tsv", header = TRUE, row.names = 1)
data <- as.matrix(data)

# Parse spot coordinates from column names (format: "ROWxCOL")
spot_coords <- parse_spot_names(colnames(data))

# Compute pairwise Moran's I
result <- pairwise_moran(data, spot_coords)
```

## Functions

### Data Loading

| Function | Description |
|----------|-------------|
| `load_visium_data()` | Load Visium data from Space Ranger output |
| `filter_in_tissue()` | Filter to in-tissue spots only |
| `get_spot_coords()` | Extract spot coordinates for Moran's I |

### Preprocessing

| Function | Description |
|----------|-------------|
| `vst_transform()` | Apply VST normalization using sctransform |
| `run_pipeline()` | Complete pipeline: load, preprocess, compute |

### Moran's I Computation

| Function | Description |
|----------|-------------|
| `pairwise_moran()` | Compute pairwise Moran's I between all genes |
| `moran_I()` | Compute Moran's I for a single gene pair |
| `create_weight_matrix()` | Create spatial weight matrix |
| `parse_spot_names()` | Parse "ROWxCOL" format spot names |

### I/O

| Function | Description |
|----------|-------------|
| `save_moran_result()` | Save results in lower triangular format |
| `load_moran_result()` | Load results from file |

## Parameters

### pairwise_moran()

| Parameter | Description | Default |
|-----------|-------------|---------|
| `data` | Gene expression matrix (genes x spots) | required |
| `spot_coords` | Spot coordinates (row, col) | required |
| `max_radius` | Maximum grid radius for spatial weights | 5 |
| `platform` | Platform type: "visium" or "old" | "visium" |
| `same_spot` | Include same-spot weights | TRUE |
| `mode` | "paired", "first", or "single" | "paired" |
| `verbose` | Print progress messages | TRUE |

## Performance

The package uses BLAS matrix operations to compute pairwise Moran's I efficiently:

- **Algorithm**: Reformulated as `Result = X @ W' @ X^T / sum(W')`
- **Typical speedup**: ~2,000x compared to naive implementation
- **Memory**: O(n_genes^2) for result matrix

## Citation

If you use this package, please cite:

```
Jiang P, Park S. sigdiscov: Pairwise Moran's I for Spatial Transcriptomics.
```

## Authors

- Peng Jiang (pengj@alumni.princeton.edu)
- Seongyong Park

## License

MIT License
