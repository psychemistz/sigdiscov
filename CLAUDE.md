# sigdiscov - Spatial Signature Discovery Package

## Project Overview

**sigdiscov** (v1.3.0) is a unified R package for spatial signature discovery in spatial transcriptomics data. It supports both Visium (10x Genomics) and single-cell resolution platforms.

**Repository**: GitHub
**Language**: R with RcppArmadillo C++ backend
**License**: MIT

## Core Functionality

### 1. Visium Analysis (`R/signature_visium.R`)
- `pairwise_moran()` - Compute pairwise Moran's I matrix across genes
- `moran_I()` - Single gene pair Moran's I
- `create_weight_matrix()` - Gaussian distance decay weight matrix

### 2. Single-Cell Analysis (`R/signature_sc.R`)
- `compute_signature_sc()` - Main entry point for single-cell spatial analysis
- `compute_batch_sc()` - Batch computation across cell type pairs
- `create_weights_sc()` - Weight matrix for single-cell data

### 3. Spatial Metrics (`R/metrics.R`)
- `compute_moran_from_lag()` - Bivariate Moran's I from pre-computed lag
- `compute_ind_from_lag()` - I_ND (cosine similarity) metric
- `compute_metric_batch()` - Batch metric computation

### 4. Statistical Testing (`R/permutation.R`)
- `permutation_test_core()` - Single gene permutation test
- `batch_permutation_test()` - Genome-wide significance testing

## Key Technical Details

### Visium Hexagonal Grid Geometry
- Array coordinates differ by 2 between adjacent spots
- Hex shift factor: `0.5 * sqrt(3)` for row spacing
- Default scale: VISIUM_DISTANCE = 100

### Weight Matrix (Critical)
```cpp
// Gaussian distance decay with sigma=100
inline double distance_decay(double d) {
    return std::exp(-(d / 100.0) * (d / 100.0) / 2.0);
}
```

**Validation**: For test dataset with 3813 spots, weight sum should be ~22079

### Bivariate Moran's I Formula
```
I = X * W * X^T / weight_sum
```
Where W is the Gaussian-weighted spatial weight matrix (NOT row-normalized for pairwise_moran).

## File Structure

```
sigdiscov/
├── R/
│   ├── signature_visium.R   # Visium analysis (pairwise_moran)
│   ├── signature_sc.R       # Single-cell analysis
│   ├── metrics.R            # Spatial correlation metrics
│   ├── weights.R            # Weight matrix creation
│   ├── permutation.R        # Statistical testing
│   ├── data_loaders.R       # Data loading functions
│   ├── utils.R              # Utilities (parse_spot_names, standardize)
│   └── RcppExports.R        # Auto-generated C++ bindings
├── src/
│   ├── moran_pairwise.cpp   # Visium pairwise Moran (Gaussian weights)
│   ├── ind_matrix.cpp       # I_ND matrix computation
│   ├── metrics.cpp          # Core metric functions
│   ├── weights.cpp          # Weight matrix C++ functions
│   ├── permutation.cpp      # Permutation test C++
│   └── RcppExports.cpp      # Auto-generated
├── tests/
│   └── testthat/            # Test files (177 tests)
├── man/                     # Documentation (auto-generated)
├── DESCRIPTION
├── NAMESPACE
└── .Rbuildignore
```

## Quick Start Examples

### Visium VST File Workflow
```r
library(sigdiscov)

# Load preprocessed VST data
data <- read.table("3_vst.tsv", header=TRUE, row.names=1)
colnames(data) <- gsub("X", "", colnames(data))
data <- as.matrix(data)

# Parse spot coordinates from "ROWxCOL" format
spot_coords <- parse_spot_names(colnames(data))

# Compute pairwise Moran's I matrix
result <- pairwise_moran(data, spot_coords, max_radius = 3)
```

### Single-Cell Analysis
```r
# Load single-cell data
sc_data <- load_data_anndata("data.h5ad")

# Create weight matrix
W <- create_weights_sc(sc_data$coords, radius = 100)

# Compute signatures
sig <- compute_signature_sc(sc_data, W, factor_genes = factor_list)
```

## Recent Changes (v1.3.0)

1. **Unified Package**: Merged sigdiscov.visium and sigdiscov.sc into single package
2. **Restored Gaussian Weights**: Fixed pairwise_moran to use original Gaussian distance decay (was incorrectly using binary weights)
3. **Added parse_spot_names()**: For "ROWxCOL" format parsing
4. **Updated defaults**: `same_spot = FALSE` in pairwise_moran

## Common Issues

### Weight Sum = 0
**Cause**: Incorrect weight matrix calculation (binary vs Gaussian)
**Solution**: Ensure using Gaussian distance decay weights, not binary

### Package Build Too Large
**Cause**: Dataset directories included
**Solution**: .Rbuildignore excludes dataset/, accelerate/, development files

## Build Commands

```bash
# Install from source
R CMD INSTALL .

# Build and check
R CMD build .
R CMD check sigdiscov_*.tar.gz

# Run tests
Rscript -e "devtools::test()"

# Generate documentation
Rscript -e "devtools::document()"
```

## Citation

Park, S. et al. "SecAct: Automated extraction of cell type-specific signatures from single-cell transcriptomics data." (2025)
