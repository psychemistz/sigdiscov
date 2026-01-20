# sigdiscov

**Spatial Correlation Analysis for Spatial Transcriptomics**

[![Version](https://img.shields.io/badge/version-1.2.0-blue.svg)](https://github.com/psychemistz/sigdiscov/releases/tag/v1.2.0)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

Compute pairwise Moran's I statistics, signed delta I signatures, and perform permutation testing for spatial transcriptomics data. Supports both spot-based (Visium) and single-cell resolution (CosMx, Xenium, MERFISH) platforms.

## What's in This Version (v1.2.0)

- **Pairwise Moran's I**: Spatial autocorrelation between all gene pairs
- **Signed Delta I Signatures**: Distance-dependent spatial correlation curves
- **Single-Cell ST Support**: Cell-type-aware analysis for CosMx, Xenium, MERFISH
- **Permutation Testing**: Statistical significance testing with vectorized multi-factor permutations
- **BLAS Optimization**: High-performance matrix operations via RcppArmadillo

> **Note**: For basic Moran's I only, see [v1.0.0](https://github.com/psychemistz/sigdiscov/tree/release/v1.0.0). For signed delta I without single-cell support, see [v1.1.0](https://github.com/psychemistz/sigdiscov/tree/release/v1.1.0).

## Installation

```r
# Install this version
devtools::install_github("psychemistz/sigdiscov@v1.2.0")

# Or install latest
devtools::install_github("psychemistz/sigdiscov")
```

## Quick Start

### Visium (Spot-based)

```r
library(sigdiscov)

# Load and process Visium data
visium <- load_visium_data("path/to/spaceranger_output")
visium <- filter_in_tissue(visium)
data_vst <- vst_transform(visium$counts)
coords <- get_spot_coords(visium)

# Basic Moran's I
moran <- pairwise_moran(data_vst, coords)

# Signed delta I signatures
delta_I <- compute_signed_delta_I(data_vst, coords, factor_gene = "IL1B")

# Permutation testing
perm_result <- permutation_test_spatial(data_vst, coords, factor_gene = "IL1B")
```

### Single-Cell ST (CosMx, Xenium, MERFISH)

```r
# Cell-type-aware I_ND analysis
result <- compute_IND_sc(
  expr_matrix = expr_data,
  cell_coords = coords,
  cell_types = cell_type_labels,
  factor_genes = c("IL1B", "TGFB1"),
  sender_type = "Macrophage",
  receiver_type = "Fibroblast"
)

# Get I_ND curves
curves <- get_IND_curve_sc(
  expr_matrix = expr_data,
  cell_coords = coords,
  cell_types = cell_type_labels,
  factor_gene = "IL1B",
  target_gene = "COL1A1",
  sender_type = "Macrophage",
  receiver_type = "Fibroblast"
)

# Permutation test with cell type structure
perm_result <- permutation_test_IND_sc(
  expr_matrix = expr_data,
  cell_coords = coords,
  cell_types = cell_type_labels,
  factor_gene = "IL1B",
  sender_type = "Macrophage",
  receiver_type = "Fibroblast",
  n_perm = 1000
)
```

## Key Features

### Signed Delta I Method

Computes distance-dependent spatial correlation curves to distinguish true responsive genes from constitutively expressed genes:

- **Bivariate mode**: Symmetric correlation between gene pairs
- **Directional mode (I_ND)**: Asymmetric correlation for ligand-receptor analysis
- **Ring vs Circular weights**: Distance-specific or cumulative spatial structure

### Single-Cell Resolution Analysis

Cell-type-aware spatial correlation:

- Sender-receiver cell type specification
- Distance-based weight matrices accounting for cell positions
- I_ND curves showing cell-type-specific spatial relationships

### Permutation Testing

Statistical significance assessment:

- Vectorized multi-factor permutations for efficiency
- Cell-type-aware permutations preserving spatial structure
- p-value and FDR computation

## Performance

- **BLAS backend**: Apple Accelerate (Mac), OpenBLAS/MKL (Linux/Windows)
- **Permutation speedup**: ~9.7x with vectorized multi-factor approach
- **Typical runtime**: 19,729 genes x 3,813 spots in ~2 minutes

## Version History

| Version | Features |
|---------|----------|
| **v1.2.0** | Single-cell ST support, permutation testing (this version) |
| **v1.1.0** | Signed delta I signatures |
| **v1.0.0** | Basic pairwise Moran's I |

## Citation

Beibei Ru, Lanqi Gong, Emily Yang, Seongyong Park, George Zaki, Kenneth Aldape, Lalage Wakefield, Peng Jiang. Inference of secreted protein activities in intercellular communication. [SecAct](https://github.com/data2intelligence/SecAct)

## Author

Seongyong Park (seongyong.park@nih.gov)

## License

MIT License
