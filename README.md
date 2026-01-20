# sigdiscov

**Spatial Correlation Analysis for Spatial Transcriptomics**

[![Version](https://img.shields.io/badge/version-1.1.0-blue.svg)](https://github.com/psychemistz/sigdiscov/releases/tag/v1.1.0)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

Compute pairwise Moran's I statistics and signed delta I signatures for spatial transcriptomics data. Optimized for 10x Visium data using BLAS matrix operations.

## What's in This Version (v1.1.0)

- **Pairwise Moran's I**: Spatial autocorrelation between gene pairs
- **Signed Delta I Signatures**: Distance-dependent spatial correlation curves
- **Moran Curve Visualization**: Plot spatial decay patterns

> **Note**: For single-cell ST support (CosMx, Xenium, MERFISH) and permutation testing, upgrade to [v1.2.0](https://github.com/psychemistz/sigdiscov/tree/release/v1.2.0).

## Installation

```r
# Install this version
devtools::install_github("psychemistz/sigdiscov@v1.1.0")

# Or install latest
devtools::install_github("psychemistz/sigdiscov")
```

## Key Features

### 1. Pairwise Moran's I
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

### 2. Signed Delta I Signatures
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

## Performance

- **BLAS backend**: Apple Accelerate (Mac), OpenBLAS/MKL (Linux/Windows)
- **Pairwise Moran's I**: ~2,000x speedup vs naive implementation

## Version History

| Version | Features |
|---------|----------|
| **v1.2.0** | Single-cell ST support, permutation testing |
| **v1.1.0** | Signed delta I signatures (this version) |
| **v1.0.0** | Basic pairwise Moran's I |

## Citation

Beibei Ru, Lanqi Gong, Emily Yang, Seongyong Park, George Zaki, Kenneth Aldape, Lalage Wakefield, Peng Jiang. Inference of secreted protein activities in intercellular communication. [SecAct](https://github.com/data2intelligence/SecAct)

## Author

Seongyong Park (seongyong.park@nih.gov)

## License

MIT License
