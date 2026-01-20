# sigdiscov

**Spatial Correlation Analysis for Spatial Transcriptomics**

[![Version](https://img.shields.io/badge/version-1.1.0-blue.svg)](https://github.com/psychemistz/sigdiscov/releases/tag/v1.1.0)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

Compute pairwise Moran's I statistics and signed delta I signatures for spatial transcriptomics data. Optimized for 10x Visium data using BLAS matrix operations.

## What's in This Version (v1.1.0)

- **Pairwise Moran's I**: Spatial autocorrelation between all gene pairs
- **Signed Delta I Signatures**: Distance-dependent spatial correlation curves
- **Bivariate and Directional Modes**: Support for both symmetric and I_ND correlation
- **BLAS Optimization**: High-performance matrix operations via RcppArmadillo
- **10x Visium Support**: Designed for spot-based spatial transcriptomics

> **Note**: For basic Moran's I only, see [v1.0.0](https://github.com/psychemistz/sigdiscov/tree/release/v1.0.0). For single-cell ST support and permutation testing, upgrade to [v1.2.0](https://github.com/psychemistz/sigdiscov/tree/release/v1.2.0).

## Installation

```r
# Install this version
devtools::install_github("psychemistz/sigdiscov@v1.1.0")

# Or install latest
devtools::install_github("psychemistz/sigdiscov")
```

## Quick Start

```r
library(sigdiscov)

# Full pipeline from Space Ranger output
result <- run_pipeline("path/to/spaceranger_output")

# Or step-by-step:
visium <- load_visium_data("path/to/spaceranger_output")
visium <- filter_in_tissue(visium)
data_vst <- vst_transform(visium$counts)
coords <- get_spot_coords(visium)

# Basic Moran's I
moran <- pairwise_moran(data_vst, coords)

# Signed delta I signatures
delta_I <- compute_signed_delta_I(data_vst, coords, target_genes,
                                   mode = "bivariate")
```

## Signed Delta I Method

The signed delta I method computes distance-dependent spatial correlation curves to distinguish true responsive genes from constitutively expressed genes:

- **Bivariate mode**: Symmetric correlation between gene pairs
- **Directional mode (I_ND)**: Asymmetric correlation for ligand-receptor analysis

## Performance

- **BLAS backend**: Apple Accelerate (Mac), OpenBLAS/MKL (Linux/Windows)
- **Typical speedup**: ~2,000x compared to naive implementation

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
