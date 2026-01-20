# sigdiscov

**Pairwise Moran's I for Spatial Transcriptomics**

[![Version](https://img.shields.io/badge/version-1.0.0-blue.svg)](https://github.com/psychemistz/sigdiscov/releases/tag/v1.0.0)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

Compute pairwise Moran's I statistics between genes in spatial transcriptomics data. Optimized for 10x Visium data using BLAS matrix operations.

## What's in This Version (v1.0.0)

- **Pairwise Moran's I**: Spatial autocorrelation between all gene pairs
- **BLAS Optimization**: High-performance matrix operations
- **10x Visium Support**: Designed for spot-based spatial transcriptomics

> **Note**: For signed delta I signatures, upgrade to [v1.1.0](https://github.com/psychemistz/sigdiscov/tree/release/v1.1.0). For single-cell ST support, upgrade to [v1.2.0](https://github.com/psychemistz/sigdiscov/tree/release/v1.2.0).

## R Package Installation

```r
# Install this version
devtools::install_github("psychemistz/sigdiscov@v1.0.0", subdir = "sigdiscov")

# Or install latest
devtools::install_github("psychemistz/sigdiscov", subdir = "sigdiscov")
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
moran <- pairwise_moran(data_vst, coords)
```

## Standalone C++ Tool

For command-line usage without R:

```bash
# Build
cd accelerate/Release
make all

# Run
./pairwise_moran_I -i input.tsv -o output.txt
```

### Options

```
-i    Input file (tab-separated, genes x spots)
-o    Output file
-r    Maximum grid radius (default: 5)
-p    Platform: 0=Visium, 1=Old ST (default: 0)
-b    Paired genes: 1=yes, 0=no (default: 1)
-g    All genes: 1=yes, 0=no (default: 1)
-s    Same spot: 1=yes, 0=no (default: 1)
```

## Performance

- **BLAS backend**: Apple Accelerate (Mac), OpenBLAS/MKL (Linux/Windows)
- **Typical speedup**: ~2,000x compared to naive implementation
- **Example**: 19,729 genes x 3,813 spots in ~15 seconds

## Version History

| Version | Features |
|---------|----------|
| **v1.2.0** | Single-cell ST support, permutation testing |
| **v1.1.0** | Signed delta I signatures |
| **v1.0.0** | Basic pairwise Moran's I (this version) |

## Authors

- Peng Jiang (peng.jiang@nih.gov)
- Seongyong Park (seongyong.park@nih.gov)

## License

MIT License
