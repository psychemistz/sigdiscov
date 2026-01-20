# sigdiscov

**Spatial Signature Discovery for Spatial Transcriptomics**

An R package for computing spatial correlation metrics between genes in spatial transcriptomics data, including pairwise Moran's I and signed delta I signatures. Uses optimized BLAS matrix operations via RcppArmadillo for high-performance computation.

## Features

- Fast pairwise Moran's I computation using matrix multiplication (BLAS)
- **Signed Delta I signatures** for identifying responsive genes vs. constitutive expression
- **Four delta I types**: ring/circular weights x Moran's I/I_ND correlation
- Support for 10x Visium spatial transcriptomics data
- Built-in VST normalization via sctransform
- Sparse matrix support for memory efficiency
- Complete pipeline from raw Space Ranger output to results

## Installation

### From GitHub

```r
# Install devtools if not already installed
install.packages("devtools")

# Install sigdiscov
devtools::install_github("psychemistz/sigdiscov")
```

### Dependencies

**Required:**
- R (>= 3.5.0)
- Rcpp (>= 1.0.0)
- RcppArmadillo
- Matrix

**Optional (for VST normalization):**
- sctransform

```r
# Install optional dependency for VST
install.packages("sctransform")
```

## Quick Start

### Option 1: Full Pipeline (Recommended)

Run the complete analysis from raw Space Ranger output:

```r
library(sigdiscov)

# Run full pipeline
result <- run_pipeline("path/to/spaceranger_output")

# Access results
moran_matrix <- result$moran
gene_names <- result$gene_names
```

### Option 2: Step-by-Step Analysis

For more control over each step:

```r
library(sigdiscov)

# Step 1: Load Visium data from Space Ranger output
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
  platform = "visium",
  same_spot = FALSE
)

# Save results
save_moran_result(result, "output.tsv")
```

## Signed Delta I Signatures

Delta I captures distance-dependent spatial correlation patterns. A gene pair showing high correlation at short distances but low correlation at long distances (decay pattern) indicates a true responsive relationship, while flat patterns suggest constitutive expression.

### Four Delta I Types

| Delta Type | Weight | Correlation | Description |
|------------|--------|-------------|-------------|
| `delta_i_ring_moran` | Ring | Moran's I | Neighbors in distance band [r_inner, r_outer) |
| `delta_i_cir_moran` | Circular | Moran's I | All neighbors in cumulative disk [0, r_outer) |
| `delta_i_ring_ind` | Ring | I_ND | Ring weights with cosine similarity |
| `delta_i_cir_ind` | Circular | I_ND | Circular weights with cosine similarity |

### Unified Interface (Recommended)

Compute delta I matrix for any combination of weight type and correlation type:

```r
library(sigdiscov)

# Load and prepare data
visium <- load_visium_data("path/to/spaceranger_output")
visium <- filter_in_tissue(visium)
expr_vst <- vst_transform(visium$counts)
coords <- get_spot_coords(visium)

# Compute ring-based Moran's I delta matrix (default)
result_ring_moran <- compute_delta_I_matrix_unified(
  expr_matrix = expr_vst,
  spot_coords = coords,
  weight_type = "ring",           # "ring" or "circular"
  correlation_type = "moran",     # "moran" or "ind"
  radii = seq(150, 650, 100),     # Distance bins (Visium ~100 unit spacing)
  coord_scale = 1,                # Scale factor for coordinates
  chunk_size = 1000               # Memory-efficient chunking
)

# Access the matrix (rows = targets, cols = factors)
delta_mat <- result_ring_moran$delta_I_signed

# Compute circular I_ND delta matrix
result_cir_ind <- compute_delta_I_matrix_unified(
  expr_matrix = expr_vst,
  spot_coords = coords,
  weight_type = "circular",
  correlation_type = "ind"
)
```

### Compute All 4 Types for Specific Factors

For focused analysis on specific factor genes:

```r
# Compute all 4 delta types for IFNG and TGFB1
results <- compute_four_deltas(
  expr_matrix = expr_vst,
  spot_coords = coords,
  factor_genes = c("IFNG", "TGFB1"),
  radii = seq(150, 650, 100)
)

# Access results for IFNG
ifng_deltas <- results$IFNG
head(ifng_deltas)
#   gene      delta_i_ring_moran  delta_i_cir_moran  delta_i_ring_ind  delta_i_cir_ind
# 1 GeneA     0.0234              0.0198             0.0312            0.0156
# ...

# Compare correlations between delta types
cor(ifng_deltas[, c("delta_i_ring_moran", "delta_i_cir_moran",
                     "delta_i_ring_ind", "delta_i_cir_ind")])
```

### Single Factor Analysis

For detailed analysis of one factor:

```r
# Compute signatures for all genes against IL1B
signatures <- compute_signed_delta_I(
  expr_matrix = expr_vst,
  spot_coords = coords,
  factor_gene = "IL1B",
  mode = "bivariate",  # or "directional" for I_ND
  radii = seq(100, 600, 100)
)

# View top responders (decay pattern = positive delta_I_signed)
head(signatures[order(-signatures$delta_I_signed), ])

# View avoiders (increase pattern = negative delta_I_signed)
head(signatures[order(signatures$delta_I_signed), ])
```

### Visualizing I(r) Curves

```r
# Get detailed curve for a gene pair
curve <- get_moran_curve(
  expr_matrix = expr_vst,
  spot_coords = coords,
  factor_gene = "IL1B",
  target_gene = "COL1A1",
  mode = "bivariate"
)

# Plot the curve
plot_moran_curve(curve)

# Manual plotting
plot(curve$radii, curve$I_raw, type = "p", pch = 19,
     xlab = "Distance", ylab = "Moran's I",
     main = paste(curve$factor_gene, "->", curve$target_gene))
lines(curve$radii, curve$I_smooth, col = "blue", lwd = 2)
abline(h = 0, lty = 2)
legend("topright",
       legend = paste("delta_I_signed =", round(curve$delta_I_signed, 4)))
```

## Functions

### Delta I Computation

| Function | Description |
|----------|-------------|
| `compute_delta_I_matrix_unified()` | **Unified interface** for computing delta I matrix with any weight/correlation type |
| `compute_four_deltas()` | Compute all 4 delta types for specific factor genes |
| `compute_signed_delta_I()` | Compute delta I for one factor against all genes |
| `compute_delta_I_matrix()` | Legacy interface for full delta I matrix |

### Curve Analysis

| Function | Description |
|----------|-------------|
| `get_moran_curve()` | Get full I(r) curve for a gene pair |
| `plot_moran_curve()` | Visualize I(r) curve with smoothing |

### Pairwise Moran's I

| Function | Description |
|----------|-------------|
| `pairwise_moran()` | Compute pairwise Moran's I between all genes |
| `moran_I()` | Compute Moran's I for a single gene pair |
| `create_weight_matrix()` | Create spatial weight matrix |

### Data Loading & Preprocessing

| Function | Description |
|----------|-------------|
| `load_visium_data()` | Load Visium data from Space Ranger output |
| `filter_in_tissue()` | Filter to in-tissue spots only |
| `get_spot_coords()` | Extract spot coordinates |
| `vst_transform()` | Apply VST normalization |
| `run_pipeline()` | Complete pipeline: load, preprocess, compute |
| `parse_spot_names()` | Parse "ROWxCOL" format spot names |

### I/O

| Function | Description |
|----------|-------------|
| `save_moran_result()` | Save results in lower triangular format |
| `load_moran_result()` | Load results from file |

## Method Details

### Weight Types

- **Ring weights**: Neighbors within distance band [r_inner, r_outer). Captures distance-SPECIFIC spatial structure at each radius.
- **Circular weights**: ALL neighbors within cumulative disk [0, r_outer). Provides distance-CUMULATIVE spatial structure.

### Correlation Types

- **Moran's I**: Bivariate spatial autocorrelation. `I = z_f' * W * z_g / n`
- **I_ND**: Directional cosine similarity. `I_ND = dot(z_f, W*z_g) / (||z_f|| * ||W*z_g||)`

### Delta I Computation

1. Compute I(r) at each distance radius
2. Smooth curve with Savitzky-Golay filter (removes noise, preserves trend)
3. Compute: `delta_I_signed = sign * (I_max - I_min)`
   - sign = +1 if I_short > I_long (decay pattern = RESPONDER)
   - sign = -1 if I_short < I_long (increase pattern = AVOIDANCE)

### Interpretation

| Pattern | delta_I_signed | Interpretation |
|---------|----------------|----------------|
| Decay | > 0 | High correlation nearby, low far away -> TRUE RESPONDER |
| Increase | < 0 | Low correlation nearby, high far away -> AVOIDANCE |
| Flat | ~ 0 | Constant correlation -> CONSTITUTIVE EXPRESSION |

## Performance

The package uses BLAS matrix operations for efficiency:

- **Algorithm**: Reformulated as `I = Z * W * Z^T / n` using matrix multiplication
- **Optimization**: Weight matrices precomputed once, reused for all gene pairs
- **Memory**: Chunked processing for large gene sets (>5000 genes)
- **Typical performance**: 19,729 genes x 3,813 spots completes in ~2 minutes (R) or ~15 seconds (C++ standalone)

## Citation

If you use this package, please cite:

Beibei Ru, Lanqi Gong, Emily Yang, Seongyong Park, George Zaki, Kenneth Aldape, Lalage Wakefield, Peng Jiang. Inference of secreted protein activities in intercellular communication. [Link](https://github.com/data2intelligence/SecAct)

## License

MIT License
