# =============================================================================
# Test: Signed Delta I with Synthetic Data
# =============================================================================
# This test creates synthetic spatial data with known patterns and verifies:
# 1. compute_signed_delta_I() works in both bivariate and directional modes
# 2. I_ND values are bounded [-1, 1]
# 3. Decaying spatial patterns give positive delta_I_signed
# 4. Increasing spatial patterns give negative delta_I_signed
# =============================================================================

library(sigdiscov)

cat("=============================================================================\n")
cat("Test: Signed Delta I with Synthetic Spatial Data\n")
cat("=============================================================================\n\n")

# -----------------------------------------------------------------------------
# Step 1: Create synthetic spatial grid (like Visium)
# -----------------------------------------------------------------------------
cat("Step 1: Creating synthetic Visium-like grid...\n")

# Create a 20x20 hexagonal grid (similar to Visium layout)
set.seed(42)
grid_size <- 20
spots <- expand.grid(row = 1:grid_size, col = 1:grid_size)

# Filter to hexagonal pattern (every other column offset)
spots <- spots[which((spots$row + spots$col) %% 2 == 0), ]
n_spots <- nrow(spots)
cat("  Grid:", grid_size, "x", grid_size, "->", n_spots, "spots\n")

# Convert to physical coordinates (100um spacing for Visium)
coords <- as.matrix(spots) * 100  # micrometers

# Define center of the tissue
center_x <- mean(coords[, 1])
center_y <- mean(coords[, 2])

# Calculate distance from center for each spot
dist_from_center <- sqrt((coords[, 1] - center_x)^2 + (coords[, 2] - center_y)^2)
dist_normalized <- dist_from_center / max(dist_from_center)  # 0-1 range

cat("  Center: (", round(center_x), ",", round(center_y), ") um\n")
cat("  Max distance from center:", round(max(dist_from_center)), "um\n\n")

# -----------------------------------------------------------------------------
# Step 2: Create synthetic gene expression with known spatial patterns
# -----------------------------------------------------------------------------
cat("Step 2: Creating synthetic gene expression patterns...\n")

n_genes <- 50
gene_names <- paste0("Gene_", 1:n_genes)

# Initialize expression matrix
expr_matrix <- matrix(0, nrow = n_genes, ncol = n_spots)
rownames(expr_matrix) <- gene_names
colnames(expr_matrix) <- paste0(spots$row, "x", spots$col)

# Gene 1: Factor gene - high expression at center, decaying outward
# This simulates a signaling source at the center
expr_matrix[1, ] <- 10 * exp(-3 * dist_normalized) + rnorm(n_spots, 0, 0.5)
expr_matrix[1, ] <- pmax(expr_matrix[1, ], 0)  # Ensure non-negative

# Genes 2-10: TRUE RESPONDERS - high near factor source, low far away
# These should have POSITIVE delta_I_signed (decay pattern)
for (i in 2:10) {
  # Expression correlates with factor at short distance, decays at long distance
  lag <- 1 + (i - 2) * 0.5  # Slight lag for variety
  response <- 5 * exp(-2 * pmax(dist_normalized - 0.1 * lag, 0)) + rnorm(n_spots, 0, 0.3)
  expr_matrix[i, ] <- pmax(response, 0)
}

# Genes 11-20: AVOIDANCE genes - low near factor source, high far away
# These should have NEGATIVE delta_I_signed (increase pattern)
for (i in 11:20) {
  # Expression inversely correlates - low near center, high at edges
  avoidance <- 5 * (1 - exp(-3 * dist_normalized)) + rnorm(n_spots, 0, 0.3)
  expr_matrix[i, ] <- pmax(avoidance, 0)
}

# Genes 21-35: CONSTITUTIVE genes - uniform expression, no spatial pattern
# These should have delta_I_signed near ZERO
for (i in 21:35) {
  expr_matrix[i, ] <- 3 + rnorm(n_spots, 0, 1)
  expr_matrix[i, ] <- pmax(expr_matrix[i, ], 0)
}

# Genes 36-50: Random spatial patterns (noise)
for (i in 36:n_genes) {
  expr_matrix[i, ] <- abs(rnorm(n_spots, 3, 2))
}

cat("  Factor gene (Gene_1): Central source pattern\n")
cat("  Genes 2-10: Responder pattern (decay from source)\n")
cat("  Genes 11-20: Avoidance pattern (increase from source)\n")
cat("  Genes 21-35: Constitutive pattern (uniform)\n")
cat("  Genes 36-50: Random noise\n\n")

# -----------------------------------------------------------------------------
# Step 3: Run compute_signed_delta_I in BIVARIATE mode
# -----------------------------------------------------------------------------
cat("Step 3: Testing BIVARIATE mode...\n")

radii <- seq(100, 600, 100)

sig_bivar <- compute_signed_delta_I(
  expr_matrix = expr_matrix,
  spot_coords = spots,
  factor_gene = 1,
  radii = radii,
  mode = "bivariate",
  coord_scale = 100,
  verbose = FALSE
)

cat("  Result: ", nrow(sig_bivar), " genes\n")

# Check responder genes (should have positive delta_I_signed)
responder_sig <- sig_bivar$delta_I_signed[2:10]
cat("  Responder genes (2-10) delta_I_signed:\n")
cat("    Mean:", round(mean(responder_sig), 4), "\n")
cat("    All positive:", all(responder_sig > 0), "\n")

# Check avoidance genes (should have negative delta_I_signed)
avoidance_sig <- sig_bivar$delta_I_signed[11:20]
cat("  Avoidance genes (11-20) delta_I_signed:\n")
cat("    Mean:", round(mean(avoidance_sig), 4), "\n")
cat("    All negative:", all(avoidance_sig < 0), "\n")

# Check constitutive genes (should have small |delta_I_signed|)
const_sig <- sig_bivar$delta_I_signed[21:35]
cat("  Constitutive genes (21-35) |delta_I_signed|:\n")
cat("    Mean abs:", round(mean(abs(const_sig)), 4), "\n")
cat("    Max abs:", round(max(abs(const_sig)), 4), "\n\n")

# -----------------------------------------------------------------------------
# Step 4: Run compute_signed_delta_I in DIRECTIONAL mode
# -----------------------------------------------------------------------------
cat("Step 4: Testing DIRECTIONAL mode...\n")

sig_dir <- compute_signed_delta_I(
  expr_matrix = expr_matrix,
  spot_coords = spots,
  factor_gene = 1,
  radii = radii,
  mode = "directional",
  sender_percentile = 75,
  coord_scale = 100,
  verbose = FALSE
)

cat("  Result: ", nrow(sig_dir), " genes\n")

# Check responder genes
responder_sig_dir <- sig_dir$delta_I_signed[2:10]
cat("  Responder genes (2-10) delta_I_signed:\n")
cat("    Mean:", round(mean(responder_sig_dir, na.rm = TRUE), 4), "\n")

# Check avoidance genes
avoidance_sig_dir <- sig_dir$delta_I_signed[11:20]
cat("  Avoidance genes (11-20) delta_I_signed:\n")
cat("    Mean:", round(mean(avoidance_sig_dir, na.rm = TRUE), 4), "\n\n")

# -----------------------------------------------------------------------------
# Step 5: Verify I_ND bounds [-1, 1] using get_moran_curve
# -----------------------------------------------------------------------------
cat("Step 5: Verifying I_ND bounds [-1, 1]...\n")

# Get raw I curves for several genes in directional mode
all_I_values <- c()
for (i in 1:20) {
  curve <- get_moran_curve(
    expr_matrix = expr_matrix,
    spot_coords = spots,
    factor_gene = 1,
    target_gene = i,
    radii = radii,
    mode = "directional",
    coord_scale = 100
  )
  all_I_values <- c(all_I_values, curve$I_raw[!is.na(curve$I_raw)])
}

I_min <- min(all_I_values, na.rm = TRUE)
I_max <- max(all_I_values, na.rm = TRUE)

cat("  I_ND range across all curves:\n")
cat("    Min:", round(I_min, 4), "\n")
cat("    Max:", round(I_max, 4), "\n")
cat("    Within [-1, 1]:", I_min >= -1.001 && I_max <= 1.001, "\n\n")

# -----------------------------------------------------------------------------
# Step 6: Visualize a few representative curves
# -----------------------------------------------------------------------------
cat("Step 6: Generating test curves...\n")

# Get curves for: responder, avoidance, and constitutive genes
curve_responder <- get_moran_curve(expr_matrix, spots, 1, 2, radii, "bivariate", coord_scale = 100)
curve_avoidance <- get_moran_curve(expr_matrix, spots, 1, 15, radii, "bivariate", coord_scale = 100)
curve_const <- get_moran_curve(expr_matrix, spots, 1, 25, radii, "bivariate", coord_scale = 100)

cat("  Responder (Gene_2):\n")
cat("    I_raw:", paste(round(curve_responder$I_raw, 4), collapse = ", "), "\n")
cat("    delta_I_signed:", round(curve_responder$delta_I_signed, 4), "(expected: positive)\n")

cat("  Avoidance (Gene_15):\n")
cat("    I_raw:", paste(round(curve_avoidance$I_raw, 4), collapse = ", "), "\n")
cat("    delta_I_signed:", round(curve_avoidance$delta_I_signed, 4), "(expected: negative)\n")

cat("  Constitutive (Gene_25):\n")
cat("    I_raw:", paste(round(curve_const$I_raw, 4), collapse = ", "), "\n")
cat("    delta_I_signed:", round(curve_const$delta_I_signed, 4), "(expected: ~0)\n\n")

# -----------------------------------------------------------------------------
# Step 7: Run assertions
# -----------------------------------------------------------------------------
cat("Step 7: Running assertions...\n")

tests_passed <- 0
tests_failed <- 0

# Test 1: Responder genes have positive delta_I_signed (bivariate)
if (mean(responder_sig) > 0) {
  cat("  [PASS] Responder genes have positive mean delta_I_signed (bivariate)\n")
  tests_passed <- tests_passed + 1
} else {
  cat("  [FAIL] Responder genes should have positive mean delta_I_signed\n")
  tests_failed <- tests_failed + 1
}

# Test 2: Avoidance genes have negative delta_I_signed (bivariate)
if (mean(avoidance_sig) < 0) {
  cat("  [PASS] Avoidance genes have negative mean delta_I_signed (bivariate)\n")
  tests_passed <- tests_passed + 1
} else {
  cat("  [FAIL] Avoidance genes should have negative mean delta_I_signed\n")
  tests_failed <- tests_failed + 1
}

# Test 3: Constitutive genes have small |delta_I_signed|
if (mean(abs(const_sig)) < mean(abs(responder_sig)) * 0.5) {
  cat("  [PASS] Constitutive genes have smaller |delta_I_signed| than responders\n")
  tests_passed <- tests_passed + 1
} else {
  cat("  [FAIL] Constitutive genes should have smaller |delta_I_signed|\n")
  tests_failed <- tests_failed + 1
}

# Test 4: I_ND bounded [-1, 1]
if (I_min >= -1.001 && I_max <= 1.001) {
  cat("  [PASS] I_ND values are bounded [-1, 1]\n")
  tests_passed <- tests_passed + 1
} else {
  cat("  [FAIL] I_ND values should be bounded [-1, 1]\n")
  tests_failed <- tests_failed + 1
}

# Test 5: Responder curve has positive delta_I_signed
if (curve_responder$delta_I_signed > 0) {
  cat("  [PASS] Single responder curve has positive delta_I_signed\n")
  tests_passed <- tests_passed + 1
} else {
  cat("  [FAIL] Responder curve should have positive delta_I_signed\n")
  tests_failed <- tests_failed + 1
}

# Test 6: Avoidance curve has negative delta_I_signed
if (curve_avoidance$delta_I_signed < 0) {
  cat("  [PASS] Single avoidance curve has negative delta_I_signed\n")
  tests_passed <- tests_passed + 1
} else {
  cat("  [FAIL] Avoidance curve should have negative delta_I_signed\n")
  tests_failed <- tests_failed + 1
}

# Test 7: DataFrame has expected columns
expected_cols <- c("gene", "delta_I", "delta_I_signed", "sign", "I_short", "I_long")
if (all(expected_cols %in% names(sig_bivar))) {
  cat("  [PASS] Result DataFrame has expected columns\n")
  tests_passed <- tests_passed + 1
} else {
  cat("  [FAIL] Result DataFrame missing expected columns\n")
  tests_failed <- tests_failed + 1
}

# Test 8: sign column matches delta_I_signed direction
signs_match <- all(sign(sig_bivar$delta_I_signed) == sig_bivar$sign | sig_bivar$delta_I_signed == 0)
if (signs_match) {
  cat("  [PASS] sign column matches delta_I_signed direction\n")
  tests_passed <- tests_passed + 1
} else {
  cat("  [FAIL] sign column should match delta_I_signed direction\n")
  tests_failed <- tests_failed + 1
}

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
cat("\n=============================================================================\n")
cat("Test Summary\n")
cat("=============================================================================\n")
cat("Passed:", tests_passed, "/", tests_passed + tests_failed, "\n")
cat("Failed:", tests_failed, "/", tests_passed + tests_failed, "\n")

if (tests_failed == 0) {
  cat("\nAll tests PASSED!\n")
} else {
  cat("\nSome tests FAILED. Please review.\n")
  quit(status = 1)
}
