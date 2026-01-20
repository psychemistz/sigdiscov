// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <vector>

using namespace Rcpp;
using namespace arma;

// =============================================================================
// Constants for spatial correlation modes
// =============================================================================
#define MODE_BIVARIATE 0    // Symmetric, all spots participate equally
#define MODE_DIRECTIONAL 1  // Sender -> Receiver directionality (I_ND)

// =============================================================================
// Helper: Compute Euclidean distance between two coordinates
// =============================================================================
inline double euclidean_dist(double x1, double y1, double x2, double y2) {
    double dx = x1 - x2;
    double dy = y1 - y2;
    return std::sqrt(dx * dx + dy * dy);
}

// =============================================================================
// Create ring weight matrix for a specific distance band
//
// Creates a weight matrix where entries are non-zero only for spot pairs
// within the specified distance ring [r_inner, r_outer).
//
// Supports two modes:
//   BIVARIATE (mode=0): N x N matrix, all spots participate symmetrically
//   DIRECTIONAL (mode=1): n_senders x n_receivers matrix for I_ND computation
//
// @param coords Nx2 matrix of spot coordinates (x, y or row, col)
// @param r_inner Inner radius of ring (inclusive)
// @param r_outer Outer radius of ring (exclusive)
// @param sender_idx Indices of sender spots (0-based, for directional mode)
// @param receiver_idx Indices of receiver spots (0-based, for directional mode)
// @param mode 0 = BIVARIATE, 1 = DIRECTIONAL
// @param row_normalize Whether to row-normalize the weights (default: true)
//
// @return Sparse weight matrix (NxN for bivariate, n_senders x n_receivers for directional)
// =============================================================================
// [[Rcpp::export]]
arma::sp_mat cpp_create_ring_weight_matrix(
    const arma::mat& coords,
    double r_inner,
    double r_outer,
    const arma::uvec& sender_idx,
    const arma::uvec& receiver_idx,
    int mode = 0,
    bool row_normalize = true
) {
    if (mode == MODE_BIVARIATE) {
        // BIVARIATE: N x N matrix, all spots
        int n = coords.n_rows;

        // Build sparse matrix using batch insert for efficiency
        std::vector<arma::uword> row_indices;
        std::vector<arma::uword> col_indices;
        std::vector<double> values;

        // First pass: find all neighbors and compute row sums
        std::vector<std::vector<int>> neighbors(n);
        std::vector<double> row_sums(n, 0.0);

        for (int i = 0; i < n; i++) {
            double xi = coords(i, 0);
            double yi = coords(i, 1);

            for (int j = 0; j < n; j++) {
                if (i == j) continue;  // Exclude self

                double dist = euclidean_dist(xi, yi, coords(j, 0), coords(j, 1));

                if (dist >= r_inner && dist < r_outer) {
                    neighbors[i].push_back(j);
                    row_sums[i] += 1.0;
                }
            }
        }

        // Second pass: build sparse matrix
        for (int i = 0; i < n; i++) {
            double norm_factor = (row_normalize && row_sums[i] > 0) ? row_sums[i] : 1.0;

            for (int j : neighbors[i]) {
                row_indices.push_back(i);
                col_indices.push_back(j);
                values.push_back(1.0 / norm_factor);
            }
        }

        // Create sparse matrix
        arma::umat locations(2, row_indices.size());
        for (size_t k = 0; k < row_indices.size(); k++) {
            locations(0, k) = row_indices[k];
            locations(1, k) = col_indices[k];
        }

        arma::vec vals(values);
        arma::sp_mat W(locations, vals, n, n);

        return W;

    } else {
        // DIRECTIONAL: n_senders x n_receivers matrix
        int n_senders = sender_idx.n_elem;
        int n_receivers = receiver_idx.n_elem;

        std::vector<arma::uword> row_indices;
        std::vector<arma::uword> col_indices;
        std::vector<double> values;

        // First pass: find neighbors and row sums
        std::vector<std::vector<int>> neighbors(n_senders);
        std::vector<double> row_sums(n_senders, 0.0);

        for (int i = 0; i < n_senders; i++) {
            int si = sender_idx(i);  // Actual spot index in coords
            double xi = coords(si, 0);
            double yi = coords(si, 1);

            for (int j = 0; j < n_receivers; j++) {
                int rj = receiver_idx(j);  // Actual spot index in coords

                double dist = euclidean_dist(xi, yi, coords(rj, 0), coords(rj, 1));

                if (dist >= r_inner && dist < r_outer) {
                    neighbors[i].push_back(j);  // j is index in receiver array
                    row_sums[i] += 1.0;
                }
            }
        }

        // Second pass: build sparse matrix
        for (int i = 0; i < n_senders; i++) {
            double norm_factor = (row_normalize && row_sums[i] > 0) ? row_sums[i] : 1.0;

            for (int j : neighbors[i]) {
                row_indices.push_back(i);
                col_indices.push_back(j);
                values.push_back(1.0 / norm_factor);
            }
        }

        // Create sparse matrix
        if (row_indices.empty()) {
            // No neighbors found - return empty sparse matrix
            arma::sp_mat W(n_senders, n_receivers);
            return W;
        }

        arma::umat locations(2, row_indices.size());
        for (size_t k = 0; k < row_indices.size(); k++) {
            locations(0, k) = row_indices[k];
            locations(1, k) = col_indices[k];
        }

        arma::vec vals(values);
        arma::sp_mat W(locations, vals, n_senders, n_receivers);

        return W;
    }
}

// =============================================================================
// Compute spatial correlation for a single distance ring
//
// Supports two measures:
//   BIVARIATE Moran's I: I = z_f' × W × z_g / n
//   DIRECTIONAL I_ND (cosine similarity): I_ND = (z_U' × W × z_V) / (||z_U|| × ||W × z_V||)
//
// @param z_f Standardized factor expression (length n for bivariate, n_senders for directional)
// @param z_g Standardized gene expression (length n for bivariate, n_receivers for directional)
// @param W Weight matrix (NxN for bivariate, n_senders x n_receivers for directional)
// @param mode 0 = BIVARIATE, 1 = DIRECTIONAL (I_ND)
//
// @return Scalar correlation value
// =============================================================================
// [[Rcpp::export]]
double cpp_compute_spatial_correlation(
    const arma::vec& z_f,
    const arma::vec& z_g,
    const arma::sp_mat& W,
    int mode = 0
) {
    // Compute spatial lag of gene: lag_g = W × z_g
    arma::vec lag_g = W * z_g;

    if (mode == MODE_BIVARIATE) {
        // BIVARIATE Moran's I: I = z_f' × W × z_g / n
        // This is equivalent to: I = dot(z_f, lag_g) / n
        double I = arma::dot(z_f, lag_g) / static_cast<double>(z_f.n_elem);
        return I;

    } else {
        // DIRECTIONAL I_ND: Cosine similarity
        //
        //              z_U^T × (W × z_V)
        // I_ND = ─────────────────────────────
        //         ||z_U|| × ||W × z_V||
        //
        double numerator = arma::dot(z_f, lag_g);
        double norm_zf = arma::norm(z_f, 2);
        double norm_lag = arma::norm(lag_g, 2);
        double denominator = norm_zf * norm_lag;

        // Avoid division by zero
        if (denominator < 1e-10) {
            return NA_REAL;
        }

        double I_ND = numerator / denominator;
        return I_ND;
    }
}

// =============================================================================
// Savitzky-Golay smoothing filter
//
// Applies polynomial smoothing to remove spikes and noise from the I(r) curve,
// revealing the underlying trend for reliable ΔI computation.
//
// The filter fits a polynomial of specified order to a window of points and
// uses the fitted value at the center point. This preserves overall trend
// direction while smoothing local fluctuations.
//
// @param y Input vector to smooth (raw I(r) curve)
// @param window_size Window size for the filter (must be odd, default: 5)
// @param poly_order Polynomial order for fitting (default: 2)
//
// @return Smoothed vector of same length (spike-free trend)
// =============================================================================
// [[Rcpp::export]]
arma::vec cpp_savgol_smooth(
    const arma::vec& y,
    int window_size = 5,
    int poly_order = 2
) {
    int n = y.n_elem;

    // Handle edge cases
    if (n <= 1) return y;

    // Ensure window is odd
    if (window_size % 2 == 0) {
        window_size++;
    }

    // Ensure window <= n and adjust if necessary
    if (window_size > n) {
        window_size = (n % 2 == 1) ? n : n - 1;  // Largest odd <= n
    }

    // Minimum window check
    if (window_size < 3) {
        return y;  // Can't smooth with less than 3 points
    }

    // Ensure poly_order is valid
    if (poly_order >= window_size) {
        poly_order = window_size - 1;
    }
    if (poly_order < 1) {
        poly_order = 1;
    }

    int half_window = window_size / 2;
    arma::vec y_smooth(n);

    // Precompute Savitzky-Golay coefficients
    // Build Vandermonde-like matrix J where J[i,j] = (i - half_window)^j
    arma::mat J(window_size, poly_order + 1);
    for (int i = 0; i < window_size; i++) {
        for (int j = 0; j <= poly_order; j++) {
            J(i, j) = std::pow(static_cast<double>(i - half_window), j);
        }
    }

    // Compute filter coefficients: first row of (J'J)^-1 J' gives smoothed center value
    arma::mat JtJ = J.t() * J;
    arma::mat JtJ_inv;

    // Check if matrix is invertible
    if (!arma::inv(JtJ_inv, JtJ)) {
        // Fallback: return original if matrix is singular
        return y;
    }

    arma::mat JtJ_inv_Jt = JtJ_inv * J.t();
    arma::vec coeffs = JtJ_inv_Jt.row(0).t();  // Coefficients for smoothed value

    // Apply filter to each point
    for (int i = 0; i < n; i++) {
        double sum = 0.0;

        for (int j = 0; j < window_size; j++) {
            int idx = i - half_window + j;

            // Handle boundaries by reflection (mirror padding)
            if (idx < 0) {
                idx = -idx;
            }
            if (idx >= n) {
                idx = 2 * n - idx - 2;
            }

            // Safety check
            if (idx < 0) idx = 0;
            if (idx >= n) idx = n - 1;

            sum += coeffs(j) * y(idx);
        }

        y_smooth(i) = sum;
    }

    return y_smooth;
}

// =============================================================================
// Compute signed delta I from smoothed I(r) curve
//
// Derives signature metrics from the distance-dependent correlation curve:
//   - delta_I: unsigned magnitude (max - min)
//   - delta_I_signed: signed based on trend direction
//   - sign: +1 (decay/responder) or -1 (increase/avoidance)
//
// Interpretation:
//   ΔI_signed > 0: Decay pattern (high correlation near, low far) → TRUE RESPONDER
//   ΔI_signed < 0: Increase pattern (low near, high far) → AVOIDANCE
//   ΔI_signed ≈ 0: Flat pattern → CONSTITUTIVE EXPRESSION
//
// @param I_smooth Smoothed I(r) curve
// @param radii Distance values for each ring
//
// @return List with all signature metrics
// =============================================================================
// [[Rcpp::export]]
Rcpp::List cpp_compute_signed_delta_I(
    const arma::vec& I_smooth,
    const arma::vec& radii
) {
    int n = I_smooth.n_elem;

    if (n == 0) {
        return Rcpp::List::create(
            Rcpp::Named("delta_I") = NA_REAL,
            Rcpp::Named("delta_I_signed") = NA_REAL,
            Rcpp::Named("sign") = NA_INTEGER,
            Rcpp::Named("I_max") = NA_REAL,
            Rcpp::Named("I_min") = NA_REAL,
            Rcpp::Named("I_short") = NA_REAL,
            Rcpp::Named("I_long") = NA_REAL,
            Rcpp::Named("peak_idx") = NA_INTEGER,
            Rcpp::Named("auc") = NA_REAL
        );
    }

    // Basic statistics
    double I_max = arma::max(I_smooth);
    double I_min = arma::min(I_smooth);
    double delta_I = I_max - I_min;

    // Short and long distance values
    double I_short = I_smooth(0);
    double I_long = I_smooth(n - 1);

    // Sign based on trend: positive if decay (short > long), negative if increase
    double trend = I_short - I_long;
    int sign = (trend >= 0) ? 1 : -1;

    // Signed delta I
    double delta_I_signed = sign * delta_I;

    // Find peak index (maximum absolute value)
    arma::vec abs_I = arma::abs(I_smooth);
    arma::uword peak_idx = arma::index_max(abs_I);

    // Area under curve (trapezoidal rule)
    double auc = 0.0;
    if (n > 1 && radii.n_elem == static_cast<arma::uword>(n)) {
        for (int i = 0; i < n - 1; i++) {
            double dr = radii(i + 1) - radii(i);
            auc += 0.5 * (I_smooth(i) + I_smooth(i + 1)) * dr;
        }
    }

    return Rcpp::List::create(
        Rcpp::Named("delta_I") = delta_I,
        Rcpp::Named("delta_I_signed") = delta_I_signed,
        Rcpp::Named("sign") = sign,
        Rcpp::Named("I_max") = I_max,
        Rcpp::Named("I_min") = I_min,
        Rcpp::Named("I_short") = I_short,
        Rcpp::Named("I_long") = I_long,
        Rcpp::Named("peak_idx") = static_cast<int>(peak_idx) + 1,  // R is 1-indexed
        Rcpp::Named("auc") = auc
    );
}

// =============================================================================
// Helper: Standardize a vector (z-score normalization)
// =============================================================================
arma::vec standardize_vec(const arma::vec& x) {
    double mean_x = arma::mean(x);
    double std_x = arma::stddev(x, 1);  // 1 = normalize by N (population)

    if (std_x < 1e-10) {
        // Return zeros if no variance
        return arma::zeros<arma::vec>(x.n_elem);
    }

    return (x - mean_x) / std_x;
}

// =============================================================================
// Compute Moran's I curve across multiple distance rings
//
// This is a convenience function that computes spatial correlation at each
// distance ring and returns the full I(r) curve.
//
// @param factor_expr Factor expression vector (all spots)
// @param gene_expr Gene expression vector (all spots)
// @param coords Spot coordinates (Nx2 matrix)
// @param radii Distance bin edges (defines rings: [0, r1), [r1, r2), ...)
// @param mode 0 = BIVARIATE, 1 = DIRECTIONAL
// @param sender_idx Sender spot indices (for directional mode)
// @param receiver_idx Receiver spot indices (for directional mode)
//
// @return Vector of correlation values (one per ring)
// =============================================================================
// [[Rcpp::export]]
arma::vec cpp_compute_moran_curve(
    const arma::vec& factor_expr,
    const arma::vec& gene_expr,
    const arma::mat& coords,
    const arma::vec& radii,
    int mode = 0,
    const arma::uvec& sender_idx = arma::uvec(),
    const arma::uvec& receiver_idx = arma::uvec()
) {
    int n_rings = radii.n_elem;
    arma::vec I_curve(n_rings);

    // Prepare expression vectors based on mode
    arma::vec z_f, z_g;
    arma::uvec s_idx, r_idx;

    if (mode == MODE_DIRECTIONAL && sender_idx.n_elem > 0 && receiver_idx.n_elem > 0) {
        // Directional mode: use sender/receiver subsets
        s_idx = sender_idx;
        r_idx = receiver_idx;

        // Extract and standardize expressions for subsets
        arma::vec factor_sub = factor_expr.elem(sender_idx);
        arma::vec gene_sub = gene_expr.elem(receiver_idx);
        z_f = standardize_vec(factor_sub);
        z_g = standardize_vec(gene_sub);
    } else {
        // Bivariate mode: use all spots
        int n = coords.n_rows;
        s_idx = arma::regspace<arma::uvec>(0, n - 1);
        r_idx = s_idx;
        z_f = standardize_vec(factor_expr);
        z_g = standardize_vec(gene_expr);
    }

    // Compute correlation for each ring
    for (int r = 0; r < n_rings; r++) {
        double r_inner = (r == 0) ? 0.0 : radii(r - 1);
        double r_outer = radii(r);

        // Create weight matrix for this ring
        arma::sp_mat W = cpp_create_ring_weight_matrix(
            coords, r_inner, r_outer, s_idx, r_idx, mode, true
        );

        // Compute correlation
        I_curve(r) = cpp_compute_spatial_correlation(z_f, z_g, W, mode);
    }

    return I_curve;
}

// =============================================================================
// Define sender and receiver spots based on factor expression threshold
//
// For Visium data where we don't have single-cell resolution, we split spots
// into senders (high factor expression) and receivers (low factor expression).
//
// @param factor_expr Factor expression vector for all spots
// @param percentile Percentile threshold for sender definition (default: 75)
//                   Spots >= this percentile are senders
//
// @return List with sender_idx, receiver_idx, threshold, n_senders, n_receivers
// =============================================================================
// [[Rcpp::export]]
Rcpp::List cpp_define_sender_receiver(
    const arma::vec& factor_expr,
    double percentile = 75.0
) {
    // Find expressed spots (> 0) for percentile calculation
    arma::uvec expressed_idx = arma::find(factor_expr > 0);

    double threshold;
    if (expressed_idx.n_elem > 0) {
        // Get values of expressed spots
        arma::vec expressed_vals = factor_expr.elem(expressed_idx);
        arma::vec sorted_vals = arma::sort(expressed_vals);

        // Calculate percentile index
        int pct_idx = static_cast<int>(sorted_vals.n_elem * percentile / 100.0);
        if (pct_idx >= static_cast<int>(sorted_vals.n_elem)) {
            pct_idx = sorted_vals.n_elem - 1;
        }
        threshold = sorted_vals(pct_idx);
    } else {
        // No expressed spots - use 0 as threshold
        threshold = 0.0;
    }

    // Split spots based on threshold
    arma::uvec sender_idx = arma::find(factor_expr >= threshold);
    arma::uvec receiver_idx = arma::find(factor_expr < threshold);

    return Rcpp::List::create(
        Rcpp::Named("sender_idx") = sender_idx + 1,  // Convert to 1-indexed for R
        Rcpp::Named("receiver_idx") = receiver_idx + 1,
        Rcpp::Named("threshold") = threshold,
        Rcpp::Named("n_senders") = static_cast<int>(sender_idx.n_elem),
        Rcpp::Named("n_receivers") = static_cast<int>(receiver_idx.n_elem)
    );
}

// =============================================================================
// Compute full signature for a single factor-gene pair
//
// Complete pipeline: compute I(r) curve, smooth it, derive signed ΔI.
//
// @param factor_expr Factor expression vector (all spots)
// @param gene_expr Gene expression vector (all spots)
// @param coords Spot coordinates (Nx2 matrix)
// @param radii Distance bin edges
// @param mode 0 = BIVARIATE, 1 = DIRECTIONAL
// @param sender_percentile Percentile for sender definition (directional mode)
// @param smooth_window Savitzky-Golay window size
// @param smooth_poly Savitzky-Golay polynomial order
//
// @return List with I_raw, I_smooth, and all signature metrics
// =============================================================================
// [[Rcpp::export]]
Rcpp::List cpp_compute_single_signature(
    const arma::vec& factor_expr,
    const arma::vec& gene_expr,
    const arma::mat& coords,
    const arma::vec& radii,
    int mode = 0,
    double sender_percentile = 75.0,
    int smooth_window = 5,
    int smooth_poly = 2
) {
    // Define sender/receiver for directional mode
    arma::uvec sender_idx, receiver_idx;
    int n_senders = 0, n_receivers = 0;
    double threshold = NA_REAL;

    if (mode == MODE_DIRECTIONAL) {
        Rcpp::List split = cpp_define_sender_receiver(factor_expr, sender_percentile);

        // Convert back to 0-indexed
        arma::uvec s_idx_r = split["sender_idx"];
        arma::uvec r_idx_r = split["receiver_idx"];
        sender_idx = s_idx_r - 1;
        receiver_idx = r_idx_r - 1;

        n_senders = split["n_senders"];
        n_receivers = split["n_receivers"];
        threshold = split["threshold"];
    }

    // Compute I(r) curve
    arma::vec I_raw = cpp_compute_moran_curve(
        factor_expr, gene_expr, coords, radii,
        mode, sender_idx, receiver_idx
    );

    // Smooth the curve
    arma::vec I_smooth = cpp_savgol_smooth(I_raw, smooth_window, smooth_poly);

    // Compute signed delta I
    Rcpp::List sig = cpp_compute_signed_delta_I(I_smooth, radii);

    // Return complete results
    return Rcpp::List::create(
        Rcpp::Named("I_raw") = I_raw,
        Rcpp::Named("I_smooth") = I_smooth,
        Rcpp::Named("radii") = radii,
        Rcpp::Named("delta_I") = sig["delta_I"],
        Rcpp::Named("delta_I_signed") = sig["delta_I_signed"],
        Rcpp::Named("sign") = sig["sign"],
        Rcpp::Named("I_max") = sig["I_max"],
        Rcpp::Named("I_min") = sig["I_min"],
        Rcpp::Named("I_short") = sig["I_short"],
        Rcpp::Named("I_long") = sig["I_long"],
        Rcpp::Named("peak_idx") = sig["peak_idx"],
        Rcpp::Named("auc") = sig["auc"],
        Rcpp::Named("mode") = mode,
        Rcpp::Named("n_senders") = n_senders,
        Rcpp::Named("n_receivers") = n_receivers,
        Rcpp::Named("sender_threshold") = threshold
    );
}

// =============================================================================
// Precompute ring weight matrices for all distance bins
//
// This function precomputes weight matrices for all distance rings, which can
// be reused when computing signatures for many genes against the same factor.
// This provides significant speedup for batch processing.
//
// @param coords Spot coordinates (Nx2 matrix)
// @param radii Distance bin edges
// @param mode 0 = BIVARIATE, 1 = DIRECTIONAL
// @param sender_idx Sender spot indices (0-based, for directional mode)
// @param receiver_idx Receiver spot indices (0-based, for directional mode)
//
// @return List of sparse weight matrices (one per ring)
// =============================================================================
// [[Rcpp::export]]
Rcpp::List cpp_precompute_ring_weights(
    const arma::mat& coords,
    const arma::vec& radii,
    int mode = 0,
    const arma::uvec& sender_idx = arma::uvec(),
    const arma::uvec& receiver_idx = arma::uvec()
) {
    int n_rings = radii.n_elem;
    Rcpp::List W_list(n_rings);

    // Determine indices for bivariate mode
    arma::uvec s_idx, r_idx;
    if (mode == MODE_BIVARIATE || sender_idx.n_elem == 0) {
        int n = coords.n_rows;
        s_idx = arma::regspace<arma::uvec>(0, n - 1);
        r_idx = s_idx;
    } else {
        s_idx = sender_idx;
        r_idx = receiver_idx;
    }

    // Precompute weight matrix for each ring
    for (int r = 0; r < n_rings; r++) {
        double r_inner = (r == 0) ? 0.0 : radii(r - 1);
        double r_outer = radii(r);

        arma::sp_mat W = cpp_create_ring_weight_matrix(
            coords, r_inner, r_outer, s_idx, r_idx, mode, true
        );

        W_list[r] = W;
    }

    return W_list;
}

// =============================================================================
// Internal helper: Compute I(r) curve using precomputed weight matrices
// =============================================================================
arma::vec compute_moran_curve_with_weights(
    const arma::vec& z_f,
    const arma::vec& z_g,
    const Rcpp::List& W_list,
    int mode
) {
    int n_rings = W_list.size();
    arma::vec I_curve(n_rings);

    for (int r = 0; r < n_rings; r++) {
        arma::sp_mat W = Rcpp::as<arma::sp_mat>(W_list[r]);
        I_curve(r) = cpp_compute_spatial_correlation(z_f, z_g, W, mode);
    }

    return I_curve;
}

// =============================================================================
// Internal helper: Compute signature from I curve (avoiding Rcpp::List overhead)
// =============================================================================
struct SignatureResult {
    double delta_I;
    double delta_I_signed;
    int sign;
    double I_max;
    double I_min;
    double I_short;
    double I_long;
    int peak_idx;
    double auc;
};

SignatureResult compute_signature_internal(
    const arma::vec& I_smooth,
    const arma::vec& radii
) {
    SignatureResult result;
    int n = I_smooth.n_elem;

    if (n == 0) {
        result.delta_I = NA_REAL;
        result.delta_I_signed = NA_REAL;
        result.sign = NA_INTEGER;
        result.I_max = NA_REAL;
        result.I_min = NA_REAL;
        result.I_short = NA_REAL;
        result.I_long = NA_REAL;
        result.peak_idx = NA_INTEGER;
        result.auc = NA_REAL;
        return result;
    }

    // Check for NAs in the smoothed curve
    bool has_na = false;
    for (int i = 0; i < n; i++) {
        if (!std::isfinite(I_smooth(i))) {
            has_na = true;
            break;
        }
    }

    if (has_na) {
        // Find valid (non-NA) values
        arma::uvec valid_idx = arma::find_finite(I_smooth);
        if (valid_idx.n_elem < 2) {
            result.delta_I = NA_REAL;
            result.delta_I_signed = NA_REAL;
            result.sign = NA_INTEGER;
            result.I_max = NA_REAL;
            result.I_min = NA_REAL;
            result.I_short = NA_REAL;
            result.I_long = NA_REAL;
            result.peak_idx = NA_INTEGER;
            result.auc = NA_REAL;
            return result;
        }

        arma::vec valid_I = I_smooth.elem(valid_idx);
        result.I_max = arma::max(valid_I);
        result.I_min = arma::min(valid_I);
        result.delta_I = result.I_max - result.I_min;
        result.I_short = valid_I(0);
        result.I_long = valid_I(valid_I.n_elem - 1);

        double trend = result.I_short - result.I_long;
        result.sign = (trend >= 0) ? 1 : -1;
        result.delta_I_signed = result.sign * result.delta_I;

        arma::vec abs_valid = arma::abs(valid_I);
        result.peak_idx = valid_idx(arma::index_max(abs_valid)) + 1;

        // AUC with valid values only
        result.auc = 0.0;
        for (arma::uword i = 0; i < valid_idx.n_elem - 1; i++) {
            int idx1 = valid_idx(i);
            int idx2 = valid_idx(i + 1);
            double dr = radii(idx2) - radii(idx1);
            result.auc += 0.5 * (valid_I(i) + valid_I(i + 1)) * dr;
        }
    } else {
        result.I_max = arma::max(I_smooth);
        result.I_min = arma::min(I_smooth);
        result.delta_I = result.I_max - result.I_min;
        result.I_short = I_smooth(0);
        result.I_long = I_smooth(n - 1);

        double trend = result.I_short - result.I_long;
        result.sign = (trend >= 0) ? 1 : -1;
        result.delta_I_signed = result.sign * result.delta_I;

        arma::vec abs_I = arma::abs(I_smooth);
        result.peak_idx = arma::index_max(abs_I) + 1;

        result.auc = 0.0;
        for (int i = 0; i < n - 1; i++) {
            double dr = radii(i + 1) - radii(i);
            result.auc += 0.5 * (I_smooth(i) + I_smooth(i + 1)) * dr;
        }
    }

    return result;
}

// =============================================================================
// Compute signatures for all genes against one factor (batch processing)
//
// This is the main batch processing function. It:
// 1. Precomputes weight matrices once for all distance rings
// 2. Standardizes the factor expression once
// 3. Loops through all genes efficiently
// 4. Returns results as a DataFrame
//
// @param expr_matrix Gene expression matrix (genes x spots), normalized
// @param gene_names Character vector of gene names
// @param factor_idx 0-based index of the factor gene in expr_matrix
// @param coords Spot coordinates (Nx2 matrix)
// @param radii Distance bin edges
// @param mode 0 = BIVARIATE, 1 = DIRECTIONAL
// @param sender_percentile Percentile for sender definition (directional mode)
// @param smooth_window Savitzky-Golay window size
// @param smooth_poly Savitzky-Golay polynomial order
// @param verbose Print progress messages
//
// @return DataFrame with signature metrics for each gene
// =============================================================================
// [[Rcpp::export]]
Rcpp::DataFrame cpp_compute_all_signatures(
    const arma::mat& expr_matrix,
    const Rcpp::CharacterVector& gene_names,
    int factor_idx,
    const arma::mat& coords,
    const arma::vec& radii,
    int mode = 0,
    double sender_percentile = 75.0,
    int smooth_window = 5,
    int smooth_poly = 2,
    bool verbose = true
) {
    int n_genes = expr_matrix.n_rows;
    int n_spots = expr_matrix.n_cols;
    int n_radii = radii.n_elem;

    if (verbose) {
        Rcpp::Rcout << "Computing signed delta I signatures" << std::endl;
        Rcpp::Rcout << n_spots << " spots, " << n_genes << " genes, "
                    << n_radii << " distance bins" << std::endl;
        Rcpp::Rcout << "Mode: " << (mode == 0 ? "bivariate" : "directional") << std::endl;
    }

    // Get factor expression
    arma::vec factor_expr = expr_matrix.row(factor_idx).t();

    // Define sender/receiver split for directional mode
    arma::uvec sender_idx, receiver_idx;
    int n_senders = n_spots, n_receivers = 0;
    double threshold = NA_REAL;

    if (mode == MODE_DIRECTIONAL) {
        Rcpp::List split = cpp_define_sender_receiver(factor_expr, sender_percentile);
        arma::uvec s_idx_r = split["sender_idx"];
        arma::uvec r_idx_r = split["receiver_idx"];
        sender_idx = s_idx_r - 1;  // Convert to 0-based
        receiver_idx = r_idx_r - 1;
        n_senders = split["n_senders"];
        n_receivers = split["n_receivers"];
        threshold = split["threshold"];

        if (verbose) {
            Rcpp::Rcout << "Sender/Receiver split: " << n_senders << " senders, "
                        << n_receivers << " receivers" << std::endl;
            Rcpp::Rcout << "Threshold: " << threshold << std::endl;
        }
    }

    // Precompute weight matrices
    if (verbose) Rcpp::Rcout << "Precomputing weight matrices..." << std::endl;
    Rcpp::List W_list = cpp_precompute_ring_weights(
        coords, radii, mode, sender_idx, receiver_idx
    );

    // Standardize factor expression (once)
    arma::vec z_f;
    if (mode == MODE_DIRECTIONAL) {
        arma::vec factor_sub = factor_expr.elem(sender_idx);
        z_f = standardize_vec(factor_sub);
    } else {
        z_f = standardize_vec(factor_expr);
    }

    // Result vectors
    Rcpp::NumericVector out_delta_I(n_genes);
    Rcpp::NumericVector out_delta_I_signed(n_genes);
    Rcpp::IntegerVector out_sign(n_genes);
    Rcpp::NumericVector out_I_short(n_genes);
    Rcpp::NumericVector out_I_long(n_genes);
    Rcpp::NumericVector out_I_max(n_genes);
    Rcpp::NumericVector out_I_min(n_genes);
    Rcpp::NumericVector out_auc(n_genes);
    Rcpp::IntegerVector out_peak_idx(n_genes);

    if (verbose) Rcpp::Rcout << "Computing signatures for " << n_genes << " genes..." << std::endl;

    // Process each gene
    for (int g = 0; g < n_genes; g++) {
        // Progress update
        if (verbose && ((g + 1) % 1000 == 0 || g == n_genes - 1)) {
            Rcpp::Rcout << "  Processed " << (g + 1) << " / " << n_genes << " genes" << std::endl;
        }

        // Check for user interrupt
        if (g % 100 == 0) {
            Rcpp::checkUserInterrupt();
        }

        // Get gene expression and standardize
        arma::vec gene_expr = expr_matrix.row(g).t();
        arma::vec z_g;
        if (mode == MODE_DIRECTIONAL) {
            arma::vec gene_sub = gene_expr.elem(receiver_idx);
            z_g = standardize_vec(gene_sub);
        } else {
            z_g = standardize_vec(gene_expr);
        }

        // Compute I(r) curve using precomputed weights
        arma::vec I_raw = compute_moran_curve_with_weights(z_f, z_g, W_list, mode);

        // Smooth the curve
        arma::vec I_smooth = cpp_savgol_smooth(I_raw, smooth_window, smooth_poly);

        // Compute signature
        SignatureResult sig = compute_signature_internal(I_smooth, radii);

        // Store results
        out_delta_I[g] = sig.delta_I;
        out_delta_I_signed[g] = sig.delta_I_signed;
        out_sign[g] = sig.sign;
        out_I_short[g] = sig.I_short;
        out_I_long[g] = sig.I_long;
        out_I_max[g] = sig.I_max;
        out_I_min[g] = sig.I_min;
        out_auc[g] = sig.auc;
        out_peak_idx[g] = sig.peak_idx;
    }

    if (verbose) Rcpp::Rcout << "Done!" << std::endl;

    // Create DataFrame
    return Rcpp::DataFrame::create(
        Rcpp::Named("gene") = gene_names,
        Rcpp::Named("delta_I") = out_delta_I,
        Rcpp::Named("delta_I_signed") = out_delta_I_signed,
        Rcpp::Named("sign") = out_sign,
        Rcpp::Named("I_short") = out_I_short,
        Rcpp::Named("I_long") = out_I_long,
        Rcpp::Named("I_max") = out_I_max,
        Rcpp::Named("I_min") = out_I_min,
        Rcpp::Named("auc") = out_auc,
        Rcpp::Named("peak_radius_idx") = out_peak_idx,
        Rcpp::Named("stringsAsFactors") = false
    );
}

// =============================================================================
// Compute full delta I matrix for ALL factor-gene pairs (highly optimized)
//
// This function is optimized for computing delta I for ALL gene pairs at once.
// Key optimization: In bivariate mode, the weight matrix W(r) is the SAME for
// all gene pairs at a given radius. So we:
//   1. Precompute weight matrices once for all radii (6 matrices)
//   2. Standardize all gene expressions once (Z matrix)
//   3. For each radius, compute full I matrix: I(r) = Z × W(r) × Z^T / n
//   4. Then compute delta I for all pairs from the I(r) cube
//
// This reduces weight matrix creation from O(n_factors × n_radii) to O(n_radii).
//
// @param expr_matrix Gene expression matrix (genes x spots), normalized
// @param coords Spot coordinates (Nx2 matrix)
// @param radii Distance bin edges
// @param smooth_window Savitzky-Golay window size
// @param smooth_poly Savitzky-Golay polynomial order
// @param verbose Print progress messages
//
// @return List with delta_I_signed matrix and additional info
// =============================================================================
// [[Rcpp::export]]
Rcpp::List cpp_compute_delta_I_matrix(
    const arma::mat& expr_matrix,
    const arma::mat& coords,
    const arma::vec& radii,
    int smooth_window = 5,
    int smooth_poly = 2,
    bool verbose = true
) {
    int n_genes = expr_matrix.n_rows;
    int n_spots = expr_matrix.n_cols;
    int n_radii = radii.n_elem;

    if (verbose) {
        Rcpp::Rcout << "Computing full delta I matrix (optimized)" << std::endl;
        Rcpp::Rcout << n_spots << " spots, " << n_genes << " genes, "
                    << n_radii << " distance bins" << std::endl;
    }

    // Step 1: Standardize all genes at once (Z matrix: n_genes x n_spots)
    if (verbose) Rcpp::Rcout << "Step 1: Standardizing expression matrix..." << std::endl;
    arma::mat Z(n_genes, n_spots);
    for (int g = 0; g < n_genes; g++) {
        arma::vec gene_expr = expr_matrix.row(g).t();
        Z.row(g) = standardize_vec(gene_expr).t();
    }

    // Step 2: Precompute weight matrices for all radii (only 6 matrices!)
    if (verbose) Rcpp::Rcout << "Step 2: Precomputing " << n_radii << " weight matrices..." << std::endl;

    // For bivariate mode, use all spots
    arma::uvec all_idx = arma::regspace<arma::uvec>(0, n_spots - 1);

    std::vector<arma::sp_mat> W_vec(n_radii);
    for (int r = 0; r < n_radii; r++) {
        double r_inner = (r == 0) ? 0.0 : radii(r - 1);
        double r_outer = radii(r);

        W_vec[r] = cpp_create_ring_weight_matrix(
            coords, r_inner, r_outer, all_idx, all_idx, 0, true
        );

        if (verbose) {
            Rcpp::Rcout << "  Ring " << (r+1) << ": [" << r_inner << ", " << r_outer << ") - "
                        << W_vec[r].n_nonzero << " non-zero weights" << std::endl;
        }
    }

    // Step 3: Compute I(r) for all pairs at each radius
    // I_cube is n_genes x n_genes x n_radii
    if (verbose) Rcpp::Rcout << "Step 3: Computing I matrices for each radius..." << std::endl;

    arma::cube I_cube(n_genes, n_genes, n_radii);

    for (int r = 0; r < n_radii; r++) {
        if (verbose) Rcpp::Rcout << "  Radius " << (r+1) << "/" << n_radii << "..." << std::endl;

        // Compute Z × W(r) first
        // Z is n_genes x n_spots, W is n_spots x n_spots
        // Z × W gives n_genes x n_spots (spatial lag of all genes)
        arma::mat ZW = Z * W_vec[r];

        // Now compute I_matrix(r) = Z × W × Z^T / n
        // = ZW × Z^T / n
        // Note: I[i,j] = z_i' * W * z_j has gene i as factor, gene j as target
        // We want rows = targets, cols = factors, so we transpose
        arma::mat I_r = (ZW * Z.t()).t() / static_cast<double>(n_spots);

        I_cube.slice(r) = I_r;

        // Check for user interrupt
        Rcpp::checkUserInterrupt();
    }

    // Step 4: Smooth I(r) curves and compute delta I for each pair
    if (verbose) Rcpp::Rcout << "Step 4: Smoothing curves and computing delta I..." << std::endl;

    arma::mat delta_I_matrix(n_genes, n_genes);

    int total_pairs = n_genes * n_genes;
    int processed = 0;

    for (int i = 0; i < n_genes; i++) {
        for (int j = 0; j < n_genes; j++) {
            // Extract I(r) curve for this pair
            arma::vec I_raw(n_radii);
            for (int r = 0; r < n_radii; r++) {
                I_raw(r) = I_cube(i, j, r);
            }

            // Smooth the curve
            arma::vec I_smooth = cpp_savgol_smooth(I_raw, smooth_window, smooth_poly);

            // Compute signed delta I
            SignatureResult sig = compute_signature_internal(I_smooth, radii);

            delta_I_matrix(i, j) = sig.delta_I_signed;

            processed++;
        }

        // Progress update and interrupt check
        if (verbose && ((i + 1) % 100 == 0 || i == n_genes - 1)) {
            Rcpp::Rcout << "  Processed " << (i + 1) << " / " << n_genes << " rows" << std::endl;
        }
        if (i % 10 == 0) {
            Rcpp::checkUserInterrupt();
        }
    }

    if (verbose) Rcpp::Rcout << "Done!" << std::endl;

    return Rcpp::List::create(
        Rcpp::Named("delta_I_signed") = delta_I_matrix,
        Rcpp::Named("n_genes") = n_genes,
        Rcpp::Named("n_spots") = n_spots,
        Rcpp::Named("n_radii") = n_radii,
        Rcpp::Named("radii") = radii
    );
}

// =============================================================================
// Memory-efficient version: Compute delta I matrix in chunks
//
// For large gene sets (>5000), the full I_cube won't fit in memory.
// This version processes genes in chunks, computing I matrices on-the-fly.
//
// @param expr_matrix Gene expression matrix (genes x spots), normalized
// @param coords Spot coordinates (Nx2 matrix)
// @param radii Distance bin edges
// @param chunk_size Number of target genes to process at once (default: 1000)
// @param smooth_window Savitzky-Golay window size
// @param smooth_poly Savitzky-Golay polynomial order
// @param verbose Print progress messages
//
// @return List with delta_I_signed matrix and additional info
// =============================================================================
// [[Rcpp::export]]
Rcpp::List cpp_compute_delta_I_matrix_chunked(
    const arma::mat& expr_matrix,
    const arma::mat& coords,
    const arma::vec& radii,
    int chunk_size = 1000,
    int smooth_window = 5,
    int smooth_poly = 2,
    bool verbose = true
) {
    int n_genes = expr_matrix.n_rows;
    int n_spots = expr_matrix.n_cols;
    int n_radii = radii.n_elem;

    if (verbose) {
        Rcpp::Rcout << "Computing full delta I matrix (memory-efficient chunked version)" << std::endl;
        Rcpp::Rcout << n_spots << " spots, " << n_genes << " genes, "
                    << n_radii << " distance bins" << std::endl;
        Rcpp::Rcout << "Chunk size: " << chunk_size << " genes" << std::endl;
    }

    // Step 1: Standardize all genes at once (Z matrix: n_genes x n_spots)
    if (verbose) Rcpp::Rcout << "Step 1: Standardizing expression matrix..." << std::endl;
    arma::mat Z(n_genes, n_spots);
    for (int g = 0; g < n_genes; g++) {
        arma::vec gene_expr = expr_matrix.row(g).t();
        Z.row(g) = standardize_vec(gene_expr).t();
    }

    // Step 2: Precompute weight matrices and ZW for all radii
    if (verbose) Rcpp::Rcout << "Step 2: Precomputing weight matrices and spatial lags..." << std::endl;

    arma::uvec all_idx = arma::regspace<arma::uvec>(0, n_spots - 1);

    // Store ZW matrices (n_genes x n_spots each) for all radii
    std::vector<arma::mat> ZW_vec(n_radii);
    for (int r = 0; r < n_radii; r++) {
        double r_inner = (r == 0) ? 0.0 : radii(r - 1);
        double r_outer = radii(r);

        arma::sp_mat W = cpp_create_ring_weight_matrix(
            coords, r_inner, r_outer, all_idx, all_idx, 0, true
        );

        if (verbose) {
            Rcpp::Rcout << "  Ring " << (r+1) << ": [" << r_inner << ", " << r_outer << ") - "
                        << W.n_nonzero << " non-zero weights" << std::endl;
        }

        // Compute ZW = Z * W (n_genes x n_spots)
        ZW_vec[r] = Z * W;
    }

    // Step 3: Allocate result matrix
    if (verbose) Rcpp::Rcout << "Step 3: Allocating result matrix..." << std::endl;
    arma::mat delta_I_matrix(n_genes, n_genes);

    // Step 4: Process in chunks of target genes
    int n_chunks = (n_genes + chunk_size - 1) / chunk_size;
    if (verbose) Rcpp::Rcout << "Step 4: Processing " << n_chunks << " chunks..." << std::endl;

    for (int chunk = 0; chunk < n_chunks; chunk++) {
        int start_idx = chunk * chunk_size;
        int end_idx = std::min(start_idx + chunk_size, n_genes);
        int chunk_genes = end_idx - start_idx;

        if (verbose) {
            Rcpp::Rcout << "  Chunk " << (chunk + 1) << "/" << n_chunks
                        << " (genes " << start_idx << "-" << (end_idx - 1) << ")..." << std::endl;
        }

        // Extract Z_chunk for this chunk of target genes
        arma::mat Z_chunk = Z.rows(start_idx, end_idx - 1);  // chunk_genes x n_spots

        // Compute I matrices for this chunk at all radii
        // I_chunk[target, factor, radius] where target is within chunk
        arma::cube I_chunk(chunk_genes, n_genes, n_radii);

        for (int r = 0; r < n_radii; r++) {
            // I_r = Z_chunk * ZW_vec[r]' / n
            // Z_chunk is chunk_genes x n_spots
            // ZW_vec[r] is n_genes x n_spots
            // Result is chunk_genes x n_genes (target x factor)
            arma::mat I_r = Z_chunk * ZW_vec[r].t() / static_cast<double>(n_spots);
            I_chunk.slice(r) = I_r;
        }

        // Smooth and compute delta I for this chunk
        for (int i = 0; i < chunk_genes; i++) {
            int global_i = start_idx + i;  // Global target index

            for (int j = 0; j < n_genes; j++) {
                // Extract I(r) curve for this pair
                arma::vec I_raw(n_radii);
                for (int r = 0; r < n_radii; r++) {
                    I_raw(r) = I_chunk(i, j, r);
                }

                // Smooth the curve
                arma::vec I_smooth = cpp_savgol_smooth(I_raw, smooth_window, smooth_poly);

                // Compute signed delta I
                SignatureResult sig = compute_signature_internal(I_smooth, radii);

                delta_I_matrix(global_i, j) = sig.delta_I_signed;
            }
        }

        // Check for user interrupt
        Rcpp::checkUserInterrupt();
    }

    if (verbose) Rcpp::Rcout << "Done!" << std::endl;

    return Rcpp::List::create(
        Rcpp::Named("delta_I_signed") = delta_I_matrix,
        Rcpp::Named("n_genes") = n_genes,
        Rcpp::Named("n_spots") = n_spots,
        Rcpp::Named("n_radii") = n_radii,
        Rcpp::Named("radii") = radii
    );
}

// =============================================================================
// Compute full delta I matrix using I_ND measure (optimized, no filtering)
//
// When all spots are used as both senders and receivers (no expression filtering),
// we can precompute weight matrices and use efficient matrix operations.
//
// I_ND (cosine similarity) = dot(z_f, W*z_g) / (||z_f|| * ||W*z_g||)
//
// We precompute:
// - Z: standardized expression matrix (n_genes x n_spots)
// - ZW: Z * W for each radius (spatial lag of all genes)
// - norms_z: ||z_i|| for all genes
// - norms_Wz: ||W*z_j|| for all genes at each radius
//
// @param expr_matrix Gene expression matrix (genes x spots), normalized
// @param coords Spot coordinates (Nx2 matrix)
// @param radii Distance bin edges
// @param chunk_size Number of target genes per chunk
// @param smooth_window Savitzky-Golay window size
// @param smooth_poly Savitzky-Golay polynomial order
// @param verbose Print progress messages
//
// @return List with delta_I_signed matrix
// =============================================================================
// [[Rcpp::export]]
Rcpp::List cpp_compute_delta_I_matrix_IND_fast(
    const arma::mat& expr_matrix,
    const arma::mat& coords,
    const arma::vec& radii,
    int chunk_size = 1000,
    int smooth_window = 5,
    int smooth_poly = 2,
    bool verbose = true
) {
    int n_genes = expr_matrix.n_rows;
    int n_spots = expr_matrix.n_cols;
    int n_radii = radii.n_elem;

    if (verbose) {
        Rcpp::Rcout << "Computing full delta I matrix (I_ND mode, optimized)" << std::endl;
        Rcpp::Rcout << n_spots << " spots, " << n_genes << " genes, "
                    << n_radii << " distance bins" << std::endl;
        Rcpp::Rcout << "Chunk size: " << chunk_size << " genes" << std::endl;
    }

    // Step 1: Standardize all genes
    if (verbose) Rcpp::Rcout << "Step 1: Standardizing expression matrix..." << std::endl;
    arma::mat Z(n_genes, n_spots);
    arma::vec norms_z(n_genes);  // ||z_i|| for each gene

    for (int g = 0; g < n_genes; g++) {
        arma::vec gene_expr = expr_matrix.row(g).t();
        arma::vec z = standardize_vec(gene_expr);
        Z.row(g) = z.t();
        norms_z(g) = arma::norm(z, 2);
    }

    // Step 2: Precompute weight matrices and ZW for all radii
    if (verbose) Rcpp::Rcout << "Step 2: Precomputing weight matrices and spatial lags..." << std::endl;

    arma::uvec all_idx = arma::regspace<arma::uvec>(0, n_spots - 1);

    std::vector<arma::mat> ZW_vec(n_radii);        // Z * W for each radius
    std::vector<arma::vec> norms_Wz_vec(n_radii);  // ||W*z_j|| for each gene at each radius

    for (int r = 0; r < n_radii; r++) {
        double r_inner = (r == 0) ? 0.0 : radii(r - 1);
        double r_outer = radii(r);

        // Create weight matrix (using directional mode = 1 for proper row normalization)
        arma::sp_mat W = cpp_create_ring_weight_matrix(
            coords, r_inner, r_outer, all_idx, all_idx, 1, true
        );

        if (verbose) {
            Rcpp::Rcout << "  Ring " << (r+1) << ": [" << r_inner << ", " << r_outer << ") - "
                        << W.n_nonzero << " non-zero weights" << std::endl;
        }

        // Compute ZW = Z * W (n_genes x n_spots)
        ZW_vec[r] = Z * W;

        // Compute ||W*z_j|| for each gene (column norms of ZW')
        // ZW[j,:] = z_j * W, so ||W*z_j|| = norm of row j of ZW
        arma::vec norms_Wz(n_genes);
        for (int g = 0; g < n_genes; g++) {
            norms_Wz(g) = arma::norm(ZW_vec[r].row(g), 2);
        }
        norms_Wz_vec[r] = norms_Wz;
    }

    // Step 3: Allocate result matrix
    if (verbose) Rcpp::Rcout << "Step 3: Allocating result matrix..." << std::endl;
    arma::mat delta_I_matrix(n_genes, n_genes);

    // Step 4: Process in chunks
    int n_chunks = (n_genes + chunk_size - 1) / chunk_size;
    if (verbose) Rcpp::Rcout << "Step 4: Processing " << n_chunks << " chunks..." << std::endl;

    for (int chunk = 0; chunk < n_chunks; chunk++) {
        int start_idx = chunk * chunk_size;
        int end_idx = std::min(start_idx + chunk_size, n_genes);
        int chunk_genes = end_idx - start_idx;

        if (verbose) {
            Rcpp::Rcout << "  Chunk " << (chunk + 1) << "/" << n_chunks
                        << " (genes " << start_idx << "-" << (end_idx - 1) << ")..." << std::endl;
        }

        // Extract Z_chunk for target genes
        arma::mat Z_chunk = Z.rows(start_idx, end_idx - 1);
        arma::vec norms_z_chunk = norms_z.subvec(start_idx, end_idx - 1);

        // Compute I_ND for this chunk at all radii
        arma::cube I_chunk(chunk_genes, n_genes, n_radii);

        for (int r = 0; r < n_radii; r++) {
            // Numerator: Z_chunk * ZW_vec[r]' = chunk_genes x n_genes
            arma::mat numerator = Z_chunk * ZW_vec[r].t();

            // Denominator: norms_z_chunk[i] * norms_Wz_vec[r][j]
            // = outer product of norms_z_chunk and norms_Wz_vec[r]
            arma::mat denominator = norms_z_chunk * norms_Wz_vec[r].t();

            // I_ND = numerator / denominator (element-wise)
            // Handle division by zero
            arma::mat I_r = numerator / (denominator + 1e-10);

            // Set to NA where denominator is too small
            for (int i = 0; i < chunk_genes; i++) {
                for (int j = 0; j < n_genes; j++) {
                    if (denominator(i, j) < 1e-10) {
                        I_r(i, j) = NA_REAL;
                    }
                }
            }

            I_chunk.slice(r) = I_r;
        }

        // Smooth and compute delta I for this chunk
        for (int i = 0; i < chunk_genes; i++) {
            int global_i = start_idx + i;

            for (int j = 0; j < n_genes; j++) {
                arma::vec I_raw(n_radii);
                for (int r = 0; r < n_radii; r++) {
                    I_raw(r) = I_chunk(i, j, r);
                }

                arma::vec I_smooth = cpp_savgol_smooth(I_raw, smooth_window, smooth_poly);
                SignatureResult sig = compute_signature_internal(I_smooth, radii);

                delta_I_matrix(global_i, j) = sig.delta_I_signed;
            }
        }

        Rcpp::checkUserInterrupt();
    }

    if (verbose) Rcpp::Rcout << "Done!" << std::endl;

    return Rcpp::List::create(
        Rcpp::Named("delta_I_signed") = delta_I_matrix,
        Rcpp::Named("n_genes") = n_genes,
        Rcpp::Named("n_spots") = n_spots,
        Rcpp::Named("n_radii") = n_radii,
        Rcpp::Named("radii") = radii,
        Rcpp::Named("mode") = "directional_IND_fast"
    );
}

// =============================================================================
// Compute full delta I matrix using directional I_ND measure (with filtering)
//
// For directional mode with expression-based filtering, each factor gene has
// different sender/receiver spots, so we cannot precompute weight matrices
// across all factors. This function loops through factors.
//
// @param expr_matrix Gene expression matrix (genes x spots), normalized
// @param coords Spot coordinates (Nx2 matrix)
// @param radii Distance bin edges
// @param sender_percentile Percentile threshold for sender definition
// @param smooth_window Savitzky-Golay window size
// @param smooth_poly Savitzky-Golay polynomial order
// @param verbose Print progress messages
//
// @return List with delta_I_signed matrix
// =============================================================================
// [[Rcpp::export]]
Rcpp::List cpp_compute_delta_I_matrix_directional(
    const arma::mat& expr_matrix,
    const arma::mat& coords,
    const arma::vec& radii,
    double sender_percentile = 75.0,
    int smooth_window = 5,
    int smooth_poly = 2,
    bool verbose = true
) {
    int n_genes = expr_matrix.n_rows;
    int n_spots = expr_matrix.n_cols;
    int n_radii = radii.n_elem;

    if (verbose) {
        Rcpp::Rcout << "Computing full delta I matrix (directional I_ND mode)" << std::endl;
        Rcpp::Rcout << n_spots << " spots, " << n_genes << " genes, "
                    << n_radii << " distance bins" << std::endl;
        Rcpp::Rcout << "Sender percentile: " << sender_percentile << std::endl;
    }

    // Allocate result matrix
    arma::mat delta_I_matrix(n_genes, n_genes);

    // Process each factor gene
    for (int f = 0; f < n_genes; f++) {
        if (verbose && ((f + 1) % 500 == 0 || f == 0 || f == n_genes - 1)) {
            Rcpp::Rcout << "  Processing factor " << (f + 1) << " / " << n_genes << std::endl;
        }

        // Get factor expression and define sender/receiver split
        arma::vec factor_expr = expr_matrix.row(f).t();
        Rcpp::List split = cpp_define_sender_receiver(factor_expr, sender_percentile);

        arma::uvec s_idx_r = split["sender_idx"];
        arma::uvec r_idx_r = split["receiver_idx"];
        arma::uvec sender_idx = s_idx_r - 1;  // Convert to 0-based
        arma::uvec receiver_idx = r_idx_r - 1;

        int n_senders = sender_idx.n_elem;
        int n_receivers = receiver_idx.n_elem;

        if (n_senders == 0 || n_receivers == 0) {
            // No valid split - fill with NA
            for (int g = 0; g < n_genes; g++) {
                delta_I_matrix(g, f) = NA_REAL;
            }
            continue;
        }

        // Standardize factor expression for senders
        arma::vec factor_sub = factor_expr.elem(sender_idx);
        arma::vec z_f = standardize_vec(factor_sub);

        // Precompute weight matrices for this factor's sender/receiver split
        std::vector<arma::sp_mat> W_vec(n_radii);
        for (int r = 0; r < n_radii; r++) {
            double r_inner = (r == 0) ? 0.0 : radii(r - 1);
            double r_outer = radii(r);

            W_vec[r] = cpp_create_ring_weight_matrix(
                coords, r_inner, r_outer, sender_idx, receiver_idx, 1, true
            );
        }

        // Compute I_ND for all target genes against this factor
        for (int g = 0; g < n_genes; g++) {
            // Get target gene expression for receivers
            arma::vec gene_expr = expr_matrix.row(g).t();
            arma::vec gene_sub = gene_expr.elem(receiver_idx);
            arma::vec z_g = standardize_vec(gene_sub);

            // Compute I_ND curve
            arma::vec I_raw(n_radii);
            for (int r = 0; r < n_radii; r++) {
                I_raw(r) = cpp_compute_spatial_correlation(z_f, z_g, W_vec[r], 1);
            }

            // Smooth and compute delta I
            arma::vec I_smooth = cpp_savgol_smooth(I_raw, smooth_window, smooth_poly);
            SignatureResult sig = compute_signature_internal(I_smooth, radii);

            delta_I_matrix(g, f) = sig.delta_I_signed;
        }

        // Check for user interrupt
        if (f % 100 == 0) {
            Rcpp::checkUserInterrupt();
        }
    }

    if (verbose) Rcpp::Rcout << "Done!" << std::endl;

    return Rcpp::List::create(
        Rcpp::Named("delta_I_signed") = delta_I_matrix,
        Rcpp::Named("n_genes") = n_genes,
        Rcpp::Named("n_spots") = n_spots,
        Rcpp::Named("n_radii") = n_radii,
        Rcpp::Named("radii") = radii,
        Rcpp::Named("mode") = "directional",
        Rcpp::Named("sender_percentile") = sender_percentile
    );
}

// =============================================================================
// Compute signatures for multiple factors (convenience wrapper)
//
// Runs batch signature computation for multiple factor genes at once.
//
// @param expr_matrix Gene expression matrix (genes x spots), normalized
// @param gene_names Character vector of gene names
// @param factor_indices 0-based indices of factor genes
// @param coords Spot coordinates (Nx2 matrix)
// @param radii Distance bin edges
// @param mode 0 = BIVARIATE, 1 = DIRECTIONAL
// @param sender_percentile Percentile for sender definition
// @param smooth_window Savitzky-Golay window size
// @param smooth_poly Savitzky-Golay polynomial order
// @param verbose Print progress messages
//
// @return List of DataFrames, one per factor
// =============================================================================
// [[Rcpp::export]]
Rcpp::List cpp_compute_signatures_multi_factor(
    const arma::mat& expr_matrix,
    const Rcpp::CharacterVector& gene_names,
    const arma::ivec& factor_indices,
    const arma::mat& coords,
    const arma::vec& radii,
    int mode = 0,
    double sender_percentile = 75.0,
    int smooth_window = 5,
    int smooth_poly = 2,
    bool verbose = true
) {
    int n_factors = factor_indices.n_elem;
    Rcpp::List results(n_factors);

    for (int f = 0; f < n_factors; f++) {
        if (verbose) {
            Rcpp::Rcout << "\n=== Factor " << (f + 1) << " / " << n_factors
                        << " (gene index " << factor_indices(f) << ") ===" << std::endl;
        }

        results[f] = cpp_compute_all_signatures(
            expr_matrix,
            gene_names,
            factor_indices(f),
            coords,
            radii,
            mode,
            sender_percentile,
            smooth_window,
            smooth_poly,
            verbose
        );
    }

    return results;
}
