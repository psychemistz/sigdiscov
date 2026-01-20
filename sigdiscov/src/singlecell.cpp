// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <vector>
#include <string>

using namespace Rcpp;
using namespace arma;

// =============================================================================
// Single-Cell Spatial Transcriptomics Analysis
//
// This module implements cell-type-aware I_ND (Normalized Directional Moran's I)
// computation for single-cell resolution spatial transcriptomics data (CosMx,
// Xenium, MERFISH, etc.).
//
// Key difference from Visium:
// - Sender/Receiver defined by CELL TYPE, not expression level
// - Physical coordinates (x, y in micrometers), not grid coordinates
// - Finer distance resolution (10-200 um vs 100-600 um)
// =============================================================================

// =============================================================================
// Helper: Compute Euclidean distance between two points
// =============================================================================
inline double sc_euclidean_dist(double x1, double y1, double x2, double y2) {
    double dx = x1 - x2;
    double dy = y1 - y2;
    return std::sqrt(dx * dx + dy * dy);
}

// =============================================================================
// Helper: Standardize a vector (z-score normalization)
// =============================================================================
arma::vec sc_standardize_vec(const arma::vec& x) {
    double mean_x = arma::mean(x);
    double std_x = arma::stddev(x, 1);  // 1 = normalize by N (population)

    if (std_x < 1e-10) {
        return arma::zeros<arma::vec>(x.n_elem);
    }

    return (x - mean_x) / std_x;
}

// =============================================================================
// Create ring weight matrix for single-cell ST with cell-type filtering
//
// Creates an n_senders x n_receivers weight matrix where entries are non-zero
// only for sender-receiver pairs within the specified distance ring.
//
// @param coords Nx2 matrix of ALL cell coordinates (x, y in micrometers)
// @param sender_idx Indices of sender cells (0-based, pointing to coords)
// @param receiver_idx Indices of receiver cells (0-based, pointing to coords)
// @param r_inner Inner radius of ring (inclusive)
// @param r_outer Outer radius of ring (exclusive)
// @param row_normalize Whether to row-normalize the weights (default: true)
//
// @return Sparse weight matrix (n_senders x n_receivers)
// =============================================================================
// [[Rcpp::export]]
arma::sp_mat cpp_create_ring_weight_matrix_sc(
    const arma::mat& coords,
    const arma::ivec& sender_idx,
    const arma::ivec& receiver_idx,
    double r_inner,
    double r_outer,
    bool row_normalize = true
) {
    int n_senders = sender_idx.n_elem;
    int n_receivers = receiver_idx.n_elem;

    std::vector<arma::uword> row_indices;
    std::vector<arma::uword> col_indices;
    std::vector<double> values;

    // First pass: find neighbors and compute row sums
    std::vector<std::vector<int>> neighbors(n_senders);
    std::vector<double> row_sums(n_senders, 0.0);

    for (int i = 0; i < n_senders; i++) {
        int si = sender_idx(i);  // Actual cell index in coords
        double xi = coords(si, 0);
        double yi = coords(si, 1);

        for (int j = 0; j < n_receivers; j++) {
            int rj = receiver_idx(j);  // Actual cell index in coords

            // Allow self-links if same cell (when sender_type == receiver_type)
            // but skip if same physical cell
            if (si == rj) continue;

            double dist = sc_euclidean_dist(xi, yi, coords(rj, 0), coords(rj, 1));

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

// =============================================================================
// Compute I_ND (cosine similarity) for single-cell ST
//
// Computes the Normalized Directional Moran's I between factor expression
// in sender cells and gene expression in receiver cells.
//
// Formula:
//              z_U^f · (W × z_V^g)
// I_ND = ─────────────────────────────
//        ||z_U^f|| × ||W × z_V^g||
//
// @param z_f Standardized factor expression in sender cells (length n_senders)
// @param z_g Standardized gene expression in receiver cells (length n_receivers)
// @param W Row-normalized weight matrix (n_senders x n_receivers)
//
// @return Scalar I_ND value in [-1, 1]
// =============================================================================
double compute_IND_sc(
    const arma::vec& z_f,
    const arma::vec& z_g,
    const arma::sp_mat& W
) {
    // Compute spatial lag: weighted average of receiver gene expression
    arma::vec lag_g = W * z_g;  // length n_senders

    // Cosine similarity (I_ND formula)
    double numerator = arma::dot(z_f, lag_g);
    double norm_zf = arma::norm(z_f, 2);
    double norm_lag = arma::norm(lag_g, 2);
    double denominator = norm_zf * norm_lag;

    if (denominator < 1e-10) {
        return NA_REAL;
    }

    return numerator / denominator;
}

// =============================================================================
// Compute I_ND curve across multiple distance rings for single-cell ST
//
// @param factor_expr Factor expression in sender cells (length n_senders)
// @param gene_expr Gene expression in receiver cells (length n_receivers)
// @param coords ALL cell coordinates (Nx2 matrix)
// @param sender_idx Indices of sender cells (0-based)
// @param receiver_idx Indices of receiver cells (0-based)
// @param radii Distance bin outer edges (rings: [0, r1), [r1, r2), ...)
// @param smooth_window Smoothing window size (must be odd)
// @param smooth_poly Polynomial order for smoothing
//
// @return List with I_raw, I_smooth, and signature metrics
// =============================================================================
// [[Rcpp::export]]
Rcpp::List cpp_compute_IND_curve_sc(
    const arma::vec& factor_expr,
    const arma::vec& gene_expr,
    const arma::mat& coords,
    const arma::ivec& sender_idx,
    const arma::ivec& receiver_idx,
    const arma::vec& radii,
    int smooth_window = 5,
    int smooth_poly = 2
) {
    int n_radii = radii.n_elem;
    arma::vec I_raw(n_radii);

    // Standardize expressions
    arma::vec z_f = sc_standardize_vec(factor_expr);
    arma::vec z_g = sc_standardize_vec(gene_expr);

    // Compute I_ND at each radius
    for (int r = 0; r < n_radii; r++) {
        double r_inner = (r == 0) ? 0.0 : radii(r - 1);
        double r_outer = radii(r);

        arma::sp_mat W = cpp_create_ring_weight_matrix_sc(
            coords, sender_idx, receiver_idx,
            r_inner, r_outer, true
        );

        I_raw(r) = compute_IND_sc(z_f, z_g, W);
    }

    // Smooth the curve
    arma::vec I_smooth = I_raw;

    // Simple moving average smoothing for spike removal
    if (n_radii >= 3) {
        if (smooth_window % 2 == 0) smooth_window++;
        if (smooth_window > n_radii) {
            smooth_window = (n_radii % 2 == 1) ? n_radii : n_radii - 1;
        }

        int half_window = smooth_window / 2;

        for (int i = 0; i < n_radii; i++) {
            double sum = 0.0;
            int count = 0;

            for (int j = -half_window; j <= half_window; j++) {
                int idx = i + j;
                if (idx < 0) idx = -idx;
                if (idx >= n_radii) idx = 2 * n_radii - idx - 2;

                if (idx >= 0 && idx < n_radii && !R_IsNA(I_raw(idx))) {
                    sum += I_raw(idx);
                    count++;
                }
            }

            I_smooth(i) = (count > 0) ? sum / count : I_raw(i);
        }
    }

    // Compute signature metrics
    double I_max = NA_REAL, I_min = NA_REAL;
    double I_short = NA_REAL, I_long = NA_REAL;
    int valid_count = 0;

    for (int i = 0; i < n_radii; i++) {
        if (!R_IsNA(I_smooth(i))) {
            if (valid_count == 0) {
                I_max = I_min = I_smooth(i);
                I_short = I_smooth(i);
            } else {
                if (I_smooth(i) > I_max) I_max = I_smooth(i);
                if (I_smooth(i) < I_min) I_min = I_smooth(i);
            }
            I_long = I_smooth(i);
            valid_count++;
        }
    }

    double delta_I = NA_REAL;
    double delta_I_signed = NA_REAL;
    int sign = NA_INTEGER;

    if (valid_count > 0) {
        delta_I = I_max - I_min;
        double trend = I_short - I_long;
        sign = (trend >= 0) ? 1 : -1;
        delta_I_signed = sign * delta_I;
    }

    return Rcpp::List::create(
        Rcpp::Named("radii") = radii,
        Rcpp::Named("I_raw") = I_raw,
        Rcpp::Named("I_smooth") = I_smooth,
        Rcpp::Named("delta_I") = delta_I,
        Rcpp::Named("delta_I_signed") = delta_I_signed,
        Rcpp::Named("sign") = sign,
        Rcpp::Named("I_max") = I_max,
        Rcpp::Named("I_min") = I_min,
        Rcpp::Named("I_short") = I_short,
        Rcpp::Named("I_long") = I_long
    );
}

// =============================================================================
// Precompute weight matrices for all radii
//
// Caches weight matrices to avoid recomputation when processing multiple genes.
// =============================================================================
std::vector<arma::sp_mat> precompute_weight_matrices_sc(
    const arma::mat& coords,
    const arma::ivec& sender_idx,
    const arma::ivec& receiver_idx,
    const arma::vec& radii
) {
    int n_radii = radii.n_elem;
    std::vector<arma::sp_mat> W_list;
    W_list.reserve(n_radii);

    for (int r = 0; r < n_radii; r++) {
        double r_inner = (r == 0) ? 0.0 : radii(r - 1);
        double r_outer = radii(r);

        W_list.push_back(cpp_create_ring_weight_matrix_sc(
            coords, sender_idx, receiver_idx,
            r_inner, r_outer, true
        ));
    }

    return W_list;
}

// =============================================================================
// Compute I_ND signatures for all genes (single-cell ST)
//
// Main batch processing function for computing I_ND signatures for all genes
// in the expression matrix against a specified factor gene, using cell-type
// based sender and receiver definitions.
//
// @param expr_matrix Gene expression matrix (genes x cells)
// @param gene_names Gene name vector
// @param factor_idx Index of factor gene (0-based)
// @param coords ALL cell coordinates (Nx2 matrix)
// @param sender_idx Indices of sender cells (0-based)
// @param receiver_idx Indices of receiver cells (0-based)
// @param radii Distance bin outer edges
// @param smooth_window Smoothing window size
// @param smooth_poly Polynomial order
// @param verbose Print progress
//
// @return DataFrame with I_ND signatures for all genes
// =============================================================================
// [[Rcpp::export]]
Rcpp::DataFrame cpp_compute_IND_signatures_sc(
    const arma::mat& expr_matrix,
    const Rcpp::CharacterVector& gene_names,
    int factor_idx,
    const arma::mat& coords,
    const arma::ivec& sender_idx,
    const arma::ivec& receiver_idx,
    const arma::vec& radii,
    int smooth_window = 5,
    int smooth_poly = 2,
    bool verbose = true
) {
    int n_genes = expr_matrix.n_rows;
    int n_radii = radii.n_elem;
    int n_senders = sender_idx.n_elem;
    int n_receivers = receiver_idx.n_elem;

    if (verbose) {
        Rcpp::Rcout << "  Precomputing weight matrices..." << std::endl;
    }

    // Precompute weight matrices for all radii
    std::vector<arma::sp_mat> W_list = precompute_weight_matrices_sc(
        coords, sender_idx, receiver_idx, radii
    );

    // Extract factor expression in sender cells
    arma::vec factor_expr(n_senders);
    for (int i = 0; i < n_senders; i++) {
        factor_expr(i) = expr_matrix(factor_idx, sender_idx(i));
    }
    arma::vec z_f = sc_standardize_vec(factor_expr);

    // Result storage
    Rcpp::CharacterVector out_genes(n_genes);
    Rcpp::NumericVector out_IND_r1(n_genes);
    Rcpp::NumericVector out_IND_max(n_genes);
    Rcpp::NumericVector out_delta_I(n_genes);

    // I_ND at each radius
    Rcpp::NumericMatrix out_IND_radii(n_genes, n_radii);

    if (verbose) {
        Rcpp::Rcout << "  Computing signatures for " << n_genes << " genes..." << std::endl;
    }

    // Process each gene
    for (int g = 0; g < n_genes; g++) {
        if (verbose && (g + 1) % 500 == 0) {
            Rcpp::Rcout << "    " << (g + 1) << " / " << n_genes << std::endl;
        }

        // Check for user interrupt
        if (g % 100 == 0) {
            Rcpp::checkUserInterrupt();
        }

        // Extract gene expression in receiver cells
        arma::vec gene_expr(n_receivers);
        for (int j = 0; j < n_receivers; j++) {
            gene_expr(j) = expr_matrix(g, receiver_idx(j));
        }
        arma::vec z_g = sc_standardize_vec(gene_expr);

        // Compute I_ND at each radius
        arma::vec I_curve(n_radii);
        for (int r = 0; r < n_radii; r++) {
            I_curve(r) = compute_IND_sc(z_f, z_g, W_list[r]);
            out_IND_radii(g, r) = I_curve(r);
        }

        // Apply simple smoothing
        arma::vec I_smooth = I_curve;
        if (n_radii >= 3) {
            int sw = smooth_window;
            if (sw % 2 == 0) sw++;
            if (sw > n_radii) sw = (n_radii % 2 == 1) ? n_radii : n_radii - 1;

            int half = sw / 2;
            for (int i = 0; i < n_radii; i++) {
                double sum = 0.0;
                int count = 0;
                for (int j = -half; j <= half; j++) {
                    int idx = i + j;
                    if (idx < 0) idx = -idx;
                    if (idx >= n_radii) idx = 2 * n_radii - idx - 2;
                    if (idx >= 0 && idx < n_radii && !R_IsNA(I_curve(idx))) {
                        sum += I_curve(idx);
                        count++;
                    }
                }
                I_smooth(i) = (count > 0) ? sum / count : I_curve(i);
            }
        }

        // Compute summary statistics from smoothed curve
        double I_max = NA_REAL, I_min = NA_REAL;
        double I_short = NA_REAL, I_long = NA_REAL;
        int valid_count = 0;

        for (int r = 0; r < n_radii; r++) {
            if (!R_IsNA(I_smooth(r))) {
                if (valid_count == 0) {
                    I_max = I_min = I_smooth(r);
                    I_short = I_smooth(r);
                } else {
                    if (I_smooth(r) > I_max) I_max = I_smooth(r);
                    if (I_smooth(r) < I_min) I_min = I_smooth(r);
                }
                I_long = I_smooth(r);
                valid_count++;
            }
        }

        out_genes[g] = gene_names[g];

        if (valid_count > 0) {
            out_IND_r1[g] = I_smooth(0);
            out_IND_max[g] = I_max;

            double delta = I_max - I_min;
            double trend = I_short - I_long;
            int sign = (trend >= 0) ? 1 : -1;
            out_delta_I[g] = sign * delta;
        } else {
            out_IND_r1[g] = NA_REAL;
            out_IND_max[g] = NA_REAL;
            out_delta_I[g] = NA_REAL;
        }
    }

    if (verbose) {
        Rcpp::Rcout << "  Done!" << std::endl;
    }

    // Build result DataFrame with all columns at once
    // First create a List, then convert to DataFrame
    Rcpp::List result_list;
    result_list["gene"] = out_genes;
    result_list["I_ND_r1"] = out_IND_r1;
    result_list["I_ND_max"] = out_IND_max;
    result_list["delta_I"] = out_delta_I;

    // Add I_ND columns for each radius
    for (int r = 0; r < n_radii; r++) {
        std::string col_name = "I_ND_" + std::to_string(static_cast<int>(radii(r)));
        result_list[col_name] = out_IND_radii(Rcpp::_, r);
    }

    // Convert to DataFrame
    Rcpp::DataFrame result(result_list);
    result.attr("class") = Rcpp::CharacterVector::create("data.frame");
    result.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER, -n_genes);

    return result;
}

// =============================================================================
// Compute I_ND signatures for multiple factors with precomputed weights
//
// OPTIMIZED VERSION: Precomputes weight matrices once per (sender, receiver)
// pair and reuses them across all factors. This is much faster when computing
// signatures for many factors.
//
// @param expr_matrix Gene expression matrix (genes x cells)
// @param gene_names Gene name vector
// @param factor_indices Vector of factor gene indices (0-based)
// @param factor_names Names of factor genes
// @param coords ALL cell coordinates (Nx2 matrix)
// @param sender_idx Indices of sender cells (0-based)
// @param receiver_idx Indices of receiver cells (0-based)
// @param radii Distance bin outer edges
// @param smooth_window Smoothing window size
// @param smooth_poly Polynomial order
// @param verbose Print progress
//
// @return List of DataFrames, one per factor
// =============================================================================
// [[Rcpp::export]]
Rcpp::List cpp_compute_IND_multi_factor(
    const arma::mat& expr_matrix,
    const Rcpp::CharacterVector& gene_names,
    const arma::ivec& factor_indices,
    const Rcpp::CharacterVector& factor_names,
    const arma::mat& coords,
    const arma::ivec& sender_idx,
    const arma::ivec& receiver_idx,
    const arma::vec& radii,
    int smooth_window = 5,
    int smooth_poly = 2,
    bool verbose = true
) {
    int n_genes = expr_matrix.n_rows;
    int n_factors = factor_indices.n_elem;
    int n_radii = radii.n_elem;
    int n_senders = sender_idx.n_elem;
    int n_receivers = receiver_idx.n_elem;

    if (verbose) {
        Rcpp::Rcout << "  Precomputing weight matrices (once for all factors)..." << std::endl;
    }

    // Precompute weight matrices ONCE - this is the key optimization
    std::vector<arma::sp_mat> W_list = precompute_weight_matrices_sc(
        coords, sender_idx, receiver_idx, radii
    );

    // Pre-extract and standardize all receiver gene expressions
    // This avoids redundant extraction for each factor
    if (verbose) {
        Rcpp::Rcout << "  Pre-extracting receiver expressions..." << std::endl;
    }

    arma::mat receiver_expr(n_genes, n_receivers);
    for (int g = 0; g < n_genes; g++) {
        for (int j = 0; j < n_receivers; j++) {
            receiver_expr(g, j) = expr_matrix(g, receiver_idx(j));
        }
    }

    // Standardize all receiver expressions
    arma::mat z_receiver(n_genes, n_receivers);
    for (int g = 0; g < n_genes; g++) {
        arma::vec expr_row = receiver_expr.row(g).t();
        z_receiver.row(g) = sc_standardize_vec(expr_row).t();
    }

    // Precompute spatial lags for all genes at all radii
    // lag[r][g] = W_list[r] * z_receiver.row(g).t()
    if (verbose) {
        Rcpp::Rcout << "  Precomputing spatial lags for all genes..." << std::endl;
    }

    std::vector<arma::mat> lag_matrices(n_radii);
    for (int r = 0; r < n_radii; r++) {
        // Each row is the spatial lag for gene g: W * z_g
        lag_matrices[r].set_size(n_genes, n_senders);
        for (int g = 0; g < n_genes; g++) {
            arma::vec z_g = z_receiver.row(g).t();
            arma::vec lag_g = W_list[r] * z_g;
            lag_matrices[r].row(g) = lag_g.t();
        }
    }

    // Now process each factor - only need to extract factor expression
    Rcpp::List results(n_factors);

    for (int f = 0; f < n_factors; f++) {
        if (verbose) {
            Rcpp::Rcout << "  Processing factor " << (f + 1) << "/" << n_factors
                        << ": " << factor_names[f] << std::endl;
        }

        int factor_idx = factor_indices(f);

        // Extract and standardize factor expression in senders
        arma::vec factor_expr(n_senders);
        for (int i = 0; i < n_senders; i++) {
            factor_expr(i) = expr_matrix(factor_idx, sender_idx(i));
        }
        arma::vec z_f = sc_standardize_vec(factor_expr);
        double norm_zf = arma::norm(z_f, 2);

        // Result storage for this factor
        Rcpp::CharacterVector out_genes(n_genes);
        Rcpp::NumericVector out_IND_r1(n_genes);
        Rcpp::NumericVector out_IND_max(n_genes);
        Rcpp::NumericVector out_delta_I(n_genes);
        Rcpp::NumericMatrix out_IND_radii(n_genes, n_radii);

        // Process each gene
        for (int g = 0; g < n_genes; g++) {
            if (g % 100 == 0) {
                Rcpp::checkUserInterrupt();
            }

            // Compute I_ND at each radius using precomputed lags
            arma::vec I_curve(n_radii);
            for (int r = 0; r < n_radii; r++) {
                arma::vec lag_g = lag_matrices[r].row(g).t();
                double norm_lag = arma::norm(lag_g, 2);
                double denominator = norm_zf * norm_lag;

                if (denominator < 1e-10) {
                    I_curve(r) = NA_REAL;
                } else {
                    I_curve(r) = arma::dot(z_f, lag_g) / denominator;
                }
                out_IND_radii(g, r) = I_curve(r);
            }

            // Apply smoothing
            arma::vec I_smooth = I_curve;
            if (n_radii >= 3) {
                int sw = smooth_window;
                if (sw % 2 == 0) sw++;
                if (sw > n_radii) sw = (n_radii % 2 == 1) ? n_radii : n_radii - 1;

                int half = sw / 2;
                for (int i = 0; i < n_radii; i++) {
                    double sum = 0.0;
                    int count = 0;
                    for (int j = -half; j <= half; j++) {
                        int idx = i + j;
                        if (idx < 0) idx = -idx;
                        if (idx >= n_radii) idx = 2 * n_radii - idx - 2;
                        if (idx >= 0 && idx < n_radii && !R_IsNA(I_curve(idx))) {
                            sum += I_curve(idx);
                            count++;
                        }
                    }
                    I_smooth(i) = (count > 0) ? sum / count : I_curve(i);
                }
            }

            // Compute summary statistics
            double I_max = NA_REAL, I_min = NA_REAL;
            double I_short = NA_REAL, I_long = NA_REAL;
            int valid_count = 0;

            for (int r = 0; r < n_radii; r++) {
                if (!R_IsNA(I_smooth(r))) {
                    if (valid_count == 0) {
                        I_max = I_min = I_smooth(r);
                        I_short = I_smooth(r);
                    } else {
                        if (I_smooth(r) > I_max) I_max = I_smooth(r);
                        if (I_smooth(r) < I_min) I_min = I_smooth(r);
                    }
                    I_long = I_smooth(r);
                    valid_count++;
                }
            }

            out_genes[g] = gene_names[g];

            if (valid_count > 0) {
                out_IND_r1[g] = I_smooth(0);
                out_IND_max[g] = I_max;
                double delta = I_max - I_min;
                double trend = I_short - I_long;
                int sign = (trend >= 0) ? 1 : -1;
                out_delta_I[g] = sign * delta;
            } else {
                out_IND_r1[g] = NA_REAL;
                out_IND_max[g] = NA_REAL;
                out_delta_I[g] = NA_REAL;
            }
        }

        // Build DataFrame for this factor
        Rcpp::List result_list;
        result_list["gene"] = out_genes;
        result_list["I_ND_r1"] = out_IND_r1;
        result_list["I_ND_max"] = out_IND_max;
        result_list["delta_I"] = out_delta_I;

        for (int r = 0; r < n_radii; r++) {
            std::string col_name = "I_ND_" + std::to_string(static_cast<int>(radii(r)));
            result_list[col_name] = out_IND_radii(Rcpp::_, r);
        }

        Rcpp::DataFrame df(result_list);
        df.attr("class") = Rcpp::CharacterVector::create("data.frame");
        df.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER, -n_genes);

        results[f] = df;
    }

    results.attr("names") = factor_names;

    if (verbose) {
        Rcpp::Rcout << "  Done!" << std::endl;
    }

    return results;
}
