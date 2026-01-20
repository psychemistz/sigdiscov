// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <vector>

using namespace Rcpp;
using namespace arma;

// =============================================================================
// Permutation Tests for Spatial Correlation
//
// Efficient permutation tests for Moran's I and I_ND metrics.
// Key optimization: Compute spatial lag ONCE, permute factor expression.
// =============================================================================

// =============================================================================
// Helper: Standardize a vector (z-score)
// =============================================================================
arma::vec perm_standardize(const arma::vec& x) {
    double mean_x = arma::mean(x);
    double std_x = arma::stddev(x, 1);
    if (std_x < 1e-10) {
        return arma::zeros<arma::vec>(x.n_elem);
    }
    return (x - mean_x) / std_x;
}

// =============================================================================
// Permutation Test for Bivariate Moran's I (Single Gene)
//
// Tests H0: No spatial association between factor and gene expression
//
// @param z_f Standardized factor expression (length n)
// @param z_g Standardized gene expression (length n)
// @param W Weight matrix (n x n, sparse)
// @param n_perm Number of permutations
// @param seed Random seed (-1 for no seed)
//
// @return List with I_obs, p_value, null_mean, null_sd
// =============================================================================
// [[Rcpp::export]]
Rcpp::List cpp_permutation_test_moran(
    const arma::vec& z_f,
    const arma::vec& z_g,
    const arma::sp_mat& W,
    int n_perm = 999,
    int seed = -1
) {
    int n = z_f.n_elem;

    if (seed >= 0) {
        arma::arma_rng::set_seed(seed);
    }

    // Compute spatial lag ONCE (expensive operation)
    arma::vec lag_g = W * z_g;

    // Observed statistic: I = (z_f' * lag_g) / n
    double I_obs = arma::dot(z_f, lag_g) / n;

    // Permutation test
    int count_extreme = 0;
    arma::vec null_dist(n_perm);

    for (int b = 0; b < n_perm; b++) {
        // Random permutation of z_f
        arma::uvec perm_idx = arma::randperm(n);
        arma::vec z_f_perm = z_f.elem(perm_idx);

        // Permuted statistic (lag_g unchanged)
        double I_perm = arma::dot(z_f_perm, lag_g) / n;
        null_dist(b) = I_perm;

        // Count extreme (two-sided)
        if (std::abs(I_perm) >= std::abs(I_obs)) {
            count_extreme++;
        }
    }

    // P-value with continuity correction
    double p_value = (1.0 + count_extreme) / (n_perm + 1.0);

    return Rcpp::List::create(
        Rcpp::Named("I_obs") = I_obs,
        Rcpp::Named("p_value") = p_value,
        Rcpp::Named("null_mean") = arma::mean(null_dist),
        Rcpp::Named("null_sd") = arma::stddev(null_dist),
        Rcpp::Named("n_perm") = n_perm
    );
}

// =============================================================================
// Permutation Test for I_ND (Single Gene)
//
// Tests H0: No directional spatial association from sender to receiver
//
// @param z_U Standardized factor expression in senders (length n_senders)
// @param z_V Standardized gene expression in receivers (length n_receivers)
// @param W Row-normalized weight matrix (n_senders x n_receivers, sparse)
// @param n_perm Number of permutations
// @param seed Random seed
//
// @return List with I_ND_obs, p_value, null_mean, null_sd
// =============================================================================
// [[Rcpp::export]]
Rcpp::List cpp_permutation_test_IND(
    const arma::vec& z_U,
    const arma::vec& z_V,
    const arma::sp_mat& W,
    int n_perm = 999,
    int seed = -1
) {
    int n_senders = z_U.n_elem;

    if (seed >= 0) {
        arma::arma_rng::set_seed(seed);
    }

    // Compute spatial lag ONCE
    arma::vec lag_V = W * z_V;

    // Pre-compute constants
    double norm_lag = arma::norm(lag_V);
    double norm_U = arma::norm(z_U);

    if (norm_lag < 1e-10 || norm_U < 1e-10) {
        return Rcpp::List::create(
            Rcpp::Named("I_ND_obs") = NA_REAL,
            Rcpp::Named("p_value") = NA_REAL,
            Rcpp::Named("null_mean") = NA_REAL,
            Rcpp::Named("null_sd") = NA_REAL,
            Rcpp::Named("n_perm") = n_perm
        );
    }

    // Observed I_ND (cosine similarity)
    double I_ND_obs = arma::dot(z_U, lag_V) / (norm_U * norm_lag);

    // Permutation test
    int count_extreme = 0;
    arma::vec null_dist(n_perm);

    for (int b = 0; b < n_perm; b++) {
        // Permute sender factor expression
        arma::uvec perm_idx = arma::randperm(n_senders);
        arma::vec z_U_perm = z_U.elem(perm_idx);

        // For standardized data: norm(z_U_perm) == norm(z_U)
        // Permutation preserves sum of squares
        double I_ND_perm = arma::dot(z_U_perm, lag_V) / (norm_U * norm_lag);
        null_dist(b) = I_ND_perm;

        // Count extreme (two-sided)
        if (std::abs(I_ND_perm) >= std::abs(I_ND_obs)) {
            count_extreme++;
        }
    }

    double p_value = (1.0 + count_extreme) / (n_perm + 1.0);

    return Rcpp::List::create(
        Rcpp::Named("I_ND_obs") = I_ND_obs,
        Rcpp::Named("p_value") = p_value,
        Rcpp::Named("null_mean") = arma::mean(null_dist),
        Rcpp::Named("null_sd") = arma::stddev(null_dist),
        Rcpp::Named("n_perm") = n_perm
    );
}

// =============================================================================
// Batch Permutation Test for All Genes vs One Factor
//
// Efficient: generates permutation indices ONCE, reuses for all genes
//
// @param expr_matrix Genes x cells expression matrix
// @param gene_names Gene name vector
// @param factor_idx Factor gene index (0-based)
// @param W Weight matrix (sparse)
// @param sender_idx Indices of sender cells (0-based)
// @param receiver_idx Indices of receiver cells (0-based)
// @param n_perm Number of permutations
// @param metric 0 = Moran's I, 1 = I_ND
// @param seed Random seed
// @param verbose Print progress
//
// @return DataFrame with gene, I_obs, p_value, z_score
// =============================================================================
// [[Rcpp::export]]
Rcpp::DataFrame cpp_batch_permutation_test(
    const arma::mat& expr_matrix,
    const Rcpp::CharacterVector& gene_names,
    int factor_idx,
    const arma::sp_mat& W,
    const arma::uvec& sender_idx,
    const arma::uvec& receiver_idx,
    int n_perm = 999,
    int metric = 1,
    int seed = -1,
    bool verbose = true
) {
    int n_genes = expr_matrix.n_rows;
    int n_senders = sender_idx.n_elem;
    int n_receivers = receiver_idx.n_elem;

    if (seed >= 0) {
        arma::arma_rng::set_seed(seed);
    }

    // Extract factor expression in senders
    arma::vec factor_full = expr_matrix.row(factor_idx).t();
    arma::vec z_f(n_senders);
    for (arma::uword i = 0; i < sender_idx.n_elem; i++) {
        z_f(i) = factor_full(sender_idx(i));
    }
    z_f = perm_standardize(z_f);
    double norm_zf = arma::norm(z_f);

    if (norm_zf < 1e-10) {
        Rcpp::stop("Factor gene has zero variance in sender cells");
    }

    // Pre-generate ALL permutation indices (reused for all genes)
    if (verbose) {
        Rcpp::Rcout << "Generating " << n_perm << " permutations..." << std::endl;
    }

    arma::umat perm_indices(n_perm, n_senders);
    for (int b = 0; b < n_perm; b++) {
        perm_indices.row(b) = arma::randperm(n_senders).t();
    }

    // Result vectors
    Rcpp::NumericVector out_I_obs(n_genes);
    Rcpp::NumericVector out_p_value(n_genes);
    Rcpp::NumericVector out_z_score(n_genes);

    if (verbose) {
        Rcpp::Rcout << "Testing " << n_genes << " genes..." << std::endl;
    }

    // Process each gene
    for (int g = 0; g < n_genes; g++) {
        if (verbose && (g + 1) % 500 == 0) {
            Rcpp::Rcout << "  " << (g + 1) << " / " << n_genes << std::endl;
        }

        // Check for user interrupt
        if (g % 100 == 0) {
            Rcpp::checkUserInterrupt();
        }

        // Extract gene expression in receivers
        arma::vec gene_full = expr_matrix.row(g).t();
        arma::vec z_g(n_receivers);
        for (arma::uword j = 0; j < receiver_idx.n_elem; j++) {
            z_g(j) = gene_full(receiver_idx(j));
        }
        z_g = perm_standardize(z_g);

        // Compute spatial lag ONCE per gene
        arma::vec lag_g = W * z_g;
        double norm_lag = arma::norm(lag_g);

        if (norm_lag < 1e-10) {
            out_I_obs[g] = NA_REAL;
            out_p_value[g] = NA_REAL;
            out_z_score[g] = NA_REAL;
            continue;
        }

        // Observed statistic
        double I_obs;
        if (metric == 0) {
            // Moran's I
            I_obs = arma::dot(z_f, lag_g) / n_senders;
        } else {
            // I_ND (cosine similarity)
            I_obs = arma::dot(z_f, lag_g) / (norm_zf * norm_lag);
        }

        // Permutation test
        int count_extreme = 0;
        double sum_perm = 0.0;
        double sum_perm_sq = 0.0;

        for (int b = 0; b < n_perm; b++) {
            // Get permuted factor using pre-computed indices
            arma::vec z_f_perm(n_senders);
            for (int i = 0; i < n_senders; i++) {
                z_f_perm(i) = z_f(perm_indices(b, i));
            }

            // Compute permuted statistic
            double I_perm;
            if (metric == 0) {
                I_perm = arma::dot(z_f_perm, lag_g) / n_senders;
            } else {
                // For standardized data, norm is preserved under permutation
                I_perm = arma::dot(z_f_perm, lag_g) / (norm_zf * norm_lag);
            }

            // Accumulate for null distribution statistics
            sum_perm += I_perm;
            sum_perm_sq += I_perm * I_perm;

            // Count extreme (two-sided)
            if (std::abs(I_perm) >= std::abs(I_obs)) {
                count_extreme++;
            }
        }

        // P-value
        double p_value = (1.0 + count_extreme) / (n_perm + 1.0);

        // Z-score
        double null_mean = sum_perm / n_perm;
        double null_var = sum_perm_sq / n_perm - null_mean * null_mean;
        double null_sd = std::sqrt(std::max(null_var, 1e-10));
        double z_score = (I_obs - null_mean) / null_sd;

        out_I_obs[g] = I_obs;
        out_p_value[g] = p_value;
        out_z_score[g] = z_score;
    }

    if (verbose) {
        Rcpp::Rcout << "Done!" << std::endl;
    }

    return Rcpp::DataFrame::create(
        Rcpp::Named("gene") = gene_names,
        Rcpp::Named("I_obs") = out_I_obs,
        Rcpp::Named("p_value") = out_p_value,
        Rcpp::Named("z_score") = out_z_score
    );
}

// =============================================================================
// Batch Permutation Test with Precomputed Weight Matrices (Multi-Radius)
//
// For signature computation: test at multiple distance radii
//
// @param expr_matrix Genes x cells expression matrix
// @param gene_names Gene name vector
// @param factor_idx Factor gene index (0-based)
// @param coords Coordinates matrix (N x 2)
// @param sender_idx Sender cell indices (0-based)
// @param receiver_idx Receiver cell indices (0-based)
// @param radii Distance radii vector
// @param n_perm Number of permutations
// @param seed Random seed
// @param verbose Print progress
//
// @return DataFrame with gene, I_ND at r1, p_value, z_score, delta_I
// =============================================================================
// [[Rcpp::export]]
Rcpp::DataFrame cpp_batch_permutation_test_radii(
    const arma::mat& expr_matrix,
    const Rcpp::CharacterVector& gene_names,
    int factor_idx,
    const arma::mat& coords,
    const arma::uvec& sender_idx,
    const arma::uvec& receiver_idx,
    const arma::vec& radii,
    int n_perm = 999,
    int seed = -1,
    bool verbose = true
) {
    int n_genes = expr_matrix.n_rows;
    int n_senders = sender_idx.n_elem;
    int n_receivers = receiver_idx.n_elem;
    int n_radii = radii.n_elem;

    if (seed >= 0) {
        arma::arma_rng::set_seed(seed);
    }

    // Extract factor expression in senders
    arma::vec factor_full = expr_matrix.row(factor_idx).t();
    arma::vec z_f(n_senders);
    for (arma::uword i = 0; i < sender_idx.n_elem; i++) {
        z_f(i) = factor_full(sender_idx(i));
    }
    z_f = perm_standardize(z_f);
    double norm_zf = arma::norm(z_f);

    if (norm_zf < 1e-10) {
        Rcpp::stop("Factor gene has zero variance in sender cells");
    }

    if (verbose) {
        Rcpp::Rcout << "Building weight matrices for " << n_radii << " radii..." << std::endl;
    }

    // Build weight matrices for all radii
    std::vector<arma::sp_mat> W_list;
    W_list.reserve(n_radii);

    for (int r = 0; r < n_radii; r++) {
        double r_inner = (r == 0) ? 0.0 : radii(r - 1);
        double r_outer = radii(r);

        // Build sparse weight matrix
        std::vector<arma::uword> row_idx_vec, col_idx_vec;
        std::vector<double> vals_vec;
        std::vector<double> row_sums(n_senders, 0.0);

        // First pass: count neighbors
        for (int i = 0; i < n_senders; i++) {
            int si = sender_idx(i);
            double xi = coords(si, 0);
            double yi = coords(si, 1);

            for (int j = 0; j < n_receivers; j++) {
                int rj = receiver_idx(j);
                if ((arma::uword)si == (arma::uword)rj) continue;

                double dx = xi - coords(rj, 0);
                double dy = yi - coords(rj, 1);
                double dist = std::sqrt(dx * dx + dy * dy);

                if (dist >= r_inner && dist < r_outer) {
                    row_sums[i] += 1.0;
                }
            }
        }

        // Second pass: build matrix
        for (int i = 0; i < n_senders; i++) {
            if (row_sums[i] < 1e-10) continue;

            int si = sender_idx(i);
            double xi = coords(si, 0);
            double yi = coords(si, 1);

            for (int j = 0; j < n_receivers; j++) {
                int rj = receiver_idx(j);
                if ((arma::uword)si == (arma::uword)rj) continue;

                double dx = xi - coords(rj, 0);
                double dy = yi - coords(rj, 1);
                double dist = std::sqrt(dx * dx + dy * dy);

                if (dist >= r_inner && dist < r_outer) {
                    row_idx_vec.push_back(i);
                    col_idx_vec.push_back(j);
                    vals_vec.push_back(1.0 / row_sums[i]);
                }
            }
        }

        if (row_idx_vec.empty()) {
            W_list.push_back(arma::sp_mat(n_senders, n_receivers));
        } else {
            arma::umat locs(2, row_idx_vec.size());
            for (size_t k = 0; k < row_idx_vec.size(); k++) {
                locs(0, k) = row_idx_vec[k];
                locs(1, k) = col_idx_vec[k];
            }
            arma::vec vals(vals_vec);
            W_list.push_back(arma::sp_mat(locs, vals, n_senders, n_receivers));
        }
    }

    // Pre-generate permutation indices
    if (verbose) {
        Rcpp::Rcout << "Generating " << n_perm << " permutations..." << std::endl;
    }

    arma::umat perm_indices(n_perm, n_senders);
    for (int b = 0; b < n_perm; b++) {
        perm_indices.row(b) = arma::randperm(n_senders).t();
    }

    // Result vectors
    Rcpp::NumericVector out_I_r1(n_genes);
    Rcpp::NumericVector out_p_value(n_genes);
    Rcpp::NumericVector out_z_score(n_genes);
    Rcpp::NumericVector out_delta_I(n_genes);

    if (verbose) {
        Rcpp::Rcout << "Testing " << n_genes << " genes at " << n_radii << " radii..." << std::endl;
    }

    // Process each gene
    for (int g = 0; g < n_genes; g++) {
        if (verbose && (g + 1) % 500 == 0) {
            Rcpp::Rcout << "  " << (g + 1) << " / " << n_genes << std::endl;
        }

        if (g % 100 == 0) {
            Rcpp::checkUserInterrupt();
        }

        // Extract gene expression in receivers
        arma::vec gene_full = expr_matrix.row(g).t();
        arma::vec z_g(n_receivers);
        for (arma::uword j = 0; j < receiver_idx.n_elem; j++) {
            z_g(j) = gene_full(receiver_idx(j));
        }
        z_g = perm_standardize(z_g);

        // Compute I_ND at each radius
        arma::vec I_curve(n_radii);
        for (int r = 0; r < n_radii; r++) {
            arma::vec lag_g = W_list[r] * z_g;
            double norm_lag = arma::norm(lag_g);

            if (norm_lag < 1e-10) {
                I_curve(r) = NA_REAL;
            } else {
                I_curve(r) = arma::dot(z_f, lag_g) / (norm_zf * norm_lag);
            }
        }

        // Compute delta_I (first radius - last radius, signed)
        double I_r1 = I_curve(0);
        double delta_I = NA_REAL;

        if (!R_IsNA(I_r1)) {
            double I_max = I_r1, I_min = I_r1;
            double I_last = I_r1;
            for (int r = 1; r < n_radii; r++) {
                if (!R_IsNA(I_curve(r))) {
                    if (I_curve(r) > I_max) I_max = I_curve(r);
                    if (I_curve(r) < I_min) I_min = I_curve(r);
                    I_last = I_curve(r);
                }
            }
            double range = I_max - I_min;
            int sign = (I_r1 >= I_last) ? 1 : -1;
            delta_I = sign * range;
        }

        // Permutation test at first radius
        if (R_IsNA(I_r1)) {
            out_I_r1[g] = NA_REAL;
            out_p_value[g] = NA_REAL;
            out_z_score[g] = NA_REAL;
            out_delta_I[g] = NA_REAL;
            continue;
        }

        arma::vec lag_r1 = W_list[0] * z_g;
        double norm_lag_r1 = arma::norm(lag_r1);

        int count_extreme = 0;
        double sum_perm = 0.0;
        double sum_perm_sq = 0.0;

        for (int b = 0; b < n_perm; b++) {
            arma::vec z_f_perm(n_senders);
            for (int i = 0; i < n_senders; i++) {
                z_f_perm(i) = z_f(perm_indices(b, i));
            }

            double I_perm = arma::dot(z_f_perm, lag_r1) / (norm_zf * norm_lag_r1);

            sum_perm += I_perm;
            sum_perm_sq += I_perm * I_perm;

            if (std::abs(I_perm) >= std::abs(I_r1)) {
                count_extreme++;
            }
        }

        double p_value = (1.0 + count_extreme) / (n_perm + 1.0);
        double null_mean = sum_perm / n_perm;
        double null_var = sum_perm_sq / n_perm - null_mean * null_mean;
        double null_sd = std::sqrt(std::max(null_var, 1e-10));
        double z_score = (I_r1 - null_mean) / null_sd;

        out_I_r1[g] = I_r1;
        out_p_value[g] = p_value;
        out_z_score[g] = z_score;
        out_delta_I[g] = delta_I;
    }

    if (verbose) {
        Rcpp::Rcout << "Done!" << std::endl;
    }

    return Rcpp::DataFrame::create(
        Rcpp::Named("gene") = gene_names,
        Rcpp::Named("I_ND_r1") = out_I_r1,
        Rcpp::Named("p_value") = out_p_value,
        Rcpp::Named("z_score") = out_z_score,
        Rcpp::Named("delta_I") = out_delta_I
    );
}

// =============================================================================
// VECTORIZED Multi-Factor Permutation Test
//
// Highly optimized version that uses matrix multiplication to compute
// ALL permuted statistics for ALL genes at once.
//
// Key optimizations:
// 1. Precompute spatial lags for ALL genes: LAG = W × Z_receiver (once)
// 2. Precompute permuted factor matrix: Z_PERM (n_senders × n_perm)
// 3. For each factor: I_perm = Z_PERM' × LAG (single BLAS matrix multiply!)
//
// Complexity: O(n_factors × n_senders × n_genes) instead of
//             O(n_factors × n_genes × n_perm × n_senders)
//
// @param expr_matrix Genes × cells expression matrix
// @param gene_names Gene name vector
// @param factor_indices Vector of factor gene indices (0-based)
// @param factor_names Factor gene names
// @param W Weight matrix (sparse, n_senders × n_receivers)
// @param sender_idx Sender cell indices (0-based)
// @param receiver_idx Receiver cell indices (0-based)
// @param n_perm Number of permutations
// @param seed Random seed
// @param verbose Print progress
//
// @return List of DataFrames, one per factor
// =============================================================================
// [[Rcpp::export]]
Rcpp::List cpp_batch_permutation_multi_factor_vectorized(
    const arma::mat& expr_matrix,
    const Rcpp::CharacterVector& gene_names,
    const arma::ivec& factor_indices,
    const Rcpp::CharacterVector& factor_names,
    const arma::sp_mat& W,
    const arma::uvec& sender_idx,
    const arma::uvec& receiver_idx,
    int n_perm = 999,
    int seed = -1,
    bool verbose = true
) {
    int n_genes = expr_matrix.n_rows;
    int n_factors = factor_indices.n_elem;
    int n_senders = sender_idx.n_elem;
    int n_receivers = receiver_idx.n_elem;

    if (seed >= 0) {
        arma::arma_rng::set_seed(seed);
    }

    // =========================================================================
    // STEP 1: Precompute spatial lags for ALL genes (expensive, but done ONCE)
    // LAG matrix: n_senders × n_genes
    // =========================================================================
    if (verbose) {
        Rcpp::Rcout << "Step 1/4: Precomputing spatial lags for " << n_genes
                    << " genes..." << std::endl;
    }

    // Extract receiver expression for all genes
    arma::mat Z_receiver(n_receivers, n_genes);
    for (int g = 0; g < n_genes; g++) {
        arma::vec gene_full = expr_matrix.row(g).t();
        for (int j = 0; j < n_receivers; j++) {
            Z_receiver(j, g) = gene_full(receiver_idx(j));
        }
        // Standardize column
        double mean_g = arma::mean(Z_receiver.col(g));
        double std_g = arma::stddev(Z_receiver.col(g), 1);
        if (std_g > 1e-10) {
            Z_receiver.col(g) = (Z_receiver.col(g) - mean_g) / std_g;
        } else {
            Z_receiver.col(g).zeros();
        }
    }

    // Compute spatial lags: LAG = W × Z_receiver (n_senders × n_genes)
    arma::mat LAG = W * Z_receiver;

    // Precompute norms of spatial lags for each gene
    arma::vec norm_lags(n_genes);
    for (int g = 0; g < n_genes; g++) {
        norm_lags(g) = arma::norm(LAG.col(g));
    }

    // =========================================================================
    // STEP 2: Generate permutation indices ONCE
    // =========================================================================
    if (verbose) {
        Rcpp::Rcout << "Step 2/4: Generating " << n_perm
                    << " permutation indices..." << std::endl;
    }

    arma::umat perm_indices(n_perm, n_senders);
    for (int b = 0; b < n_perm; b++) {
        perm_indices.row(b) = arma::randperm(n_senders).t();
    }

    // =========================================================================
    // STEP 3: Extract and standardize ALL factor expressions
    // =========================================================================
    if (verbose) {
        Rcpp::Rcout << "Step 3/4: Extracting " << n_factors
                    << " factor expressions..." << std::endl;
    }

    arma::mat Z_factors(n_senders, n_factors);
    arma::vec norm_factors(n_factors);

    for (int f = 0; f < n_factors; f++) {
        int fidx = factor_indices(f);
        arma::vec factor_full = expr_matrix.row(fidx).t();
        for (int i = 0; i < n_senders; i++) {
            Z_factors(i, f) = factor_full(sender_idx(i));
        }
        // Standardize
        double mean_f = arma::mean(Z_factors.col(f));
        double std_f = arma::stddev(Z_factors.col(f), 1);
        if (std_f > 1e-10) {
            Z_factors.col(f) = (Z_factors.col(f) - mean_f) / std_f;
        } else {
            Z_factors.col(f).zeros();
        }
        norm_factors(f) = arma::norm(Z_factors.col(f));
    }

    // =========================================================================
    // STEP 4: Process each factor using VECTORIZED permutation test
    // =========================================================================
    if (verbose) {
        Rcpp::Rcout << "Step 4/4: Running vectorized permutation tests..."
                    << std::endl;
    }

    Rcpp::List results(n_factors);

    for (int f = 0; f < n_factors; f++) {
        if (verbose) {
            Rcpp::Rcout << "  Factor " << (f + 1) << "/" << n_factors
                        << ": " << factor_names[f] << std::endl;
        }

        Rcpp::checkUserInterrupt();

        arma::vec z_f = Z_factors.col(f);
        double norm_f = norm_factors(f);

        if (norm_f < 1e-10) {
            // Factor has zero variance - skip
            Rcpp::NumericVector na_vec(n_genes, NA_REAL);
            results[f] = Rcpp::DataFrame::create(
                Rcpp::Named("gene") = gene_names,
                Rcpp::Named("I_obs") = na_vec,
                Rcpp::Named("p_value") = na_vec,
                Rcpp::Named("p_adj") = na_vec,
                Rcpp::Named("z_score") = na_vec
            );
            continue;
        }

        // -----------------------------------------------------------------
        // Compute OBSERVED statistics for ALL genes at once
        // I_obs = z_f' × LAG / (norm_f × norm_lags)
        // This is a single matrix-vector multiply!
        // -----------------------------------------------------------------
        arma::vec dot_obs = LAG.t() * z_f;  // n_genes × 1
        arma::vec I_obs = dot_obs / (norm_f * norm_lags);

        // Handle genes with zero spatial lag
        for (int g = 0; g < n_genes; g++) {
            if (norm_lags(g) < 1e-10) {
                I_obs(g) = NA_REAL;
            }
        }

        // -----------------------------------------------------------------
        // Build permuted factor matrix: Z_PERM (n_senders × n_perm)
        // -----------------------------------------------------------------
        arma::mat Z_f_perm(n_senders, n_perm);
        for (int b = 0; b < n_perm; b++) {
            for (int i = 0; i < n_senders; i++) {
                Z_f_perm(i, b) = z_f(perm_indices(b, i));
            }
        }

        // -----------------------------------------------------------------
        // Compute ALL permuted statistics with ONE matrix multiply!
        // I_perm_all = Z_PERM' × LAG gives (n_perm × n_genes) matrix
        // -----------------------------------------------------------------
        arma::mat dot_perm = Z_f_perm.t() * LAG;  // n_perm × n_genes

        // Normalize by norms
        // For standardized data, norm of permuted vector equals original
        arma::mat I_perm_all(n_perm, n_genes);
        for (int g = 0; g < n_genes; g++) {
            if (norm_lags(g) > 1e-10) {
                I_perm_all.col(g) = dot_perm.col(g) / (norm_f * norm_lags(g));
            } else {
                I_perm_all.col(g).fill(NA_REAL);
            }
        }

        // -----------------------------------------------------------------
        // Compute p-values and z-scores for all genes
        // -----------------------------------------------------------------
        Rcpp::NumericVector out_I_obs(n_genes);
        Rcpp::NumericVector out_p_value(n_genes);
        Rcpp::NumericVector out_z_score(n_genes);

        for (int g = 0; g < n_genes; g++) {
            out_I_obs[g] = I_obs(g);

            if (R_IsNA(I_obs(g))) {
                out_p_value[g] = NA_REAL;
                out_z_score[g] = NA_REAL;
                continue;
            }

            // Count extreme values (vectorized)
            arma::vec perm_col = I_perm_all.col(g);
            int count_extreme = arma::sum(arma::abs(perm_col) >= std::abs(I_obs(g)));

            // P-value
            out_p_value[g] = (1.0 + count_extreme) / (n_perm + 1.0);

            // Z-score
            double null_mean = arma::mean(perm_col);
            double null_sd = arma::stddev(perm_col);
            out_z_score[g] = (I_obs(g) - null_mean) / (null_sd + 1e-10);
        }

        // Create DataFrame for this factor
        results[f] = Rcpp::DataFrame::create(
            Rcpp::Named("gene") = gene_names,
            Rcpp::Named("I_obs") = out_I_obs,
            Rcpp::Named("p_value") = out_p_value,
            Rcpp::Named("z_score") = out_z_score
        );
    }

    results.attr("names") = factor_names;

    if (verbose) {
        Rcpp::Rcout << "Done!" << std::endl;
    }

    return results;
}
