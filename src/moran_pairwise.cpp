// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>

using namespace Rcpp;
using namespace arma;

// Constants
#define VISIUM 0
#define OLD 1
#define VISIUM_DISTANCE 100
#define OLD_DISTANCE 200

// Distance decay function
inline double distance_decay(double d) {
    return std::exp(-(d / 100.0) * (d / 100.0) / 2.0);
}

// Create distance lookup table
// [[Rcpp::export]]
NumericMatrix cpp_create_distance(int max_shift, int mode) {
    const double visum_shift = 0.5 * std::sqrt(3.0);
    const double old_shift = 0.5;

    NumericMatrix distance(max_shift, 2 * max_shift);

    for (int i = 0; i < max_shift; i++) {
        for (int j = 0; j < 2 * max_shift; j++) {
            double x, y;

            if (mode == VISIUM) {
                x = 0.5 * j * VISIUM_DISTANCE;
                y = i * visum_shift * VISIUM_DISTANCE;
            } else if (mode == OLD) {
                x = 0.5 * j * OLD_DISTANCE;
                y = i * old_shift * OLD_DISTANCE;
            } else {
                stop("Unknown platform mode");
            }

            double d = std::sqrt(x * x + y * y);
            distance(i, j) = distance_decay(d);
        }
    }

    return distance;
}

// Z-normalize rows of a matrix (genes are rows)
// [[Rcpp::export]]
arma::mat cpp_z_normalize(arma::mat data) {
    int nrow = data.n_rows;

    for (int i = 0; i < nrow; i++) {
        // Compute mean
        double mean = arma::mean(data.row(i));

        // Compute std
        double var = arma::var(data.row(i), 1);  // 1 = normalize by N (population)
        double std_dev = std::sqrt(var);

        // Normalize
        if (std_dev > 0) {
            data.row(i) = (data.row(i) - mean) / std_dev;
        }
    }

    return data;
}

// Create weight matrix W' = W + diag(W)
// [[Rcpp::export]]
List cpp_create_weight_matrix(
    IntegerVector spot_row,
    IntegerVector spot_col,
    NumericMatrix distance,
    int max_shift,
    bool flag_samespot) {

    int n_spots = spot_row.size();
    arma::mat W(n_spots, n_spots, arma::fill::zeros);
    double weight_sum = 0.0;

    for (int i = 0; i < n_spots; i++) {
        for (int j = i; j < n_spots; j++) {
            int row_diff = spot_row[j] - spot_row[i];
            int col_diff = spot_col[j] - spot_col[i];

            // Only consider upper triangle in spatial terms
            if (row_diff < 0) continue;
            if (row_diff == 0 && col_diff < 0) continue;

            int row_offset = row_diff;
            int col_offset = std::abs(col_diff);

            // Check if within radius
            if (row_offset >= max_shift || col_offset >= 2 * max_shift) continue;

            double w = distance(row_offset, col_offset);

            // Handle same spot case
            if (i == j && !flag_samespot) {
                w = 0.0;
            }

            if (w > 0.0) {
                if (i == j) {
                    // Diagonal: W'[i,i] = 2 * W[i,i]
                    W(i, i) = 2.0 * w;
                    weight_sum += 2.0 * w;
                } else {
                    // Off-diagonal: symmetric
                    W(i, j) = w;
                    W(j, i) = w;
                    weight_sum += 2.0 * w;
                }
            }
        }
    }

    return List::create(
        Named("W") = W,
        Named("weight_sum") = weight_sum
    );
}

// Main pairwise Moran's I computation using BLAS
// [[Rcpp::export]]
arma::mat cpp_pairwise_moran(
    arma::mat data,           // genes x spots (already z-normalized)
    arma::mat W,              // weight matrix (spots x spots)
    double weight_sum,
    bool paired_genes = true,
    bool all_genes = true) {

    if (weight_sum == 0.0) {
        stop("Weight sum is zero. Check radius and data.");
    }

    // Step 1: temp = data @ W  (n x ncol) @ (ncol x ncol) = (n x ncol)
    arma::mat temp = data * W;

    // Step 2: result = temp @ data^T  (n x ncol) @ (ncol x n) = (n x n)
    arma::mat result = temp * data.t();

    // Normalize by weight_sum
    result /= weight_sum;

    if (paired_genes) {
        if (all_genes) {
            // Return full matrix (will extract lower triangular in R)
            return result;
        } else {
            // Only pairs with first gene
            return result.col(0);
        }
    } else {
        // Single gene mode (diagonal only)
        return arma::diagvec(result);
    }
}

// Sparse matrix version - converts to dense for computation
// [[Rcpp::export]]
arma::mat cpp_pairwise_moran_sparse(
    arma::sp_mat data_sparse,  // sparse genes x spots
    arma::mat W,               // weight matrix (spots x spots)
    double weight_sum,
    bool paired_genes = true,
    bool all_genes = true) {

    // Convert sparse to dense (BLAS is more efficient for dense operations)
    arma::mat data = arma::mat(data_sparse);

    return cpp_pairwise_moran(data, W, weight_sum, paired_genes, all_genes);
}

// All-in-one function for convenience
// [[Rcpp::export]]
arma::mat cpp_compute_moran_full(
    arma::mat data,           // genes x spots (raw, will be normalized)
    IntegerVector spot_row,
    IntegerVector spot_col,
    int max_radius = 5,
    int platform = 0,         // 0 = VISIUM, 1 = OLD
    bool flag_samespot = true,
    bool paired_genes = true,
    bool all_genes = true,
    bool verbose = true) {

    int n_genes = data.n_rows;
    int n_spots = data.n_cols;

    if (verbose) {
        Rcout << "Pairwise Moran's I" << std::endl;
        Rcout << n_spots << " spots, " << n_genes << " genes" << std::endl;
    }

    // Step 1: Z-normalize
    if (verbose) Rcout << "Z-normalizing data..." << std::endl;
    data = cpp_z_normalize(data);

    // Step 2: Create distance lookup
    NumericMatrix distance = cpp_create_distance(max_radius, platform);
    if (!flag_samespot) {
        distance(0, 0) = 0.0;
    }

    // Step 3: Create weight matrix
    if (verbose) Rcout << "Creating weight matrix..." << std::endl;
    List W_result = cpp_create_weight_matrix(
        spot_row, spot_col, distance, max_radius, flag_samespot);
    arma::mat W = as<arma::mat>(W_result["W"]);
    double weight_sum = W_result["weight_sum"];

    if (verbose) Rcout << "Weight sum: " << weight_sum << std::endl;

    // Step 4: Compute pairwise Moran's I
    if (verbose) Rcout << "Computing pairwise Moran's I..." << std::endl;
    arma::mat result = cpp_pairwise_moran(data, W, weight_sum, paired_genes, all_genes);

    if (verbose) Rcout << "Done." << std::endl;

    return result;
}

// Sparse matrix version of all-in-one function
// [[Rcpp::export]]
arma::mat cpp_compute_moran_full_sparse(
    arma::sp_mat data_sparse,
    IntegerVector spot_row,
    IntegerVector spot_col,
    int max_radius = 5,
    int platform = 0,
    bool flag_samespot = true,
    bool paired_genes = true,
    bool all_genes = true,
    bool verbose = true) {

    // Convert to dense and call dense version
    arma::mat data = arma::mat(data_sparse);
    return cpp_compute_moran_full(
        data, spot_row, spot_col, max_radius, platform,
        flag_samespot, paired_genes, all_genes, verbose);
}

// =============================================================================
// SPARSE WEIGHT MATRIX IMPLEMENTATION
// =============================================================================

// Create SPARSE weight matrix - memory efficient for large datasets
// [[Rcpp::export]]
List cpp_create_weight_matrix_sparse(
    IntegerVector spot_row,
    IntegerVector spot_col,
    NumericMatrix distance,
    int max_shift,
    bool flag_samespot) {

    int n_spots = spot_row.size();

    // Pre-count entries to estimate capacity
    std::vector<arma::uword> row_indices;
    std::vector<arma::uword> col_indices;
    std::vector<double> values;

    // Reserve space (estimate ~20 neighbors per spot on average)
    size_t est_nnz = n_spots * 20;
    row_indices.reserve(est_nnz);
    col_indices.reserve(est_nnz);
    values.reserve(est_nnz);

    double weight_sum = 0.0;

    for (int i = 0; i < n_spots; i++) {
        for (int j = i; j < n_spots; j++) {
            int row_diff = spot_row[j] - spot_row[i];
            int col_diff = spot_col[j] - spot_col[i];

            // Only consider upper triangle in spatial terms
            if (row_diff < 0) continue;
            if (row_diff == 0 && col_diff < 0) continue;

            int row_offset = row_diff;
            int col_offset = std::abs(col_diff);

            // Check if within radius
            if (row_offset >= max_shift || col_offset >= 2 * max_shift) continue;

            double w = distance(row_offset, col_offset);

            // Handle same spot case
            if (i == j && !flag_samespot) {
                w = 0.0;
            }

            if (w > 0.0) {
                if (i == j) {
                    // Diagonal: W'[i,i] = 2 * W[i,i]
                    row_indices.push_back(i);
                    col_indices.push_back(i);
                    values.push_back(2.0 * w);
                    weight_sum += 2.0 * w;
                } else {
                    // Off-diagonal: symmetric (store both)
                    row_indices.push_back(i);
                    col_indices.push_back(j);
                    values.push_back(w);

                    row_indices.push_back(j);
                    col_indices.push_back(i);
                    values.push_back(w);

                    weight_sum += 2.0 * w;
                }
            }
        }
    }

    // Create sparse matrix from triplets
    arma::umat locations(2, row_indices.size());
    arma::vec vals(values.size());

    for (size_t k = 0; k < row_indices.size(); k++) {
        locations(0, k) = row_indices[k];
        locations(1, k) = col_indices[k];
        vals(k) = values[k];
    }

    arma::sp_mat W(locations, vals, n_spots, n_spots);

    return List::create(
        Named("W") = W,
        Named("weight_sum") = weight_sum,
        Named("nnz") = W.n_nonzero,
        Named("sparsity") = 1.0 - (double)W.n_nonzero / ((double)n_spots * n_spots)
    );
}

// Pairwise Moran's I with SPARSE weight matrix
// [[Rcpp::export]]
arma::mat cpp_pairwise_moran_W_sparse(
    arma::mat data,           // genes x spots (already z-normalized)
    arma::sp_mat W,           // SPARSE weight matrix (spots x spots)
    double weight_sum,
    bool paired_genes = true,
    bool all_genes = true) {

    if (weight_sum == 0.0) {
        stop("Weight sum is zero. Check radius and data.");
    }

    // Step 1: temp = data @ W_sparse  (dense x sparse = dense)
    // Armadillo handles this efficiently
    arma::mat temp = data * W;

    // Step 2: result = temp @ data^T  (dense x dense)
    arma::mat result = temp * data.t();

    // Normalize by weight_sum
    result /= weight_sum;

    if (paired_genes) {
        if (all_genes) {
            return result;
        } else {
            return result.col(0);
        }
    } else {
        return arma::diagvec(result);
    }
}

// All-in-one with SPARSE weight matrix
// [[Rcpp::export]]
arma::mat cpp_compute_moran_full_W_sparse(
    arma::mat data,           // genes x spots (raw, will be normalized)
    IntegerVector spot_row,
    IntegerVector spot_col,
    int max_radius = 5,
    int platform = 0,
    bool flag_samespot = true,
    bool paired_genes = true,
    bool all_genes = true,
    bool verbose = true) {

    int n_genes = data.n_rows;
    int n_spots = data.n_cols;

    if (verbose) {
        Rcout << "Pairwise Moran's I (Sparse W)" << std::endl;
        Rcout << n_spots << " spots, " << n_genes << " genes" << std::endl;
    }

    // Step 1: Z-normalize
    if (verbose) Rcout << "Z-normalizing data..." << std::endl;
    data = cpp_z_normalize(data);

    // Step 2: Create distance lookup
    NumericMatrix distance = cpp_create_distance(max_radius, platform);
    if (!flag_samespot) {
        distance(0, 0) = 0.0;
    }

    // Step 3: Create SPARSE weight matrix
    if (verbose) Rcout << "Creating sparse weight matrix..." << std::endl;
    List W_result = cpp_create_weight_matrix_sparse(
        spot_row, spot_col, distance, max_radius, flag_samespot);
    arma::sp_mat W = as<arma::sp_mat>(W_result["W"]);
    double weight_sum = W_result["weight_sum"];
    int nnz = W_result["nnz"];
    double sparsity = W_result["sparsity"];

    if (verbose) {
        Rcout << "Weight sum: " << weight_sum << std::endl;
        Rcout << "W non-zeros: " << nnz << " (sparsity: " <<
            (sparsity * 100) << "%)" << std::endl;
        Rcout << "W memory: ~" << (nnz * 16 / 1024) << " KB (vs dense: " <<
            (n_spots * n_spots * 8 / 1024 / 1024) << " MB)" << std::endl;
    }

    // Step 4: Compute pairwise Moran's I with sparse W
    if (verbose) Rcout << "Computing pairwise Moran's I..." << std::endl;
    arma::mat result = cpp_pairwise_moran_W_sparse(data, W, weight_sum, paired_genes, all_genes);

    if (verbose) Rcout << "Done." << std::endl;

    return result;
}
