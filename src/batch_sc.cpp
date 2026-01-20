#include <RcppArmadillo.h>
#include <map>
#include <set>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// Forward declaration
arma::sp_mat create_gaussian_weights_cpp(
    const arma::mat& sender_coords,
    const arma::mat& receiver_coords,
    double radius,
    double inner_radius,
    double sigma,
    double min_weight
);

//' Compute I_ND Matrix for One Cell Type Pair at One Radius
//'
//' Matches Python MoransIMatrixComputer.compute_morans_matrix_batched()
//'
//' @param factor_expr_norm Normalized factor expression (n_sender_cells x n_factors)
//' @param factor_expr_raw Raw factor expression (for quantile threshold)
//' @param gene_expr_norm Normalized gene expression (n_receiver_cells x n_genes)
//' @param W Weight matrix (n_senders x n_receivers, row-normalized)
//' @param quantile_prob Quantile threshold (default: 0.25 = keep top 75 percent)
//' @param min_cells Minimum high-expression cells
//' @return Matrix of I_ND values (n_factors x n_genes)
//' @export
// [[Rcpp::export]]
arma::mat compute_ind_matrix_cpp(
    const arma::mat& factor_expr_norm,   // n_sender_cells x n_factors
    const arma::mat& factor_expr_raw,    // n_sender_cells x n_factors (RAW)
    const arma::mat& gene_expr_norm,     // n_receiver_cells x n_genes
    const arma::sp_mat& W,
    double quantile_prob = 0.25,
    int min_cells = 10
) {
    int n_sender_cells = factor_expr_norm.n_rows;
    int n_factors = factor_expr_norm.n_cols;
    int n_genes = gene_expr_norm.n_cols;

    // Output matrix
    arma::mat IND_matrix(n_factors, n_genes, arma::fill::zeros);

    // Step 1: Pre-compute spatial lags for ALL genes (done ONCE)
    // This is the expensive operation - O(nnz(W) x n_genes)
    arma::mat spatial_lags_all = W * gene_expr_norm;  // n_sender_cells x n_genes

    // Step 2: For each factor
    for (int f = 0; f < n_factors; f++) {
        // Get RAW factor expression for quantile threshold
        arma::vec factor_raw = factor_expr_raw.col(f);

        // Compute quantile threshold on RAW expression
        arma::vec sorted_expr = arma::sort(factor_raw);
        int quantile_idx = static_cast<int>(n_sender_cells * quantile_prob);
        if (quantile_idx >= n_sender_cells) quantile_idx = n_sender_cells - 1;
        double threshold = sorted_expr(quantile_idx);

        // Get high-expression indices (> threshold, matching Python)
        arma::uvec high_idx = arma::find(factor_raw > threshold);
        int n_high = high_idx.n_elem;

        if (n_high < min_cells) {
            IND_matrix.row(f).fill(NA_REAL);
            continue;
        }

        // Extract NORMALIZED factor for high-expression cells
        arma::vec factor_high = factor_expr_norm.col(f);
        factor_high = factor_high.elem(high_idx);

        // Normalize factor vector
        double factor_norm = arma::norm(factor_high);
        if (factor_norm < 1e-10) {
            IND_matrix.row(f).fill(NA_REAL);
            continue;
        }
        arma::vec factor_normalized = factor_high / factor_norm;

        // Extract spatial lags for high-expression senders
        arma::mat spatial_lags_high = spatial_lags_all.rows(high_idx);

        // Compute I_ND for all genes at once (vectorized)
        arma::rowvec correlations = factor_normalized.t() * spatial_lags_high;
        arma::rowvec spatial_norms = arma::sqrt(arma::sum(arma::square(spatial_lags_high), 0));

        for (int g = 0; g < n_genes; g++) {
            if (spatial_norms(g) > 1e-10) {
                IND_matrix(f, g) = correlations(g) / spatial_norms(g);
            } else {
                IND_matrix(f, g) = NA_REAL;
            }
        }
    }

    return IND_matrix;
}

//' Batch Compute All Cell Type Pairs at All Radii
//'
//' Matches Python run_matrix_analysis() output format.
//'
//' @param expr_matrix Expression matrix (genes x cells)
//' @param coords Cell coordinates (cells x 2)
//' @param cell_types Cell type vector
//' @param factor_indices Factor gene indices (0-based)
//' @param radii Distance radii vector
//' @param pairs DataFrame with sender, receiver columns
//' @param quantile_prob Quantile threshold
//' @param min_cells Minimum cells per group
//' @param min_connections Minimum neighbors in weight matrix
//' @param verbose Print progress
//' @return List of cubes (one per radius), each is (n_pairs x n_factors x n_genes)
//' @export
// [[Rcpp::export]]
List batch_compute_all_pairs_cpp(
    const arma::mat& expr_matrix,        // genes x cells
    const arma::mat& coords,             // cells x 2
    const std::vector<std::string>& cell_types,
    const arma::uvec& factor_indices,    // 0-based
    const arma::vec& radii,
    const DataFrame& pairs,
    double quantile_prob = 0.25,
    int min_cells = 10,
    int min_connections = 10,
    bool verbose = true
) {
    CharacterVector senders = pairs["sender"];
    CharacterVector receivers = pairs["receiver"];

    int n_cells = expr_matrix.n_cols;
    int n_genes = expr_matrix.n_rows;
    int n_factors = factor_indices.n_elem;
    int n_radii = radii.n_elem;
    int n_pairs = senders.size();

    // Gene-wise standardization (matching Python)
    arma::mat expr_norm = expr_matrix;
    for (int g = 0; g < n_genes; g++) {
        double mean_g = arma::mean(expr_matrix.row(g));
        double std_g = arma::stddev(expr_matrix.row(g));
        if (std_g > 1e-10) {
            expr_norm.row(g) = (expr_matrix.row(g) - mean_g) / std_g;
        } else {
            expr_norm.row(g).zeros();
        }
    }

    // Cache cell indices by type
    std::map<std::string, arma::uvec> cell_idx_by_type;
    std::set<std::string> unique_types(cell_types.begin(), cell_types.end());

    for (const auto& ct : unique_types) {
        std::vector<arma::uword> indices;
        for (int i = 0; i < n_cells; i++) {
            if (cell_types[i] == ct) {
                indices.push_back(i);
            }
        }
        cell_idx_by_type[ct] = arma::conv_to<arma::uvec>::from(indices);
    }

    // Output: list of cubes
    List result(n_radii);

    // Process each radius
    for (int r = 0; r < n_radii; r++) {
        double radius = radii(r);
        double inner_radius = (r == 0) ? 0.0 : radii(r - 1);

        if (verbose) {
            Rcpp::Rcout << "Processing radius " << radius << " um..." << std::endl;
        }

        // Cube for this radius
        arma::cube radius_cube(n_pairs, n_factors, n_genes, arma::fill::zeros);

        // Process each pair
        for (int p = 0; p < n_pairs; p++) {
            std::string sender_type = Rcpp::as<std::string>(senders[p]);
            std::string receiver_type = Rcpp::as<std::string>(receivers[p]);

            arma::uvec sender_type_idx = cell_idx_by_type[sender_type];
            arma::uvec receiver_idx = cell_idx_by_type[receiver_type];

            int n_sender_cells = sender_type_idx.n_elem;
            int n_receiver_cells = receiver_idx.n_elem;

            if (n_sender_cells < min_cells || n_receiver_cells < min_cells) {
                for (int f = 0; f < n_factors; f++) {
                    for (int g = 0; g < n_genes; g++) {
                        radius_cube(p, f, g) = NA_REAL;
                    }
                }
                continue;
            }

            // Extract coordinates
            arma::mat sender_coords = coords.rows(sender_type_idx);
            arma::mat receiver_coords = coords.rows(receiver_idx);

            // Build Gaussian weight matrix
            arma::sp_mat W = create_gaussian_weights_cpp(
                sender_coords, receiver_coords, radius, inner_radius, -1.0, 1e-6
            );

            if (W.n_nonzero < static_cast<size_t>(min_connections)) {
                for (int f = 0; f < n_factors; f++) {
                    for (int g = 0; g < n_genes; g++) {
                        radius_cube(p, f, g) = NA_REAL;
                    }
                }
                continue;
            }

            // Extract expression for factors
            arma::mat factor_expr_norm_local(n_sender_cells, n_factors);
            arma::mat factor_expr_raw(n_sender_cells, n_factors);
            for (int f = 0; f < n_factors; f++) {
                int g_idx = factor_indices(f);
                arma::rowvec row_norm = expr_norm.row(g_idx);
                arma::rowvec row_raw = expr_matrix.row(g_idx);
                for (int k = 0; k < n_sender_cells; k++) {
                    factor_expr_norm_local(k, f) = row_norm(sender_type_idx(k));
                    factor_expr_raw(k, f) = row_raw(sender_type_idx(k));
                }
            }

            // Extract expression for all genes
            arma::mat gene_expr_norm_local(n_receiver_cells, n_genes);
            for (int g = 0; g < n_genes; g++) {
                arma::rowvec row_g = expr_norm.row(g);
                for (int k = 0; k < n_receiver_cells; k++) {
                    gene_expr_norm_local(k, g) = row_g(receiver_idx(k));
                }
            }

            // Compute I_ND matrix
            arma::mat IND_pair = compute_ind_matrix_cpp(
                factor_expr_norm_local, factor_expr_raw, gene_expr_norm_local, W,
                quantile_prob, min_cells
            );

            // Store in cube
            for (int f = 0; f < n_factors; f++) {
                for (int g = 0; g < n_genes; g++) {
                    radius_cube(p, f, g) = IND_pair(f, g);
                }
            }
        }

        result[r] = radius_cube;
    }

    return result;
}
