#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

//' Create Gaussian Weight Matrix for Single-Cell Data
//'
//' Matches Python implementation: sigma = radius / 3
//'
//' @param sender_coords Sender cell coordinates (n_s x 2)
//' @param receiver_coords Receiver cell coordinates (n_r x 2)
//' @param radius Outer radius (micrometers)
//' @param inner_radius Inner radius for ring (default: 0)
//' @param sigma Gaussian sigma (default: radius/3)
//' @param min_weight Minimum weight threshold (default: 1e-6)
//' @return Sparse weight matrix (n_s x n_r), row-normalized
//' @export
// [[Rcpp::export]]
arma::sp_mat create_gaussian_weights_cpp(
    const arma::mat& sender_coords,
    const arma::mat& receiver_coords,
    double radius,
    double inner_radius = 0.0,
    double sigma = -1.0,
    double min_weight = 1e-6
) {
    // Default sigma = radius / 3 (matching Python)
    if (sigma < 0) {
        sigma = radius / 3.0;
    }

    int n_senders = sender_coords.n_rows;
    int n_receivers = receiver_coords.n_rows;

    double radius_sq = radius * radius;
    double inner_radius_sq = inner_radius * inner_radius;
    double gaussian_factor = -1.0 / (2.0 * sigma * sigma);

    // Triplet storage
    std::vector<arma::uword> rows, cols;
    std::vector<double> values;

    // Compute weights
    for (int i = 0; i < n_senders; i++) {
        double xi = sender_coords(i, 0);
        double yi = sender_coords(i, 1);

        std::vector<int> neighbors;
        std::vector<double> weights;
        double row_sum = 0.0;

        for (int j = 0; j < n_receivers; j++) {
            double dx = xi - receiver_coords(j, 0);
            double dy = yi - receiver_coords(j, 1);
            double dist_sq = dx * dx + dy * dy;

            // Check distance bounds
            if (dist_sq <= radius_sq && dist_sq > inner_radius_sq) {
                double w = std::exp(dist_sq * gaussian_factor);

                // Threshold small weights for sparsity
                if (w > min_weight) {
                    neighbors.push_back(j);
                    weights.push_back(w);
                    row_sum += w;
                }
            }
        }

        // Row normalize and store
        if (row_sum > 0) {
            for (size_t k = 0; k < neighbors.size(); k++) {
                rows.push_back(i);
                cols.push_back(neighbors[k]);
                values.push_back(weights[k] / row_sum);
            }
        }
    }

    // Build sparse matrix
    if (rows.empty()) {
        return arma::sp_mat(n_senders, n_receivers);
    }

    arma::umat locations(2, rows.size());
    for (size_t k = 0; k < rows.size(); k++) {
        locations(0, k) = rows[k];
        locations(1, k) = cols[k];
    }

    return arma::sp_mat(locations, arma::vec(values), n_senders, n_receivers);
}

//' Create Gaussian Ring Weight Matrix
//'
//' Weight matrix for annular region between inner and outer radius.
//'
//' @param coords Cell coordinates (n x 2)
//' @param outer_radius Outer radius
//' @param inner_radius Inner radius
//' @param sigma Gaussian sigma (default: outer_radius/3)
//' @return Sparse weight matrix (n x n), row-normalized
//' @export
// [[Rcpp::export]]
arma::sp_mat create_gaussian_ring_weights_cpp(
    const arma::mat& coords,
    double outer_radius,
    double inner_radius,
    double sigma = -1.0
) {
    return create_gaussian_weights_cpp(
        coords, coords, outer_radius, inner_radius, sigma, 1e-6
    );
}


//' Create Sparse Weight Matrix from Pre-computed Neighbors
//'
//' Efficiently creates weight matrix from pre-computed neighbor indices
//' (e.g., from RANN k-NN search). Scales to millions of cells.
//'
//' @param nn_idx Neighbor index matrix (n x k) from RANN::nn2, 1-indexed
//' @param nn_dist Neighbor distance matrix (n x k) from RANN::nn2
//' @param sigma Gaussian sigma for weight computation
//' @param radius Maximum distance cutoff (neighbors beyond this are ignored)
//' @param row_normalize Whether to row-normalize weights (default: TRUE)
//' @param self_weight Weight for self-connections (default: 0, no self)
//' @return List with sparse weight matrix W and weight_sum
//' @export
// [[Rcpp::export]]
List create_weights_from_neighbors_cpp(
    const arma::imat& nn_idx,
    const arma::mat& nn_dist,
    double sigma,
    double radius = -1.0,
    bool row_normalize = true,
    double self_weight = 0.0
) {
    arma::uword n = nn_idx.n_rows;
    arma::uword k = nn_idx.n_cols;

    double gaussian_factor = -1.0 / (2.0 * sigma * sigma);

    // If no radius specified, use all neighbors
    if (radius < 0) {
        radius = std::numeric_limits<double>::infinity();
    }

    // Triplet storage for sparse matrix
    std::vector<arma::uword> rows, cols;
    std::vector<double> values;
    rows.reserve(n * k);
    cols.reserve(n * k);
    values.reserve(n * k);

    // Row sums for normalization
    arma::vec row_sums(n, arma::fill::zeros);

    // First pass: compute weights and row sums
    std::vector<std::vector<std::pair<arma::uword, double>>> row_entries(n);

    #pragma omp parallel for schedule(dynamic)
    for (arma::uword i = 0; i < n; i++) {
        std::vector<std::pair<arma::uword, double>>& entries = row_entries[i];
        double local_sum = 0.0;

        for (arma::uword j = 0; j < k; j++) {
            int idx = nn_idx(i, j);
            double dist = nn_dist(i, j);

            // Skip invalid neighbors (RANN returns 0 for no neighbor)
            if (idx <= 0) continue;

            // Convert to 0-indexed
            arma::uword neighbor = static_cast<arma::uword>(idx - 1);

            // Skip if beyond radius
            if (dist > radius) continue;

            // Skip self-connections if self_weight is 0
            if (neighbor == i && self_weight == 0.0) continue;

            // Compute Gaussian weight
            double w;
            if (neighbor == i) {
                w = self_weight;
            } else if (dist > 0) {
                w = std::exp(dist * dist * gaussian_factor);
            } else {
                continue;  // Skip zero-distance non-self
            }

            if (w > 1e-10) {
                entries.push_back(std::make_pair(neighbor, w));
                local_sum += w;
            }
        }

        row_sums(i) = local_sum;
    }

    // Second pass: normalize and collect triplets
    double total_weight = 0.0;
    for (arma::uword i = 0; i < n; i++) {
        double norm = row_normalize && row_sums(i) > 0 ? row_sums(i) : 1.0;

        for (const auto& entry : row_entries[i]) {
            rows.push_back(i);
            cols.push_back(entry.first);
            double w = entry.second / norm;
            values.push_back(w);
            total_weight += w;
        }
    }

    // Build sparse matrix
    arma::sp_mat W;
    if (!rows.empty()) {
        arma::umat locations(2, rows.size());
        for (size_t i = 0; i < rows.size(); i++) {
            locations(0, i) = rows[i];
            locations(1, i) = cols[i];
        }
        W = arma::sp_mat(locations, arma::vec(values), n, n);
    } else {
        W = arma::sp_mat(n, n);
    }

    return List::create(
        Named("W") = W,
        Named("weight_sum") = total_weight,
        Named("n_edges") = rows.size()
    );
}


//' Compute Pairwise Moran's I from Pre-computed Neighbors (Scalable)
//'
//' Computes pairwise Moran's I for large single-cell datasets using
//' pre-computed neighbors. Avoids O(n²) distance computation.
//'
//' @param data Gene expression matrix (genes x cells), will be z-normalized
//' @param nn_idx Neighbor index matrix (n_cells x k) from RANN, 1-indexed
//' @param nn_dist Neighbor distance matrix (n_cells x k) from RANN
//' @param sigma Gaussian sigma for weights
//' @param radius Maximum distance cutoff
//' @param verbose Print progress messages
//' @return List with moran matrix and weight_sum
//' @export
// [[Rcpp::export]]
List pairwise_moran_from_neighbors_cpp(
    arma::mat data,
    const arma::imat& nn_idx,
    const arma::mat& nn_dist,
    double sigma,
    double radius = -1.0,
    bool verbose = true
) {
    arma::uword n_genes = data.n_rows;
    arma::uword n_cells = data.n_cols;

    if (verbose) {
        Rcpp::Rcout << "Pairwise Moran's I from Neighbors\n";
        Rcpp::Rcout << n_cells << " cells, " << n_genes << " genes\n";
    }

    // Create weight matrix from neighbors
    if (verbose) Rcpp::Rcout << "Creating weight matrix from neighbors...\n";
    List W_result = create_weights_from_neighbors_cpp(
        nn_idx, nn_dist, sigma, radius, true, 0.0
    );
    arma::sp_mat W = as<arma::sp_mat>(W_result["W"]);
    double S0 = as<double>(W_result["weight_sum"]);
    arma::uword n_edges = as<arma::uword>(W_result["n_edges"]);

    if (verbose) {
        Rcpp::Rcout << "W non-zeros: " << n_edges << "\n";
        Rcpp::Rcout << "S0: " << S0 << "\n";
    }

    if (S0 < 1e-10) {
        Rcpp::stop("Weight matrix has no non-zero weights");
    }

    // Z-normalize data (population SD)
    if (verbose) Rcpp::Rcout << "Z-normalizing data...\n";
    for (arma::uword i = 0; i < n_genes; i++) {
        arma::rowvec row = data.row(i);
        double mean_val = arma::mean(row);
        double var_val = arma::var(row, 1);  // population variance (N)
        double sd_val = std::sqrt(var_val);

        if (sd_val > 1e-10) {
            data.row(i) = (row - mean_val) / sd_val;
        } else {
            data.row(i).zeros();
        }
    }

    // Compute Moran's I: M = Z * W * Z^T / S0
    if (verbose) Rcpp::Rcout << "Computing pairwise Moran's I...\n";

    // Compute spatial lag: lag = data * W^T
    arma::mat lag = data * W.t();

    // Compute Moran matrix: M = data * lag^T / S0
    arma::mat moran = (data * lag.t()) / S0;

    if (verbose) Rcpp::Rcout << "Done.\n";

    return List::create(
        Named("moran") = moran,
        Named("weight_sum") = S0,
        Named("n_edges") = n_edges
    );
}
