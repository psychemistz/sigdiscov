// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "nanoflann.hpp"
#include <vector>
#include <cmath>

using namespace Rcpp;

// Forward declaration of PointCloud from neighbor_search.cpp
struct PointCloudDirectional {
    const arma::mat& pts;

    PointCloudDirectional(const arma::mat& points) : pts(points) {}

    inline size_t kdtree_get_point_count() const { return pts.n_rows; }

    inline double kdtree_get_pt(const size_t idx, const size_t dim) const {
        return pts(idx, dim);
    }

    template <class BBOX>
    bool kdtree_get_bbox(BBOX&) const { return false; }
};

typedef nanoflann::KDTreeSingleIndexAdaptor<
    nanoflann::L2_Simple_Adaptor<double, PointCloudDirectional>,
    PointCloudDirectional,
    2,
    size_t
> KDTreeDirectional;


//' Create Directional Weight Matrix for Single-Cell Cell Type Pair Analysis
//'
//' Creates a sparse Gaussian weight matrix for sender->receiver cell type pairs.
//' Weights are computed from sender cells to receiver cells within radius.
//' Uses KD-tree for efficient neighbor search on large datasets.
//' Supports annular (ring) weights via inner_radius parameter.
//'
//' @param sender_coords Sender cell coordinates (n_sender x 2)
//' @param receiver_coords Receiver cell coordinates (n_receiver x 2)
//' @param radius Distance radius (outer radius) in same units as coordinates
//' @param inner_radius Inner radius for annular weights (default: 0, circular)
//' @param sigma Gaussian sigma (default: radius/3)
//' @param use_kdtree Use KD-tree for large datasets (default: TRUE)
//' @param verbose Print progress messages (default: FALSE)
//' @return List with W (sparse n_sender x n_receiver), weight_sum, n_edges
//' @export
// [[Rcpp::export]]
List create_directional_weights_sc_cpp(
    const arma::mat& sender_coords,
    const arma::mat& receiver_coords,
    double radius,
    double inner_radius = 0.0,
    double sigma = -1.0,
    bool use_kdtree = true,
    bool verbose = false
) {
    arma::uword n_sender = sender_coords.n_rows;
    arma::uword n_receiver = receiver_coords.n_rows;

    if (sigma < 0) {
        sigma = radius / 3.0;
    }

    double radius_sq = radius * radius;
    double inner_radius_sq = inner_radius * inner_radius;
    double gaussian_factor = -1.0 / (2.0 * sigma * sigma);

    if (verbose && inner_radius > 0) {
        Rcpp::Rcout << "Annular weights: inner=" << inner_radius
                    << ", outer=" << radius << "\n";
    }

    // Triplet storage
    std::vector<arma::uword> rows, cols;
    std::vector<double> values;

    double total_weight = 0.0;
    arma::uword n_edges = 0;

    if (use_kdtree && n_receiver > 1000) {
        // Use KD-tree for large receiver sets
        PointCloudDirectional cloud(receiver_coords);
        KDTreeDirectional index(2, cloud, nanoflann::KDTreeSingleIndexAdaptorParams(10));
        index.buildIndex();

        nanoflann::SearchParameters params;
        params.sorted = false;

        for (arma::uword i = 0; i < n_sender; i++) {
            double query_pt[2] = {sender_coords(i, 0), sender_coords(i, 1)};

            std::vector<nanoflann::ResultItem<size_t, double>> matches;
            index.radiusSearch(query_pt, radius_sq, matches, params);

            // First pass: compute row sum for normalization
            double row_sum = 0.0;
            std::vector<std::pair<arma::uword, double>> cell_weights;
            cell_weights.reserve(matches.size());

            for (const auto& match : matches) {
                double dist_sq = match.second;  // dist_sq is already in match
                // Check inner radius for annular weights
                if (dist_sq < inner_radius_sq) continue;

                arma::uword j = static_cast<arma::uword>(match.first);
                double w = std::exp(dist_sq * gaussian_factor);

                if (w > 1e-10) {
                    cell_weights.push_back(std::make_pair(j, w));
                    row_sum += w;
                }
            }

            // Second pass: normalize and store
            if (row_sum > 0) {
                for (const auto& cw : cell_weights) {
                    rows.push_back(i);
                    cols.push_back(cw.first);
                    double w_norm = cw.second / row_sum;
                    values.push_back(w_norm);
                    total_weight += w_norm;
                    n_edges++;
                }
            }
        }
    } else {
        // Brute force for small datasets
        for (arma::uword i = 0; i < n_sender; i++) {
            double xi = sender_coords(i, 0);
            double yi = sender_coords(i, 1);

            double row_sum = 0.0;
            std::vector<std::pair<arma::uword, double>> cell_weights;

            for (arma::uword j = 0; j < n_receiver; j++) {
                double dx = xi - receiver_coords(j, 0);
                double dy = yi - receiver_coords(j, 1);
                double dist_sq = dx * dx + dy * dy;

                if (dist_sq <= radius_sq && dist_sq >= inner_radius_sq) {
                    double w = std::exp(dist_sq * gaussian_factor);
                    if (w > 1e-10) {
                        cell_weights.push_back(std::make_pair(j, w));
                        row_sum += w;
                    }
                }
            }

            if (row_sum > 0) {
                for (const auto& cw : cell_weights) {
                    rows.push_back(i);
                    cols.push_back(cw.first);
                    double w_norm = cw.second / row_sum;
                    values.push_back(w_norm);
                    total_weight += w_norm;
                    n_edges++;
                }
            }
        }
    }

    // Build sparse matrix
    arma::sp_mat W;
    if (!rows.empty()) {
        arma::umat locations(2, rows.size());
        for (size_t k = 0; k < rows.size(); k++) {
            locations(0, k) = rows[k];
            locations(1, k) = cols[k];
        }
        W = arma::sp_mat(locations, arma::vec(values), n_sender, n_receiver);
    } else {
        W = arma::sp_mat(n_sender, n_receiver);
    }

    return List::create(
        Named("W") = W,
        Named("weight_sum") = total_weight,
        Named("n_edges") = n_edges,
        Named("n_sender") = n_sender,
        Named("n_receiver") = n_receiver
    );
}


//' Compute Directional Pairwise Moran's I for Cell Type Pairs
//'
//' Computes gene x gene Moran's I matrix for sender->receiver cell type pairs.
//' Formula: I[i,j] = Z_sender[i,:] * W * Z_receiver[j,:]^T / S0
//'
//' @param sender_data Gene expression for sender cells (genes x n_sender), pre-normalized
//' @param receiver_data Gene expression for receiver cells (genes x n_receiver), pre-normalized
//' @param W Sparse directional weight matrix (n_sender x n_receiver), row-normalized
//' @param verbose Print progress messages (default: FALSE)
//' @return List with moran matrix (genes x genes) and weight_sum
//' @export
// [[Rcpp::export]]
List pairwise_moran_directional_cpp(
    const arma::mat& sender_data,
    const arma::mat& receiver_data,
    const arma::sp_mat& W,
    bool verbose = false
) {
    arma::uword n_genes = sender_data.n_rows;
    arma::uword n_sender = sender_data.n_cols;
    arma::uword n_receiver = receiver_data.n_cols;

    if (verbose) {
        Rcpp::Rcout << "Directional Pairwise Moran's I\n";
        Rcpp::Rcout << n_genes << " genes, " << n_sender << " senders, "
                    << n_receiver << " receivers\n";
    }

    // Validate dimensions
    if (W.n_rows != n_sender || W.n_cols != n_receiver) {
        Rcpp::stop("Weight matrix dimensions (%d x %d) don't match data (%d senders, %d receivers)",
                   W.n_rows, W.n_cols, n_sender, n_receiver);
    }

    if (sender_data.n_rows != receiver_data.n_rows) {
        Rcpp::stop("Sender and receiver data must have same number of genes");
    }

    // Compute S0 (sum of weights)
    double S0 = arma::accu(W);

    if (S0 < 1e-10) {
        Rcpp::stop("Weight matrix has no non-zero weights");
    }

    if (verbose) {
        Rcpp::Rcout << "S0: " << S0 << ", W nnz: " << W.n_nonzero << "\n";
    }

    // Compute spatial lag: lag = sender_data * W (genes x n_receiver)
    // This is the weighted average of receiver expression around each sender
    arma::mat lag = sender_data * W;  // (genes x n_sender) * (n_sender x n_receiver) = genes x n_receiver

    // Compute Moran's I: M = lag * receiver_data^T / S0
    // Each M[i,j] = sum_k(lag[i,k] * receiver_data[j,k]) / S0
    arma::mat moran = (lag * receiver_data.t()) / S0;

    if (verbose) {
        Rcpp::Rcout << "Done. Moran matrix: " << moran.n_rows << " x " << moran.n_cols << "\n";
    }

    return List::create(
        Named("moran") = moran,
        Named("weight_sum") = S0
    );
}


//' Compute Directional Pairwise Moran's I with Streaming (Memory Efficient)
//'
//' Computes directional Moran's I without storing the full weight matrix.
//' Uses KD-tree for efficient neighbor search. Best for large cell counts.
//'
//' @param sender_data Gene expression for sender cells (genes x n_sender)
//' @param receiver_data Gene expression for receiver cells (genes x n_receiver)
//' @param sender_coords Sender cell coordinates (n_sender x 2)
//' @param receiver_coords Receiver cell coordinates (n_receiver x 2)
//' @param radius Distance radius
//' @param sigma Gaussian sigma (default: radius/3)
//' @param normalize_data Z-normalize data before computation (default: TRUE)
//' @param verbose Print progress messages (default: TRUE)
//' @return List with moran matrix, weight_sum, n_edges
//' @export
// [[Rcpp::export]]
List pairwise_moran_directional_streaming_cpp(
    arma::mat sender_data,
    arma::mat receiver_data,
    const arma::mat& sender_coords,
    const arma::mat& receiver_coords,
    double radius,
    double sigma = -1.0,
    bool normalize_data = true,
    bool verbose = true
) {
    arma::uword n_genes = sender_data.n_rows;
    arma::uword n_sender = sender_data.n_cols;
    arma::uword n_receiver = receiver_data.n_cols;

    if (sigma < 0) {
        sigma = radius / 3.0;
    }

    if (verbose) {
        Rcpp::Rcout << "Directional Pairwise Moran's I (Streaming)\n";
        Rcpp::Rcout << n_genes << " genes, " << n_sender << " senders, "
                    << n_receiver << " receivers\n";
        Rcpp::Rcout << "radius=" << radius << ", sigma=" << sigma << "\n";
    }

    // Z-normalize data if requested
    if (normalize_data) {
        if (verbose) Rcpp::Rcout << "Z-normalizing data...\n";

        // Normalize sender data
        for (arma::uword i = 0; i < n_genes; i++) {
            arma::rowvec row = sender_data.row(i);
            double mean_val = arma::mean(row);
            double var_val = arma::var(row, 1);  // Population variance
            double sd_val = std::sqrt(var_val);

            if (sd_val > 1e-10) {
                sender_data.row(i) = (row - mean_val) / sd_val;
            } else {
                sender_data.row(i).zeros();
            }
        }

        // Normalize receiver data
        for (arma::uword i = 0; i < n_genes; i++) {
            arma::rowvec row = receiver_data.row(i);
            double mean_val = arma::mean(row);
            double var_val = arma::var(row, 1);
            double sd_val = std::sqrt(var_val);

            if (sd_val > 1e-10) {
                receiver_data.row(i) = (row - mean_val) / sd_val;
            } else {
                receiver_data.row(i).zeros();
            }
        }
    }

    // Build KD-tree on receiver coordinates
    if (verbose) Rcpp::Rcout << "Building KD-tree on receiver cells...\n";

    PointCloudDirectional cloud(receiver_coords);
    KDTreeDirectional index(2, cloud, nanoflann::KDTreeSingleIndexAdaptorParams(10));
    index.buildIndex();

    double radius_sq = radius * radius;
    double gaussian_factor = -1.0 / (2.0 * sigma * sigma);

    nanoflann::SearchParameters params;
    params.sorted = false;

    // Compute spatial lag on-the-fly
    // lag[:, j] = sum over receivers k of (w_jk * receiver_data[:, k])
    // where j is sender index, k is receiver index
    if (verbose) Rcpp::Rcout << "Computing spatial lag...\n";

    arma::mat lag(n_genes, n_sender, arma::fill::zeros);
    double total_weight = 0.0;
    arma::uword total_edges = 0;

    for (arma::uword i = 0; i < n_sender; i++) {
        double query_pt[2] = {sender_coords(i, 0), sender_coords(i, 1)};

        // Find receiver neighbors
        std::vector<nanoflann::ResultItem<size_t, double>> matches;
        index.radiusSearch(query_pt, radius_sq, matches, params);

        // First pass: compute row sum for normalization
        double row_sum = 0.0;
        std::vector<std::pair<arma::uword, double>> neighbors;
        neighbors.reserve(matches.size());

        for (const auto& match : matches) {
            arma::uword k = static_cast<arma::uword>(match.first);
            double w = std::exp(match.second * gaussian_factor);

            if (w > 1e-10) {
                neighbors.push_back(std::make_pair(k, w));
                row_sum += w;
            }
        }

        // Second pass: compute weighted lag
        if (row_sum > 0) {
            arma::vec lag_i(n_genes, arma::fill::zeros);

            for (const auto& nb : neighbors) {
                arma::uword k = nb.first;
                double w_norm = nb.second / row_sum;
                lag_i += w_norm * receiver_data.col(k);
                total_weight += w_norm;
                total_edges++;
            }

            lag.col(i) = lag_i;
        }

        // Progress reporting
        if (verbose && (i + 1) % 5000 == 0) {
            Rcpp::Rcout << "  Processed " << (i + 1) << "/" << n_sender << " senders\n";
        }

        // Check for interrupt
        if ((i + 1) % 1000 == 0) {
            Rcpp::checkUserInterrupt();
        }
    }

    if (verbose) {
        Rcpp::Rcout << "Total edges: " << total_edges << "\n";
        Rcpp::Rcout << "Avg receivers/sender: " << static_cast<double>(total_edges) / n_sender << "\n";
        Rcpp::Rcout << "S0: " << total_weight << "\n";
    }

    if (total_weight < 1e-10) {
        Rcpp::stop("No valid neighbors found within radius");
    }

    // Compute Moran's I matrix
    // M = sender_data * lag^T / S0
    // where lag is genes x n_sender
    if (verbose) Rcpp::Rcout << "Computing pairwise Moran's I...\n";

    arma::mat moran = (sender_data * lag.t()) / total_weight;

    if (verbose) Rcpp::Rcout << "Done.\n";

    return List::create(
        Named("moran") = moran,
        Named("weight_sum") = total_weight,
        Named("n_edges") = total_edges
    );
}


//' Compute Delta I from I(r) Curves (Batch)
//'
//' Efficiently computes delta I for multiple radii curves.
//' Delta I = sign(I_short - I_long) * (I_max - I_min)
//'
//' @param I_curves 3D array of I values (n_pairs x n_genes x n_genes x n_radii) stored as
//'        a matrix where each row is a flattened curve for one gene pair
//' @return List with delta_I, I_max, argmax matrices
//' @export
// [[Rcpp::export]]
List compute_delta_i_batch_cpp(
    const arma::mat& I_curves  // (n_pairs * n_genes^2) x n_radii
) {
    arma::uword n_entries = I_curves.n_rows;
    arma::uword n_radii = I_curves.n_cols;

    arma::vec delta_I(n_entries);
    arma::vec I_max(n_entries);
    arma::vec I_min(n_entries);
    arma::uvec argmax(n_entries);
    arma::ivec sign_vec(n_entries);

    #pragma omp parallel for schedule(static)
    for (arma::uword i = 0; i < n_entries; i++) {
        arma::rowvec curve = I_curves.row(i);

        // Handle all-NA/zero case
        if (arma::accu(arma::abs(curve)) < 1e-10) {
            delta_I(i) = 0.0;
            I_max(i) = 0.0;
            I_min(i) = 0.0;
            argmax(i) = 0;
            sign_vec(i) = 0;
            continue;
        }

        double max_val = curve.max();
        double min_val = curve.min();
        arma::uword max_idx = curve.index_max();

        double I_short = curve(0);
        double I_long = curve(n_radii - 1);

        int sign = (I_short >= I_long) ? 1 : -1;

        delta_I(i) = sign * (max_val - min_val);
        I_max(i) = max_val;
        I_min(i) = min_val;
        argmax(i) = max_idx;
        sign_vec(i) = sign;
    }

    return List::create(
        Named("delta_I") = delta_I,
        Named("I_max") = I_max,
        Named("I_min") = I_min,
        Named("argmax") = argmax,
        Named("sign") = sign_vec
    );
}


//' Compute Cell Type Pair Analysis for Multiple Radii
//'
//' Main function for computing Moran's I across multiple radii for a single
//' cell type pair. Returns I(r) curves for delta I computation.
//'
//' @param sender_data Gene expression for sender cells (genes x n_sender)
//' @param receiver_data Gene expression for receiver cells (genes x n_receiver)
//' @param sender_coords Sender cell coordinates (n_sender x 2)
//' @param receiver_coords Receiver cell coordinates (n_receiver x 2)
//' @param radii Vector of distance radii to compute
//' @param sigma_factor Sigma = radius * sigma_factor (default: 1/3)
//' @param verbose Print progress messages (default: TRUE)
//' @return List with I_curves (n_radii x n_genes x n_genes), delta_I, I_max, argmax
//' @export
// [[Rcpp::export]]
List compute_celltype_pair_moran_cpp(
    arma::mat sender_data,
    arma::mat receiver_data,
    const arma::mat& sender_coords,
    const arma::mat& receiver_coords,
    const arma::vec& radii,
    double sigma_factor = 0.333333,
    bool verbose = true
) {
    arma::uword n_genes = sender_data.n_rows;
    arma::uword n_sender = sender_data.n_cols;
    arma::uword n_receiver = receiver_data.n_cols;
    arma::uword n_radii = radii.n_elem;

    if (verbose) {
        Rcpp::Rcout << "Cell Type Pair Moran's I Analysis\n";
        Rcpp::Rcout << n_genes << " genes, " << n_sender << " senders, "
                    << n_receiver << " receivers, " << n_radii << " radii\n";
    }

    // Z-normalize data
    if (verbose) Rcpp::Rcout << "Z-normalizing data...\n";

    for (arma::uword i = 0; i < n_genes; i++) {
        arma::rowvec row = sender_data.row(i);
        double mean_val = arma::mean(row);
        double var_val = arma::var(row, 1);
        double sd_val = std::sqrt(var_val);

        if (sd_val > 1e-10) {
            sender_data.row(i) = (row - mean_val) / sd_val;
        } else {
            sender_data.row(i).zeros();
        }
    }

    for (arma::uword i = 0; i < n_genes; i++) {
        arma::rowvec row = receiver_data.row(i);
        double mean_val = arma::mean(row);
        double var_val = arma::var(row, 1);
        double sd_val = std::sqrt(var_val);

        if (sd_val > 1e-10) {
            receiver_data.row(i) = (row - mean_val) / sd_val;
        } else {
            receiver_data.row(i).zeros();
        }
    }

    // Build KD-tree on receiver coordinates
    if (verbose) Rcpp::Rcout << "Building KD-tree...\n";

    PointCloudDirectional cloud(receiver_coords);
    KDTreeDirectional index(2, cloud, nanoflann::KDTreeSingleIndexAdaptorParams(10));
    index.buildIndex();

    nanoflann::SearchParameters params;
    params.sorted = false;

    // Storage for I(r) curves
    // Using a 3D representation: n_radii matrices of n_genes x n_genes
    std::vector<arma::mat> I_curves(n_radii);
    for (arma::uword r = 0; r < n_radii; r++) {
        I_curves[r] = arma::mat(n_genes, n_genes, arma::fill::zeros);
    }

    // Process each radius
    for (arma::uword r = 0; r < n_radii; r++) {
        double radius = radii(r);
        double sigma = radius * sigma_factor;
        double radius_sq = radius * radius;
        double gaussian_factor = -1.0 / (2.0 * sigma * sigma);

        if (verbose) {
            Rcpp::Rcout << "  Radius " << (r + 1) << "/" << n_radii
                        << " (" << radius << ")...\n";
        }

        // Compute spatial lag
        arma::mat lag(n_genes, n_sender, arma::fill::zeros);
        double total_weight = 0.0;

        for (arma::uword i = 0; i < n_sender; i++) {
            double query_pt[2] = {sender_coords(i, 0), sender_coords(i, 1)};

            std::vector<nanoflann::ResultItem<size_t, double>> matches;
            index.radiusSearch(query_pt, radius_sq, matches, params);

            double row_sum = 0.0;
            std::vector<std::pair<arma::uword, double>> neighbors;
            neighbors.reserve(matches.size());

            for (const auto& match : matches) {
                arma::uword k = static_cast<arma::uword>(match.first);
                double w = std::exp(match.second * gaussian_factor);

                if (w > 1e-10) {
                    neighbors.push_back(std::make_pair(k, w));
                    row_sum += w;
                }
            }

            if (row_sum > 0) {
                arma::vec lag_i(n_genes, arma::fill::zeros);

                for (const auto& nb : neighbors) {
                    double w_norm = nb.second / row_sum;
                    lag_i += w_norm * receiver_data.col(nb.first);
                    total_weight += w_norm;
                }

                lag.col(i) = lag_i;
            }
        }

        // Compute Moran's I matrix for this radius
        if (total_weight > 1e-10) {
            I_curves[r] = (sender_data * lag.t()) / total_weight;
        }

        Rcpp::checkUserInterrupt();
    }

    // Compute delta I from curves
    if (verbose) Rcpp::Rcout << "Computing delta I...\n";

    arma::mat delta_I(n_genes, n_genes, arma::fill::zeros);
    arma::mat I_max(n_genes, n_genes, arma::fill::zeros);
    arma::umat argmax(n_genes, n_genes, arma::fill::zeros);

    for (arma::uword i = 0; i < n_genes; i++) {
        for (arma::uword j = 0; j < n_genes; j++) {
            // Extract I(r) curve for gene pair (i, j)
            arma::vec curve(n_radii);
            for (arma::uword r = 0; r < n_radii; r++) {
                curve(r) = I_curves[r](i, j);
            }

            // Compute delta I
            double max_val = curve.max();
            double min_val = curve.min();
            arma::uword max_idx = curve.index_max();

            double I_short = curve(0);
            double I_long = curve(n_radii - 1);
            int sign = (I_short >= I_long) ? 1 : -1;

            delta_I(i, j) = sign * (max_val - min_val);
            I_max(i, j) = max_val;
            argmax(i, j) = max_idx;
        }
    }

    if (verbose) Rcpp::Rcout << "Done.\n";

    // Return as list of matrices
    List I_curves_list(n_radii);
    for (arma::uword r = 0; r < n_radii; r++) {
        I_curves_list[r] = I_curves[r];
    }

    return List::create(
        Named("I_curves") = I_curves_list,
        Named("delta_I") = delta_I,
        Named("I_max") = I_max,
        Named("argmax") = argmax,
        Named("radii") = radii,
        Named("n_sender") = n_sender,
        Named("n_receiver") = n_receiver
    );
}
