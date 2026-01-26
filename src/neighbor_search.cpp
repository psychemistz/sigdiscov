#include <RcppArmadillo.h>
#include "nanoflann.hpp"
#ifdef _OPENMP
#include <omp.h>
#endif
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// Adaptor for nanoflann to work with arma::mat
struct PointCloud {
    const arma::mat& pts;

    PointCloud(const arma::mat& points) : pts(points) {}

    // Must return the number of data points
    inline size_t kdtree_get_point_count() const { return pts.n_rows; }

    // Returns the dim'th component of the idx'th point
    inline double kdtree_get_pt(const size_t idx, const size_t dim) const {
        return pts(idx, dim);
    }

    // Optional bounding-box computation (not used)
    template <class BBOX>
    bool kdtree_get_bbox(BBOX&) const { return false; }
};

// KD-tree type definition with explicit size_t index type
typedef nanoflann::KDTreeSingleIndexAdaptor<
    nanoflann::L2_Simple_Adaptor<double, PointCloud>,
    PointCloud,
    2,      // 2D points
    size_t  // Index type
> KDTree;


//' Find Neighbors Within Radius Using KD-Tree (Native C++)
//'
//' Uses nanoflann for O(n log n) neighbor search. No R dependencies required.
//'
//' @param coords Cell coordinates (n x 2 matrix)
//' @param radius Search radius
//' @param max_neighbors Maximum neighbors to return per point (default: 200)
//' @return List with nn_idx (1-indexed) and nn_dist matrices
//' @export
// [[Rcpp::export]]
List find_neighbors_radius_cpp(
    const arma::mat& coords,
    double radius,
    int max_neighbors = 200
) {
    size_t n = coords.n_rows;

    // Build KD-tree
    PointCloud cloud(coords);
    KDTree index(2, cloud, nanoflann::KDTreeSingleIndexAdaptorParams(10));
    index.buildIndex();

    // Output matrices (1-indexed for R compatibility)
    arma::imat nn_idx(n, max_neighbors, arma::fill::zeros);
    arma::mat nn_dist(n, max_neighbors, arma::fill::zeros);

    // Search parameters
    nanoflann::SearchParameters params;
    params.sorted = true;

    double radius_sq = radius * radius;

    // Find neighbors for each point
    #pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < n; i++) {
        double query_pt[2] = {coords(i, 0), coords(i, 1)};

        // Radius search returns pairs of (index, distance_squared)
        std::vector<nanoflann::ResultItem<size_t, double>> matches;
        size_t n_matches = index.radiusSearch(query_pt, radius_sq, matches, params);

        // Store results (up to max_neighbors)
        size_t n_store = std::min(n_matches, static_cast<size_t>(max_neighbors));
        for (size_t j = 0; j < n_store; j++) {
            nn_idx(i, j) = static_cast<int>(matches[j].first + 1);  // 1-indexed
            nn_dist(i, j) = std::sqrt(matches[j].second);  // Convert to distance
        }
    }

    return List::create(
        Named("nn.idx") = nn_idx,
        Named("nn.dists") = nn_dist
    );
}


//' Compute Pairwise Moran's I for Large Datasets (Native C++)
//'
//' Complete native C++ implementation using nanoflann for neighbor search.
//' No R package dependencies (RANN not required).
//'
//' @param data Gene expression matrix (genes x cells)
//' @param coords Cell coordinates (n x 2 matrix)
//' @param radius Radius for neighbor search
//' @param sigma Gaussian sigma (default: radius/3)
//' @param max_neighbors Maximum neighbors per cell (default: 200)
//' @param verbose Print progress messages (default: TRUE)
//' @return List with moran matrix, weight_sum, and n_edges
//' @export
// [[Rcpp::export]]
List pairwise_moran_native_cpp(
    arma::mat data,
    const arma::mat& coords,
    double radius,
    double sigma = -1.0,
    int max_neighbors = 200,
    bool verbose = true
) {
    arma::uword n_genes = data.n_rows;
    arma::uword n_cells = data.n_cols;

    if (sigma < 0) {
        sigma = radius / 3.0;
    }

    if (verbose) {
        Rcpp::Rcout << "Pairwise Moran's I (Native C++)\n";
        Rcpp::Rcout << n_cells << " cells, " << n_genes << " genes\n";
        Rcpp::Rcout << "radius=" << radius << ", sigma=" << sigma << "\n";
    }

    // Step 1: Build KD-tree
    if (verbose) Rcpp::Rcout << "Building KD-tree...\n";

    PointCloud cloud(coords);
    KDTree index(2, cloud, nanoflann::KDTreeSingleIndexAdaptorParams(10));
    index.buildIndex();

    if (verbose) Rcpp::Rcout << "Finding neighbors and building weight matrix...\n";

    double radius_sq = radius * radius;
    double gaussian_factor = -1.0 / (2.0 * sigma * sigma);

    nanoflann::SearchParameters params;
    params.sorted = false;

    // Step 2: Build sparse weight matrix directly (memory efficient)
    // Process in chunks to avoid storing all neighbors in memory
    std::vector<arma::uword> all_rows, all_cols;
    std::vector<double> all_values;

    // Estimate and reserve (conservative estimate)
    size_t estimated_edges = std::min(static_cast<size_t>(n_cells) * max_neighbors,
                                       static_cast<size_t>(n_cells) * 100);
    all_rows.reserve(estimated_edges);
    all_cols.reserve(estimated_edges);
    all_values.reserve(estimated_edges);

    double total_weight = 0.0;
    arma::uword n_edges = 0;

    // Process cells in batches for memory efficiency
    const arma::uword batch_size = 10000;

    for (arma::uword batch_start = 0; batch_start < n_cells; batch_start += batch_size) {
        arma::uword batch_end = std::min(batch_start + batch_size, n_cells);

        if (verbose && batch_start > 0) {
            Rcpp::Rcout << "  Processing cells " << batch_start << "-" << batch_end << "...\n";
        }

        // Process this batch
        for (arma::uword i = batch_start; i < batch_end; i++) {
            double query_pt[2] = {coords(i, 0), coords(i, 1)};

            std::vector<nanoflann::ResultItem<size_t, double>> matches;
            index.radiusSearch(query_pt, radius_sq, matches, params);

            // Compute weights and row sum
            double row_sum = 0.0;
            std::vector<std::pair<arma::uword, double>> cell_weights;
            cell_weights.reserve(std::min(matches.size(), static_cast<size_t>(max_neighbors)));

            for (const auto& match : matches) {
                arma::uword j = static_cast<arma::uword>(match.first);
                if (j == i) continue;  // Skip self

                double w = std::exp(match.second * gaussian_factor);

                if (w > 1e-10) {
                    cell_weights.push_back(std::make_pair(j, w));
                    row_sum += w;

                    // Limit to max_neighbors if specified (0 = no limit)
                    if (max_neighbors > 0 &&
                        cell_weights.size() >= static_cast<size_t>(max_neighbors)) {
                        break;
                    }
                }
            }

            // Normalize and store
            if (row_sum > 0) {
                for (const auto& cw : cell_weights) {
                    all_rows.push_back(i);
                    all_cols.push_back(cw.first);
                    double w_norm = cw.second / row_sum;
                    all_values.push_back(w_norm);
                    total_weight += w_norm;
                    n_edges++;
                }
            }
        }

        // Allow R to check for interrupts
        Rcpp::checkUserInterrupt();
    }

    if (verbose) {
        Rcpp::Rcout << "W non-zeros: " << n_edges << "\n";
        Rcpp::Rcout << "S0: " << total_weight << "\n";
        Rcpp::Rcout << "Avg neighbors/cell: " << static_cast<double>(n_edges) / n_cells << "\n";
    }

    if (total_weight < 1e-10) {
        Rcpp::stop("Weight matrix has no non-zero weights");
    }

    // Build sparse matrix
    if (verbose) Rcpp::Rcout << "Building sparse matrix...\n";
    arma::sp_mat W;
    if (!all_rows.empty()) {
        arma::umat locations(2, all_rows.size());
        for (size_t k = 0; k < all_rows.size(); k++) {
            locations(0, k) = all_rows[k];
            locations(1, k) = all_cols[k];
        }
        W = arma::sp_mat(locations, arma::vec(all_values), n_cells, n_cells);
    } else {
        W = arma::sp_mat(n_cells, n_cells);
    }

    // Free memory from triplet vectors
    all_rows.clear(); all_rows.shrink_to_fit();
    all_cols.clear(); all_cols.shrink_to_fit();
    all_values.clear(); all_values.shrink_to_fit();

    // Step 3: Z-normalize data (population SD)
    if (verbose) Rcpp::Rcout << "Z-normalizing data...\n";

    for (arma::uword i = 0; i < n_genes; i++) {
        arma::rowvec row = data.row(i);
        double mean_val = arma::mean(row);
        double var_val = arma::var(row, 1);  // Population variance (N)
        double sd_val = std::sqrt(var_val);

        if (sd_val > 1e-10) {
            data.row(i) = (row - mean_val) / sd_val;
        } else {
            data.row(i).zeros();
        }
    }

    // Step 4: Compute Moran's I matrix
    if (verbose) Rcpp::Rcout << "Computing pairwise Moran's I...\n";

    // Spatial lag: lag = data * W^T
    arma::mat lag = data * W.t();

    // Moran matrix: M = data * lag^T / S0
    arma::mat moran = (data * lag.t()) / total_weight;

    if (verbose) Rcpp::Rcout << "Done.\n";

    return List::create(
        Named("moran") = moran,
        Named("weight_sum") = total_weight,
        Named("n_edges") = n_edges
    );
}


//' Compute Pairwise Moran's I Using Chunked Dense Approach (No Neighbor Limit)
//'
//' Similar to Python GPU implementation: processes senders in chunks,
//' computes dense distance matrix for each chunk, extracts sparse entries.
//' No neighbor limit - finds ALL neighbors within radius.
//'
//' Two-pass approach for memory efficiency:
//' - Pass 1: Compute row sums only (for normalization)
//' - Pass 2: Extract normalized triplets and build sparse matrix
//'
//' @param data Gene expression matrix (genes x cells)
//' @param coords Cell coordinates (n x 2 matrix)
//' @param radius Radius for neighbor search
//' @param sigma Gaussian sigma (default: radius/3)
//' @param chunk_size Number of cells per chunk (default: 1000)
//' @param verbose Print progress messages (default: TRUE)
//' @return List with moran matrix, weight_sum, and n_edges
//' @export
// [[Rcpp::export]]
List pairwise_moran_chunked_cpp(
    arma::mat data,
    const arma::mat& coords,
    double radius,
    double sigma = -1.0,
    int chunk_size = 1000,
    bool verbose = true
) {
    arma::uword n_genes = data.n_rows;
    arma::uword n_cells = data.n_cols;

    if (sigma < 0) {
        sigma = radius / 3.0;
    }

    if (verbose) {
        Rcpp::Rcout << "Pairwise Moran's I (Chunked Dense - No Neighbor Limit)\n";
        Rcpp::Rcout << n_cells << " cells, " << n_genes << " genes\n";
        Rcpp::Rcout << "radius=" << radius << ", sigma=" << sigma << "\n";
        Rcpp::Rcout << "chunk_size=" << chunk_size << "\n";
    }

    double radius_sq = radius * radius;
    double gaussian_factor = -1.0 / (2.0 * sigma * sigma);
    arma::uword n_chunks = (n_cells + chunk_size - 1) / chunk_size;

    // ============================================================
    // PASS 1: Compute row sums only (no neighbor storage)
    // ============================================================
    if (verbose) Rcpp::Rcout << "Pass 1: Computing row sums...\n";

    arma::vec row_sums(n_cells, arma::fill::zeros);
    arma::uword total_edges = 0;

    for (arma::uword chunk_idx = 0; chunk_idx < n_chunks; chunk_idx++) {
        arma::uword chunk_start = chunk_idx * chunk_size;
        arma::uword chunk_end = std::min(chunk_start + static_cast<arma::uword>(chunk_size), n_cells);
        arma::uword chunk_len = chunk_end - chunk_start;

        if (verbose && chunk_idx % 20 == 0) {
            Rcpp::Rcout << "  Chunk " << chunk_idx + 1 << "/" << n_chunks << "\n";
        }

        // Extract chunk coordinates
        arma::mat chunk_coords = coords.rows(chunk_start, chunk_end - 1);

        // Compute dense distance matrix: chunk_len x n_cells
        arma::mat diff_x = arma::repmat(chunk_coords.col(0), 1, n_cells) -
                           arma::repmat(coords.col(0).t(), chunk_len, 1);
        arma::mat diff_y = arma::repmat(chunk_coords.col(1), 1, n_cells) -
                           arma::repmat(coords.col(1).t(), chunk_len, 1);
        arma::mat dist_sq = diff_x % diff_x + diff_y % diff_y;

        // Compute Gaussian weights with radius mask
        arma::mat weights = arma::exp(dist_sq * gaussian_factor);
        weights.elem(arma::find(dist_sq > radius_sq)).zeros();

        // Zero out self-connections and small weights
        for (arma::uword i = 0; i < chunk_len; i++) {
            weights(i, chunk_start + i) = 0.0;
        }
        weights.elem(arma::find(weights < 1e-10)).zeros();

        // Accumulate row sums and count edges
        for (arma::uword i = 0; i < chunk_len; i++) {
            arma::uword global_i = chunk_start + i;
            for (arma::uword j = 0; j < n_cells; j++) {
                double w = weights(i, j);
                if (w > 0) {
                    row_sums(global_i) += w;
                    total_edges++;
                }
            }
        }

        Rcpp::checkUserInterrupt();
    }

    if (verbose) {
        Rcpp::Rcout << "Total edges: " << total_edges << "\n";
        Rcpp::Rcout << "Avg neighbors/cell: " << static_cast<double>(total_edges) / n_cells << "\n";
    }

    // ============================================================
    // PASS 2: Extract normalized triplets directly
    // ============================================================
    if (verbose) Rcpp::Rcout << "Pass 2: Extracting normalized triplets...\n";

    std::vector<arma::uword> all_rows, all_cols;
    std::vector<double> all_values;
    all_rows.reserve(total_edges);
    all_cols.reserve(total_edges);
    all_values.reserve(total_edges);

    double total_weight = 0.0;

    for (arma::uword chunk_idx = 0; chunk_idx < n_chunks; chunk_idx++) {
        arma::uword chunk_start = chunk_idx * chunk_size;
        arma::uword chunk_end = std::min(chunk_start + static_cast<arma::uword>(chunk_size), n_cells);
        arma::uword chunk_len = chunk_end - chunk_start;

        if (verbose && chunk_idx % 20 == 0) {
            Rcpp::Rcout << "  Chunk " << chunk_idx + 1 << "/" << n_chunks << "\n";
        }

        // Recompute weights for this chunk
        arma::mat chunk_coords = coords.rows(chunk_start, chunk_end - 1);

        arma::mat diff_x = arma::repmat(chunk_coords.col(0), 1, n_cells) -
                           arma::repmat(coords.col(0).t(), chunk_len, 1);
        arma::mat diff_y = arma::repmat(chunk_coords.col(1), 1, n_cells) -
                           arma::repmat(coords.col(1).t(), chunk_len, 1);
        arma::mat dist_sq = diff_x % diff_x + diff_y % diff_y;

        arma::mat weights = arma::exp(dist_sq * gaussian_factor);
        weights.elem(arma::find(dist_sq > radius_sq)).zeros();

        for (arma::uword i = 0; i < chunk_len; i++) {
            weights(i, chunk_start + i) = 0.0;
        }
        weights.elem(arma::find(weights < 1e-10)).zeros();

        // Extract and normalize
        for (arma::uword i = 0; i < chunk_len; i++) {
            arma::uword global_i = chunk_start + i;
            double norm = row_sums(global_i) > 0 ? row_sums(global_i) : 1.0;

            for (arma::uword j = 0; j < n_cells; j++) {
                double w = weights(i, j);
                if (w > 0) {
                    all_rows.push_back(global_i);
                    all_cols.push_back(j);
                    double w_norm = w / norm;
                    all_values.push_back(w_norm);
                    total_weight += w_norm;
                }
            }
        }

        Rcpp::checkUserInterrupt();
    }

    if (total_weight < 1e-10) {
        Rcpp::stop("Weight matrix has no non-zero weights");
    }

    // Build sparse matrix
    if (verbose) Rcpp::Rcout << "Building sparse matrix...\n";
    arma::sp_mat W;
    if (!all_rows.empty()) {
        arma::umat locations(2, all_rows.size());
        for (size_t k = 0; k < all_rows.size(); k++) {
            locations(0, k) = all_rows[k];
            locations(1, k) = all_cols[k];
        }
        W = arma::sp_mat(locations, arma::vec(all_values), n_cells, n_cells);
    } else {
        W = arma::sp_mat(n_cells, n_cells);
    }

    // Free triplet vectors
    all_rows.clear(); all_rows.shrink_to_fit();
    all_cols.clear(); all_cols.shrink_to_fit();
    all_values.clear(); all_values.shrink_to_fit();

    // Z-normalize data (population SD)
    if (verbose) Rcpp::Rcout << "Z-normalizing data...\n";

    for (arma::uword i = 0; i < n_genes; i++) {
        arma::rowvec row = data.row(i);
        double mean_val = arma::mean(row);
        double var_val = arma::var(row, 1);
        double sd_val = std::sqrt(var_val);

        if (sd_val > 1e-10) {
            data.row(i) = (row - mean_val) / sd_val;
        } else {
            data.row(i).zeros();
        }
    }

    // Compute Moran's I matrix
    if (verbose) Rcpp::Rcout << "Computing pairwise Moran's I...\n";

    arma::mat lag = data * W.t();
    arma::mat moran = (data * lag.t()) / total_weight;

    if (verbose) Rcpp::Rcout << "Done.\n";

    return List::create(
        Named("moran") = moran,
        Named("weight_sum") = total_weight,
        Named("n_edges") = total_edges
    );
}


//' Compute Pairwise Moran's I with Streaming Lag (No Neighbor Limit, No W Storage)
//'
//' Best approach for CPU: uses KD-tree for O(n log n) neighbor search,
//' computes spatial lag on-the-fly without storing weight matrix W.
//' Unlimited neighbors without memory explosion.
//'
//' Memory usage: O(n_genes × n_cells) for data + lag matrices only.
//' Uses OpenMP for parallel processing across cells when available.
//'
//' @param data Gene expression matrix (genes x cells)
//' @param coords Cell coordinates (n x 2 matrix)
//' @param radius Radius for neighbor search
//' @param sigma Gaussian sigma (default: radius/3)
//' @param n_threads Number of threads (default: 0 = auto, ignored if OpenMP unavailable)
//' @param verbose Print progress messages (default: TRUE)
//' @return List with moran matrix, weight_sum, and n_edges
//' @export
// [[Rcpp::export]]
List pairwise_moran_streaming_cpp(
    arma::mat data,
    const arma::mat& coords,
    double radius,
    double sigma = -1.0,
    int n_threads = 0,
    bool verbose = true
) {
    arma::uword n_genes = data.n_rows;
    arma::uword n_cells = data.n_cols;

    if (sigma < 0) {
        sigma = radius / 3.0;
    }

    // Set number of threads (OpenMP only)
    #ifdef _OPENMP
    if (n_threads > 0) {
        omp_set_num_threads(n_threads);
    }
    int actual_threads = omp_get_max_threads();
    #endif

    if (verbose) {
        Rcpp::Rcout << "Pairwise Moran's I (Streaming";
        #ifdef _OPENMP
        Rcpp::Rcout << " - " << actual_threads << " threads";
        #else
        Rcpp::Rcout << " - single thread";
        #endif
        Rcpp::Rcout << ")\n";
        Rcpp::Rcout << n_cells << " cells, " << n_genes << " genes\n";
        Rcpp::Rcout << "radius=" << radius << ", sigma=" << sigma << "\n";
    }

    // Step 1: Z-normalize data first (population SD) - parallelized
    if (verbose) Rcpp::Rcout << "Z-normalizing data...\n";

    #pragma omp parallel for schedule(static)
    for (arma::uword i = 0; i < n_genes; i++) {
        arma::rowvec row = data.row(i);
        double mean_val = arma::mean(row);
        double var_val = arma::var(row, 1);
        double sd_val = std::sqrt(var_val);

        if (sd_val > 1e-10) {
            data.row(i) = (row - mean_val) / sd_val;
        } else {
            data.row(i).zeros();
        }
    }

    // Step 2: Build KD-tree
    if (verbose) Rcpp::Rcout << "Building KD-tree...\n";

    PointCloud cloud(coords);
    KDTree index(2, cloud, nanoflann::KDTreeSingleIndexAdaptorParams(10));
    index.buildIndex();

    double radius_sq = radius * radius;
    double gaussian_factor = -1.0 / (2.0 * sigma * sigma);

    // Step 3: Compute spatial lag ON-THE-FLY
    // lag[:, i] = sum_j(w_ij * data[:, j]) for all j in radius
    // where w_ij is row-normalized
    if (verbose) Rcpp::Rcout << "Computing spatial lag...\n";

    arma::mat lag(n_genes, n_cells, arma::fill::zeros);

    // Accumulators for statistics
    double total_weight = 0.0;
    arma::uword total_edges = 0;

    // Progress tracking
    const arma::uword report_interval = 20000;
    arma::uword last_report = 0;

    #ifdef _OPENMP
    // OpenMP parallel version
    #pragma omp parallel reduction(+:total_weight, total_edges)
    {
        // Thread-local search parameters
        nanoflann::SearchParameters params;
        params.sorted = false;

        #pragma omp for schedule(dynamic, 100)
        for (arma::uword i = 0; i < n_cells; i++) {
            double query_pt[2] = {coords(i, 0), coords(i, 1)};

            // KD-tree radius search - finds ALL neighbors (no limit)
            std::vector<nanoflann::ResultItem<size_t, double>> matches;
            index.radiusSearch(query_pt, radius_sq, matches, params);

            // First pass: compute row sum for normalization
            double row_sum = 0.0;
            std::vector<std::pair<arma::uword, double>> neighbors;
            neighbors.reserve(matches.size());

            for (const auto& match : matches) {
                arma::uword j = static_cast<arma::uword>(match.first);
                if (j == i) continue;

                double w = std::exp(match.second * gaussian_factor);
                if (w > 1e-10) {
                    neighbors.push_back(std::make_pair(j, w));
                    row_sum += w;
                }
            }

            // Second pass: compute weighted sum for spatial lag
            if (row_sum > 0) {
                arma::vec lag_i(n_genes, arma::fill::zeros);

                for (const auto& neighbor : neighbors) {
                    arma::uword j = neighbor.first;
                    double w_norm = neighbor.second / row_sum;
                    lag_i += w_norm * data.col(j);
                    total_weight += w_norm;
                    total_edges++;
                }

                lag.col(i) = lag_i;
            }
        }
    }
    #else
    // Serial version with progress reporting and interrupt checking
    nanoflann::SearchParameters params;
    params.sorted = false;

    for (arma::uword i = 0; i < n_cells; i++) {
        // Progress reporting
        if (verbose && i >= last_report + report_interval) {
            Rcpp::Rcout << "  Processed " << i << "/" << n_cells << " cells ("
                        << (100 * i / n_cells) << "%)\n";
            last_report = i;
        }

        // Check for R interrupt periodically
        if (i % 1000 == 0) {
            Rcpp::checkUserInterrupt();
        }

        double query_pt[2] = {coords(i, 0), coords(i, 1)};

        std::vector<nanoflann::ResultItem<size_t, double>> matches;
        index.radiusSearch(query_pt, radius_sq, matches, params);

        double row_sum = 0.0;
        std::vector<std::pair<arma::uword, double>> neighbors;
        neighbors.reserve(matches.size());

        for (const auto& match : matches) {
            arma::uword j = static_cast<arma::uword>(match.first);
            if (j == i) continue;

            double w = std::exp(match.second * gaussian_factor);
            if (w > 1e-10) {
                neighbors.push_back(std::make_pair(j, w));
                row_sum += w;
            }
        }

        if (row_sum > 0) {
            arma::vec lag_i(n_genes, arma::fill::zeros);

            for (const auto& neighbor : neighbors) {
                arma::uword j = neighbor.first;
                double w_norm = neighbor.second / row_sum;
                lag_i += w_norm * data.col(j);
                total_weight += w_norm;
                total_edges++;
            }

            lag.col(i) = lag_i;
        }
    }
    #endif

    // Check for R interrupt after processing
    Rcpp::checkUserInterrupt();

    if (verbose) {
        Rcpp::Rcout << "Total edges: " << total_edges << "\n";
        Rcpp::Rcout << "Avg neighbors/cell: " << static_cast<double>(total_edges) / n_cells << "\n";
        Rcpp::Rcout << "S0: " << total_weight << "\n";
    }

    if (total_weight < 1e-10) {
        Rcpp::stop("No valid neighbors found within radius");
    }

    // Step 4: Compute Moran's I matrix
    // M = Z * lag^T / S0
    if (verbose) Rcpp::Rcout << "Computing pairwise Moran's I...\n";

    arma::mat moran = (data * lag.t()) / total_weight;

    if (verbose) Rcpp::Rcout << "Done.\n";

    return List::create(
        Named("moran") = moran,
        Named("weight_sum") = total_weight,
        Named("n_edges") = total_edges
    );
}


//' Compute Pairwise Moran's I with Chunked Dense BLAS (Best for Large Radii)
//'
//' Uses chunked dense matrix operations with BLAS for efficient computation
//' when many neighbors per cell. Processes cells in batches to limit memory.
//'
//' For each chunk of cells, computes dense weight matrix and uses optimized
//' BLAS matrix multiplication for spatial lag. Much faster than streaming
//' when density > 5% (i.e., > n/20 neighbors per cell).
//'
//' @param data Gene expression matrix (genes x cells)
//' @param coords Cell coordinates (n x 2 matrix)
//' @param radius Radius for neighbor search
//' @param sigma Gaussian sigma (default: radius/3)
//' @param chunk_size Cells per chunk (default: 500, tune for memory)
//' @param verbose Print progress messages (default: TRUE)
//' @return List with moran matrix, weight_sum, and n_edges
//' @export
// [[Rcpp::export]]
List pairwise_moran_dense_cpp(
    arma::mat data,
    const arma::mat& coords,
    double radius,
    double sigma = -1.0,
    int chunk_size = 500,
    bool verbose = true
) {
    arma::uword n_genes = data.n_rows;
    arma::uword n_cells = data.n_cols;

    if (sigma < 0) {
        sigma = radius / 3.0;
    }

    if (verbose) {
        Rcpp::Rcout << "Pairwise Moran's I (Chunked Dense BLAS)\n";
        Rcpp::Rcout << n_cells << " cells, " << n_genes << " genes\n";
        Rcpp::Rcout << "radius=" << radius << ", sigma=" << sigma << "\n";
        Rcpp::Rcout << "chunk_size=" << chunk_size << "\n";
    }

    // Step 1: Z-normalize data (population SD)
    if (verbose) Rcpp::Rcout << "Z-normalizing data...\n";

    #pragma omp parallel for schedule(static)
    for (arma::uword i = 0; i < n_genes; i++) {
        arma::rowvec row = data.row(i);
        double mean_val = arma::mean(row);
        double var_val = arma::var(row, 1);
        double sd_val = std::sqrt(var_val);

        if (sd_val > 1e-10) {
            data.row(i) = (row - mean_val) / sd_val;
        } else {
            data.row(i).zeros();
        }
    }

    // Transpose data for efficient column access: data_t is n_cells x n_genes
    arma::mat data_t = data.t();

    double radius_sq = radius * radius;
    double gaussian_factor = -1.0 / (2.0 * sigma * sigma);

    // Step 2: Compute spatial lag using chunked dense BLAS
    if (verbose) Rcpp::Rcout << "Computing spatial lag (chunked dense)...\n";

    arma::mat lag(n_genes, n_cells, arma::fill::zeros);

    arma::uword n_chunks = (n_cells + chunk_size - 1) / chunk_size;
    double total_weight = 0.0;
    arma::uword total_edges = 0;

    // Row sums for normalization (need to compute first)
    arma::vec row_sums(n_cells, arma::fill::zeros);

    for (arma::uword chunk_idx = 0; chunk_idx < n_chunks; chunk_idx++) {
        arma::uword chunk_start = chunk_idx * chunk_size;
        arma::uword chunk_end = std::min(chunk_start + static_cast<arma::uword>(chunk_size), n_cells);
        arma::uword chunk_len = chunk_end - chunk_start;

        if (verbose && (chunk_idx % 10 == 0 || chunk_idx == n_chunks - 1)) {
            Rcpp::Rcout << "  Chunk " << chunk_idx + 1 << "/" << n_chunks
                        << " (cells " << chunk_start << "-" << chunk_end - 1 << ")\n";
        }

        // Compute pairwise distances: chunk_len x n_cells
        arma::mat chunk_x = coords.rows(chunk_start, chunk_end - 1).col(0);
        arma::mat chunk_y = coords.rows(chunk_start, chunk_end - 1).col(1);

        // Broadcast difference computation using outer product pattern
        arma::mat diff_x = arma::repmat(chunk_x, 1, n_cells) -
                           arma::repmat(coords.col(0).t(), chunk_len, 1);
        arma::mat diff_y = arma::repmat(chunk_y, 1, n_cells) -
                           arma::repmat(coords.col(1).t(), chunk_len, 1);

        // Squared distances
        arma::mat dist_sq = diff_x % diff_x + diff_y % diff_y;

        // Gaussian weights with radius cutoff
        arma::mat W_chunk = arma::exp(dist_sq * gaussian_factor);
        W_chunk.elem(arma::find(dist_sq > radius_sq)).zeros();

        // Zero out self-connections
        for (arma::uword i = 0; i < chunk_len; i++) {
            W_chunk(i, chunk_start + i) = 0.0;
        }

        // Threshold small weights
        W_chunk.elem(arma::find(W_chunk < 1e-10)).zeros();

        // Compute row sums and count edges
        for (arma::uword i = 0; i < chunk_len; i++) {
            double rs = arma::accu(W_chunk.row(i));
            row_sums(chunk_start + i) = rs;

            // Count non-zero entries for this row
            arma::uword nnz = arma::accu(W_chunk.row(i) > 0);
            total_edges += nnz;
        }

        // Row normalize
        for (arma::uword i = 0; i < chunk_len; i++) {
            if (row_sums(chunk_start + i) > 0) {
                W_chunk.row(i) /= row_sums(chunk_start + i);
            }
        }

        // Compute lag contribution using BLAS: (chunk_len x n_cells) × (n_cells x n_genes)
        // lag_chunk = W_chunk × data_t, result is chunk_len x n_genes
        arma::mat lag_chunk = W_chunk * data_t;

        // Store in lag matrix (transposed)
        lag.cols(chunk_start, chunk_end - 1) = lag_chunk.t();

        // Accumulate total weight
        total_weight += arma::accu(W_chunk);

        Rcpp::checkUserInterrupt();
    }

    if (verbose) {
        Rcpp::Rcout << "Total edges: " << total_edges << "\n";
        Rcpp::Rcout << "Avg neighbors/cell: " << static_cast<double>(total_edges) / n_cells << "\n";
        Rcpp::Rcout << "S0: " << total_weight << "\n";
    }

    if (total_weight < 1e-10) {
        Rcpp::stop("No valid neighbors found within radius");
    }

    // Step 3: Compute Moran's I matrix using BLAS
    // M = data × lag^T / S0
    if (verbose) Rcpp::Rcout << "Computing pairwise Moran's I...\n";

    arma::mat moran = (data * lag.t()) / total_weight;

    if (verbose) Rcpp::Rcout << "Done.\n";

    return List::create(
        Named("moran") = moran,
        Named("weight_sum") = total_weight,
        Named("n_edges") = total_edges
    );
}
