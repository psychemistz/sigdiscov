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
