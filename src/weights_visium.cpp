// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <vector>

using namespace Rcpp;

// =============================================================================
// VISIUM-SPECIFIC WEIGHT MATRIX FUNCTIONS
//
// For Visium data, we use binary weights (all neighbors within radius get equal
// weight) because:
// - Grid structure has regular spacing (~100μm center-to-center)
// - Spot density is relatively low (~4K spots per sample)
// - Distance-decay effects are less pronounced at this resolution
// =============================================================================

//' Create Binary Weight Matrix for Visium
//'
//' Computes a row-normalized binary weight matrix where all spots within
//' the specified radius receive equal weight.
//'
//' @param coords Spot coordinates matrix (n x 2)
//' @param radius Distance threshold for neighbor definition
//' @param include_self Include self-connections (default: FALSE)
//'
//' @return Sparse weight matrix (n x n), row-normalized so each row sums to 1
//'
//' @details
//' The weight is computed as:
//' \deqn{w_{ij} = \frac{1}{n_i} \text{ if } d(i,j) \leq r, \text{ else } 0}
//'
//' where \eqn{n_i} is the number of neighbors for spot i.
//'
//' @export
// [[Rcpp::export]]
arma::sp_mat create_binary_weights_cpp(
    const arma::mat& coords,
    double radius,
    bool include_self = false
) {
    int n = coords.n_rows;
    double radius_sq = radius * radius;

    // Triplet storage for sparse matrix
    std::vector<arma::uword> rows, cols;
    std::vector<double> values;

    for (int i = 0; i < n; i++) {
        double xi = coords(i, 0);
        double yi = coords(i, 1);

        std::vector<int> neighbors;

        for (int j = 0; j < n; j++) {
            if (!include_self && i == j) continue;

            double dx = xi - coords(j, 0);
            double dy = yi - coords(j, 1);
            double dist_sq = dx * dx + dy * dy;

            if (dist_sq <= radius_sq) {
                neighbors.push_back(j);
            }
        }

        // Row normalize: each neighbor gets 1/n_neighbors
        int n_neighbors = neighbors.size();
        if (n_neighbors > 0) {
            double weight = 1.0 / n_neighbors;
            for (int j : neighbors) {
                rows.push_back(i);
                cols.push_back(j);
                values.push_back(weight);
            }
        }
    }

    // Build sparse matrix
    if (rows.empty()) {
        return arma::sp_mat(n, n);
    }

    arma::umat locations(2, rows.size());
    for (size_t k = 0; k < rows.size(); k++) {
        locations(0, k) = rows[k];
        locations(1, k) = cols[k];
    }

    return arma::sp_mat(locations, arma::vec(values), n, n);
}

//' Create Ring Weight Matrix for Visium
//'
//' Computes a weight matrix where only spots in a ring (annulus) between
//' inner and outer radii receive weight. Useful for excluding immediate
//' neighbors or analyzing specific distance bands.
//'
//' @param coords Spot coordinates matrix (n x 2)
//' @param inner_radius Inner radius (exclusive)
//' @param outer_radius Outer radius (inclusive)
//'
//' @return Sparse ring weight matrix (n x n), row-normalized
//'
//' @details
//' Weight is assigned to spots where \eqn{r_{inner} < d(i,j) \leq r_{outer}}.
//'
//' @export
// [[Rcpp::export]]
arma::sp_mat create_ring_weights_cpp(
    const arma::mat& coords,
    double inner_radius,
    double outer_radius
) {
    int n = coords.n_rows;
    double inner_sq = inner_radius * inner_radius;
    double outer_sq = outer_radius * outer_radius;

    std::vector<arma::uword> rows, cols;
    std::vector<double> values;

    for (int i = 0; i < n; i++) {
        double xi = coords(i, 0);
        double yi = coords(i, 1);

        std::vector<int> neighbors;

        for (int j = 0; j < n; j++) {
            if (i == j) continue;

            double dx = xi - coords(j, 0);
            double dy = yi - coords(j, 1);
            double dist_sq = dx * dx + dy * dy;

            if (dist_sq > inner_sq && dist_sq <= outer_sq) {
                neighbors.push_back(j);
            }
        }

        int n_neighbors = neighbors.size();
        if (n_neighbors > 0) {
            double weight = 1.0 / n_neighbors;
            for (int j : neighbors) {
                rows.push_back(i);
                cols.push_back(j);
                values.push_back(weight);
            }
        }
    }

    if (rows.empty()) {
        return arma::sp_mat(n, n);
    }

    arma::umat locations(2, rows.size());
    for (size_t k = 0; k < rows.size(); k++) {
        locations(0, k) = rows[k];
        locations(1, k) = cols[k];
    }

    return arma::sp_mat(locations, arma::vec(values), n, n);
}

//' Create Directional Weight Matrix (Sender to Receiver)
//'
//' Computes a weight matrix for directional analysis where senders and
//' receivers are different sets of spots. Used for expression-based
//' sender/receiver splitting in Visium analysis.
//'
//' @param sender_coords Sender spot coordinates (n_s x 2)
//' @param receiver_coords Receiver spot coordinates (n_r x 2)
//' @param radius Distance threshold
//'
//' @return Sparse weight matrix (n_senders x n_receivers), row-normalized
//'
//' @details
//' For directional mode, the weight matrix is not square. Each sender spot
//' gets connections to receiver spots within the radius, with row normalization
//' so that weights sum to 1 per sender.
//'
//' @export
// [[Rcpp::export]]
arma::sp_mat create_directional_weights_cpp(
    const arma::mat& sender_coords,
    const arma::mat& receiver_coords,
    double radius
) {
    int n_senders = sender_coords.n_rows;
    int n_receivers = receiver_coords.n_rows;
    double radius_sq = radius * radius;

    std::vector<arma::uword> rows, cols;
    std::vector<double> values;

    for (int i = 0; i < n_senders; i++) {
        double xi = sender_coords(i, 0);
        double yi = sender_coords(i, 1);

        std::vector<int> neighbors;

        for (int j = 0; j < n_receivers; j++) {
            double dx = xi - receiver_coords(j, 0);
            double dy = yi - receiver_coords(j, 1);
            double dist_sq = dx * dx + dy * dy;

            if (dist_sq <= radius_sq) {
                neighbors.push_back(j);
            }
        }

        int n_neighbors = neighbors.size();
        if (n_neighbors > 0) {
            double weight = 1.0 / n_neighbors;
            for (int j : neighbors) {
                rows.push_back(i);
                cols.push_back(j);
                values.push_back(weight);
            }
        }
    }

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
