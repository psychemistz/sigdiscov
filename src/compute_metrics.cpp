// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>

using namespace Rcpp;

// =============================================================================
// CORE METRIC COMPUTATION FUNCTIONS
//
// These functions compute spatial correlation metrics from PRE-COMPUTED
// spatial lags. The key insight is that the expensive operation (W * z_g)
// is done once, then metrics can be computed cheaply.
//
// Mathematical formulations:
//   Bivariate Moran's I: I = z_f' * lag_g / n
//   I_ND (cosine similarity): I_ND = z_f' * lag_g / (||z_f|| * ||lag_g||)
// =============================================================================

//' Compute Bivariate Moran's I from Pre-computed Spatial Lag
//'
//' Calculates Moran's I statistic using the formula:
//'   I = z_f' * lag_g / n
//'
//' where lag_g = W * z_g has been pre-computed.
//'
//' @param z_f Standardized factor expression vector (length n)
//' @param lag_g Pre-computed spatial lag of gene (W * z_g), same length as z_f
//'
//' @return Scalar Moran's I value (unbounded)
//'
//' @details
//' Moran's I measures spatial autocorrelation. Positive values indicate
//' spatial clustering (similar values near each other), negative values
//' indicate dispersion (dissimilar values near each other).
//'
//' @export
// [[Rcpp::export]]
double compute_moran_from_lag_cpp(const arma::vec& z_f, const arma::vec& lag_g) {
    // I = z_f' * lag_g / n
    // This is the bivariate Moran's I formula with pre-computed spatial lag
    return arma::dot(z_f, lag_g) / static_cast<double>(z_f.n_elem);
}

//' Compute I_ND (Cosine Similarity) from Pre-computed Spatial Lag
//'
//' Calculates the normalized directional Moran's I (I_ND) using:
//'   I_ND = z_f' * lag_g / (||z_f|| * ||lag_g||)
//'
//' This is equivalent to the cosine similarity between z_f and lag_g.
//'
//' @param z_f Standardized factor expression vector (length n)
//' @param lag_g Pre-computed spatial lag of gene (W * z_g)
//'
//' @return Scalar I_ND value bounded between -1 and 1, or NA if norms are too small
//'
//' @details
//' I_ND is the cosine of the angle between z_f and the spatial lag of z_g.
//' It is bounded between -1 and 1, making it interpretable as a correlation:
//'   +1: Perfect positive spatial association
//'    0: No spatial association
//'   -1: Perfect negative spatial association
//'
//' Returns NA if either vector has near-zero norm (constant or zero expression).
//'
//' @export
// [[Rcpp::export]]
double compute_ind_from_lag_cpp(const arma::vec& z_f, const arma::vec& lag_g) {
    // Compute norms
    double norm_f = arma::norm(z_f, 2);
    double norm_lag = arma::norm(lag_g, 2);

    // Check for degenerate cases (constant or zero expression)
    if (norm_f < 1e-10 || norm_lag < 1e-10) {
        return NA_REAL;
    }

    // I_ND = dot(z_f, lag_g) / (||z_f|| * ||lag_g||)
    double I_ND = arma::dot(z_f, lag_g) / (norm_f * norm_lag);

    return I_ND;
}

//' Batch Compute Metrics for All Genes
//'
//' Efficiently computes Moran's I or I_ND for one factor against all genes
//' using pre-computed spatial lag matrix.
//'
//' This is the workhorse function for genome-wide analysis: compute the
//' spatial lag matrix once (W * Z_g), then call this function to get
//' metrics for all genes in one vectorized operation.
//'
//' @param z_f Standardized factor expression (length n)
//' @param lag_G Pre-computed spatial lag matrix (n x n_genes)
//'        Each column is W * z_g for one gene
//' @param metric String: "moran" or "ind"
//'
//' @return Numeric vector of metric values (length n_genes)
//'
//' @details
//' For "moran": I_g = z_f' * lag_G_g / n for each gene g
//' For "ind": I_g = z_f' * lag_G_g / (||z_f|| * ||lag_G_g||) for each gene g
//'
//' The function is optimized using vectorized operations where possible.
//'
//' @export
// [[Rcpp::export]]
arma::vec compute_metric_batch_cpp(
    const arma::vec& z_f,
    const arma::mat& lag_G,
    const std::string& metric
) {
    int n_genes = lag_G.n_cols;
    int n = z_f.n_elem;
    arma::vec result(n_genes);

    if (metric == "moran") {
        // Moran's I: I = z_f' * lag_g / n
        // Vectorized: result = (z_f' * lag_G) / n
        // z_f' * lag_G gives a 1 x n_genes row vector
        arma::rowvec products = z_f.t() * lag_G;
        result = products.t() / static_cast<double>(n);

    } else if (metric == "ind") {
        // I_ND: I = z_f' * lag_g / (||z_f|| * ||lag_g||)
        double norm_f = arma::norm(z_f, 2);

        // Check for degenerate factor
        if (norm_f < 1e-10) {
            result.fill(NA_REAL);
            return result;
        }

        // Normalize factor once
        arma::vec f_normalized = z_f / norm_f;

        // Compute dot products: f_normalized' * lag_G (1 x n_genes)
        arma::rowvec correlations = f_normalized.t() * lag_G;

        // Compute column norms of lag_G: ||lag_G[,g]|| for each gene
        // sum(lag_G.^2, 0) gives row vector of column sums of squares
        arma::rowvec lag_norms = arma::sqrt(arma::sum(arma::square(lag_G), 0));

        // Compute I_ND for each gene
        for (int g = 0; g < n_genes; g++) {
            if (lag_norms(g) > 1e-10) {
                result(g) = correlations(g) / lag_norms(g);
            } else {
                result(g) = NA_REAL;
            }
        }
    } else {
        Rcpp::stop("Unknown metric: '%s'. Use 'moran' or 'ind'.", metric.c_str());
    }

    return result;
}

// =============================================================================
// STANDARDIZATION FUNCTIONS
//
// Gene-wise standardization (z-score) is a preprocessing step for spatial
// correlation analysis. These functions provide efficient C++ implementations.
// =============================================================================

//' Standardize a Vector (Z-score)
//'
//' @param x Numeric vector
//' @return Standardized vector (mean=0, sd=1)
//'
//' @details
//' Computes (x - mean(x)) / sd(x). Uses population standard deviation (N, not N-1).
//' If sd is near zero (constant vector), returns zeros.
//'
//' @keywords internal
// [[Rcpp::export]]
arma::vec standardize_vec_cpp(const arma::vec& x) {
    double mean_x = arma::mean(x);
    double var_x = arma::var(x, 1);  // 1 = normalize by N (population variance)
    double sd_x = std::sqrt(var_x);

    if (sd_x < 1e-10) {
        // Constant vector - return zeros
        return arma::vec(x.n_elem, arma::fill::zeros);
    }

    return (x - mean_x) / sd_x;
}

//' Standardize Matrix Row-wise (Gene-wise Z-score)
//'
//' @param X Numeric matrix (genes x observations)
//' @return Row-wise standardized matrix
//'
//' @details
//' Each row is independently standardized to have mean 0 and sd 1.
//' This is the standard preprocessing for expression matrices where
//' genes are rows and observations (spots/cells) are columns.
//'
//' @keywords internal
// [[Rcpp::export]]
arma::mat standardize_matrix_cpp(const arma::mat& X) {
    int n_rows = X.n_rows;
    arma::mat result = X;

    for (int i = 0; i < n_rows; i++) {
        double mean_i = arma::mean(X.row(i));
        double var_i = arma::var(X.row(i), 1);  // Population variance
        double sd_i = std::sqrt(var_i);

        if (sd_i > 1e-10) {
            result.row(i) = (X.row(i) - mean_i) / sd_i;
        } else {
            result.row(i).zeros();
        }
    }

    return result;
}
