// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <random>
#include <algorithm>
#include <numeric>

using namespace Rcpp;

// =============================================================================
// PERMUTATION TESTS FOR SPATIAL CORRELATION
//
// These functions implement efficient permutation-based significance testing
// for spatial correlation metrics (Moran's I and I_ND).
//
// Key optimization: The expensive spatial lag computation (W * z_g) is done
// ONCE, then factor expression is permuted. This exploits the fact that under
// the null hypothesis of no spatial association, factor expression values are
// randomly distributed across locations.
//
// For batch testing: Permutation indices are pre-generated ONCE and reused
// across all genes, ensuring:
//   1. Computational efficiency (single random number generation)
//   2. Consistent null distribution across genes
// =============================================================================

//' Permutation Test for Spatial Correlation (Single Gene)
//'
//' Tests H0: No spatial association between factor and gene expression.
//' The test permutes factor expression while keeping the spatial lag fixed.
//'
//' @param z_f Standardized factor expression vector (length n)
//' @param lag_g Pre-computed spatial lag of gene (W * z_g), same length as z_f
//' @param metric String: "moran" for bivariate Moran's I, "ind" for I_ND
//' @param n_perm Number of permutations (default: 999)
//'
//' @return List with components:
//'   - I_obs: Observed statistic value
//'   - p_value: Two-sided p-value with continuity correction
//'   - z_score: Z-score relative to null distribution
//'   - null_mean: Mean of null distribution
//'   - null_sd: Standard deviation of null distribution
//'   - n_perm: Number of permutations performed
//'
//' @details
//' The permutation strategy is:
//' 1. Compute observed statistic from z_f and lag_g
//' 2. For b = 1 to n_perm:
//'    a. Randomly permute z_f
//'    b. Compute permuted statistic (lag_g stays fixed)
//' 3. p-value = (1 + count(|I_perm| >= |I_obs|)) / (n_perm + 1)
//'
//' The continuity correction (+1 in numerator and denominator) ensures
//' p-values are never exactly 0.
//'
//' @export
// [[Rcpp::export]]
List permutation_test_cpp(
    const arma::vec& z_f,
    const arma::vec& lag_g,
    const std::string& metric,
    int n_perm
) {
    int n = z_f.n_elem;

    // Pre-compute norms for efficiency
    double norm_f = arma::norm(z_f, 2);
    double norm_lag = arma::norm(lag_g, 2);

    // Check for degenerate cases (for I_ND)
    if (metric == "ind" && (norm_f < 1e-10 || norm_lag < 1e-10)) {
        return List::create(
            Named("I_obs") = NA_REAL,
            Named("p_value") = NA_REAL,
            Named("z_score") = NA_REAL,
            Named("null_mean") = NA_REAL,
            Named("null_sd") = NA_REAL,
            Named("n_perm") = n_perm
        );
    }

    // Compute observed statistic
    double I_obs;
    if (metric == "moran") {
        // Moran's I: I = z_f' * lag_g / n
        I_obs = arma::dot(z_f, lag_g) / static_cast<double>(n);
    } else {
        // I_ND: cosine similarity
        I_obs = arma::dot(z_f, lag_g) / (norm_f * norm_lag);
    }

    // Permutation test
    int count_extreme = 0;
    double sum_perm = 0.0;
    double sum_perm_sq = 0.0;

    for (int b = 0; b < n_perm; b++) {
        // Generate random permutation
        arma::uvec perm_idx = arma::randperm(n);
        arma::vec z_f_perm = z_f.elem(perm_idx);

        // Compute permuted statistic
        double I_perm;
        if (metric == "moran") {
            I_perm = arma::dot(z_f_perm, lag_g) / static_cast<double>(n);
        } else {
            // For standardized data, norm(z_f_perm) == norm(z_f)
            // because permutation preserves sum of squares
            I_perm = arma::dot(z_f_perm, lag_g) / (norm_f * norm_lag);
        }

        // Accumulate for null distribution statistics
        sum_perm += I_perm;
        sum_perm_sq += I_perm * I_perm;

        // Two-sided test: count permutations with |I_perm| >= |I_obs|
        if (std::abs(I_perm) >= std::abs(I_obs)) {
            count_extreme++;
        }
    }

    // Compute p-value with continuity correction
    double p_value = (1.0 + count_extreme) / (n_perm + 1.0);

    // Compute null distribution statistics
    double null_mean = sum_perm / n_perm;
    double null_var = sum_perm_sq / n_perm - null_mean * null_mean;
    double null_sd = std::sqrt(std::max(null_var, 1e-10));

    // Z-score
    double z_score = (I_obs - null_mean) / null_sd;

    return List::create(
        Named("I_obs") = I_obs,
        Named("p_value") = p_value,
        Named("z_score") = z_score,
        Named("null_mean") = null_mean,
        Named("null_sd") = null_sd,
        Named("n_perm") = n_perm
    );
}

//' Batch Permutation Test for Multiple Genes
//'
//' Efficiently tests all genes against one factor. Pre-generates permutation
//' indices once and reuses them across all genes.
//'
//' @param z_f Standardized factor expression vector (length n_obs)
//' @param Z_g Gene expression matrix (n_obs x n_genes), genes in columns
//' @param W Sparse weight matrix (n_obs x n_obs) for computing spatial lags
//' @param metric String: "moran" or "ind"
//' @param n_perm Number of permutations (default: 999)
//' @param weight_sum Sum of weights for normalization. If <= 0, uses n_obs.
//'
//' @return DataFrame with columns:
//'   - gene_idx: Gene index (1-based)
//'   - I_obs: Observed statistic (Moran's I or I_ND)
//'   - p_value: Two-sided p-value
//'   - z_score: Z-score relative to null
//'   - null_sd: Standard deviation of the null distribution
//'
//' @details
//' Key optimizations:
//' 1. Permutation indices are generated ONCE and reused for all genes
//' 2. Spatial lag is computed once per gene
//' 3. User interrupt is checked every 100 genes
//'
//' This function is the workhorse for genome-wide significance testing.
//' For a typical analysis with 20,000 genes and 999 permutations, this
//' processes approximately 20 million permutation statistics.
//'
//' @export
// [[Rcpp::export]]
DataFrame batch_permutation_test_cpp(
    const arma::vec& z_f,
    const arma::mat& Z_g,
    const arma::sp_mat& W,
    const std::string& metric,
    int n_perm,
    double weight_sum = -1.0
) {
    int n_genes = Z_g.n_cols;
    int n_obs = z_f.n_elem;

    // Use n_obs for normalization if weight_sum not provided
    double normalizer = (weight_sum > 0) ? weight_sum : static_cast<double>(n_obs);

    // =========================================================================
    // STEP 1: Pre-generate ALL permutation indices (key optimization!)
    // This ensures consistency across genes and avoids repeated RNG calls
    // =========================================================================
    arma::umat perm_indices(n_perm, n_obs);
    for (int b = 0; b < n_perm; b++) {
        perm_indices.row(b) = arma::randperm(n_obs).t();
    }

    // Pre-compute factor norm (for I_ND)
    double norm_f = arma::norm(z_f, 2);

    // Result vectors
    NumericVector out_I_obs(n_genes);
    NumericVector out_p_value(n_genes);
    NumericVector out_z_score(n_genes);
    NumericVector out_null_sd(n_genes);

    // =========================================================================
    // STEP 2: Process each gene
    // =========================================================================
    for (int g = 0; g < n_genes; g++) {
        // Check for user interrupt every 100 genes
        if (g % 100 == 0) {
            Rcpp::checkUserInterrupt();
        }

        // Compute spatial lag for this gene
        arma::vec lag_g = W * Z_g.col(g);
        double norm_lag = arma::norm(lag_g, 2);

        // Check for degenerate cases
        if (norm_lag < 1e-10 || (metric == "ind" && norm_f < 1e-10)) {
            out_I_obs[g] = NA_REAL;
            out_p_value[g] = NA_REAL;
            out_z_score[g] = NA_REAL;
            out_null_sd[g] = NA_REAL;
            continue;
        }

        // Compute observed statistic
        double I_obs;
        if (metric == "moran") {
            I_obs = arma::dot(z_f, lag_g) / normalizer;
        } else {
            I_obs = arma::dot(z_f, lag_g) / (norm_f * norm_lag);
        }

        // Permutation test using pre-generated indices
        int count_extreme = 0;
        double sum_perm = 0.0;
        double sum_perm_sq = 0.0;

        for (int b = 0; b < n_perm; b++) {
            // Get permuted factor using pre-computed indices
            arma::vec z_f_perm = z_f.elem(perm_indices.row(b).t());

            // Compute permuted statistic
            double I_perm;
            if (metric == "moran") {
                I_perm = arma::dot(z_f_perm, lag_g) / normalizer;
            } else {
                I_perm = arma::dot(z_f_perm, lag_g) / (norm_f * norm_lag);
            }

            // Accumulate for null distribution
            sum_perm += I_perm;
            sum_perm_sq += I_perm * I_perm;

            // Two-sided test
            if (std::abs(I_perm) >= std::abs(I_obs)) {
                count_extreme++;
            }
        }

        // P-value with continuity correction
        double p_value = (1.0 + count_extreme) / (n_perm + 1.0);

        // Z-score
        double null_mean = sum_perm / n_perm;
        double null_var = sum_perm_sq / n_perm - null_mean * null_mean;
        double null_sd = std::sqrt(std::max(null_var, 1e-10));

        out_I_obs[g] = I_obs;
        out_p_value[g] = p_value;
        out_z_score[g] = (I_obs - null_mean) / null_sd;
        out_null_sd[g] = null_sd;
    }

    return DataFrame::create(
        Named("gene_idx") = seq(1, n_genes),
        Named("I_obs") = out_I_obs,
        Named("p_value") = out_p_value,
        Named("z_score") = out_z_score,
        Named("null_sd") = out_null_sd
    );
}


// =============================================================================
// ALL-PAIRS PERMUTATION TESTING (Pairwise Moran's I Matrix)
// =============================================================================

//' All-Pairs Permutation Test for Pairwise Moran's I
//'
//' Fast C++ implementation of permutation testing for the full pairwise
//' Moran's I matrix. This function is called by the R wrapper moran_permutation_test().
//'
//' @param data_z Z-normalized expression matrix (genes x spots)
//' @param W Dense weight matrix (spots x spots)
//' @param S0 Sum of all weights
//' @param n_perm Number of permutations
//' @param seed Random seed for reproducibility (0 = use random seed)
//' @param verbose Print progress messages
//'
//' @return List containing:
//'   - sum: Sum of Moran's I across permutations
//'   - sum_sq: Sum of squared Moran's I
//'   - count_extreme: Count of permutations with |I_perm| >= |I_obs|
//'   - observed: Observed Moran's I matrix
//'
//' @keywords internal
// [[Rcpp::export]]
List allpairs_permutation_test_cpp(
    const arma::mat& data_z,
    const arma::mat& W,
    double S0,
    int n_perm,
    unsigned int seed = 0,
    bool verbose = true
) {
    int n_genes = data_z.n_rows;
    int n_spots = data_z.n_cols;

    // Compute observed Moran's I: M = Z * W * Z^T / S0
    arma::mat observed = (data_z * W * data_z.t()) / S0;
    arma::mat abs_observed = arma::abs(observed);

    // Initialize accumulators
    arma::mat perm_sum(n_genes, n_genes, arma::fill::zeros);
    arma::mat perm_sum_sq(n_genes, n_genes, arma::fill::zeros);
    arma::mat count_extreme(n_genes, n_genes, arma::fill::zeros);

    // Random number generator
    std::mt19937 gen;
    if (seed == 0) {
        std::random_device rd;
        gen.seed(rd());
    } else {
        gen.seed(seed);
    }

    // Create permutation indices
    std::vector<int> indices(n_spots);
    std::iota(indices.begin(), indices.end(), 0);

    // Progress reporting
    int report_interval = std::max(1, n_perm / 10);

    // Temporary matrix for permuted data
    arma::mat data_perm(n_genes, n_spots);

    for (int p = 0; p < n_perm; p++) {
        // Progress reporting
        if (verbose && ((p + 1) % report_interval == 0)) {
            Rcpp::Rcout << "  Permutation " << (p + 1)
                        << " / " << n_perm << std::endl;
        }

        // Check for user interrupt
        if ((p + 1) % 100 == 0) {
            Rcpp::checkUserInterrupt();
        }

        // Shuffle indices
        std::shuffle(indices.begin(), indices.end(), gen);

        // Permute data columns
        for (int j = 0; j < n_spots; j++) {
            data_perm.col(j) = data_z.col(indices[j]);
        }

        // Compute permuted Moran's I
        arma::mat moran_perm = (data_perm * W * data_perm.t()) / S0;

        // Update running sum and sum of squares
        perm_sum += moran_perm;
        perm_sum_sq += moran_perm % moran_perm;  // element-wise square

        // Count extreme values (two-sided test)
        arma::mat abs_perm = arma::abs(moran_perm);
        count_extreme += arma::conv_to<arma::mat>::from(abs_perm >= abs_observed);
    }

    return List::create(
        Named("sum") = perm_sum,
        Named("sum_sq") = perm_sum_sq,
        Named("count_extreme") = count_extreme,
        Named("observed") = observed
    );
}


//' Fast Pairwise Moran's I Computation
//'
//' Computes pairwise Moran's I for all gene pairs using optimized BLAS
//' matrix operations.
//'
//' @param data_z Z-normalized expression matrix (genes x spots)
//' @param W Weight matrix (spots x spots)
//' @param S0 Sum of all weights
//'
//' @return Moran's I matrix (genes x genes)
//'
//' @keywords internal
// [[Rcpp::export]]
arma::mat pairwise_moran_matrix_cpp(
    const arma::mat& data_z,
    const arma::mat& W,
    double S0
) {
    // M = Z * W * Z^T / S0
    arma::mat ZW = data_z * W;
    arma::mat M = ZW * data_z.t() / S0;

    return M;
}


//' Z-Normalize Expression Matrix
//'
//' Z-normalizes each row (gene) of the expression matrix.
//'
//' @param X Expression matrix (genes x spots)
//'
//' @return Z-normalized matrix
//'
//' @keywords internal
// [[Rcpp::export]]
arma::mat z_normalize_matrix_cpp(const arma::mat& X) {
    int n_genes = X.n_rows;
    int n_spots = X.n_cols;

    arma::mat Z(n_genes, n_spots);

    for (int g = 0; g < n_genes; g++) {
        arma::rowvec row = X.row(g);
        double mu = arma::mean(row);
        double sd = arma::stddev(row, 0);  // 0 = normalize by n-1

        if (sd < 1e-10) {
            // Constant expression - set to zeros
            Z.row(g).zeros();
        } else {
            Z.row(g) = (row - mu) / sd;
        }
    }

    return Z;
}


//' Univariate Moran's I for Each Gene
//'
//' Computes univariate Moran's I for each gene (diagonal of pairwise matrix).
//'
//' @param data_z Z-normalized expression matrix (genes x spots)
//' @param W Weight matrix (spots x spots)
//' @param S0 Sum of all weights
//'
//' @return Vector of Moran's I values
//'
//' @keywords internal
// [[Rcpp::export]]
arma::vec single_gene_moran_cpp(
    const arma::mat& data_z,
    const arma::mat& W,
    double S0
) {
    int n_genes = data_z.n_rows;
    arma::vec moran_i(n_genes);

    for (int g = 0; g < n_genes; g++) {
        arma::rowvec z = data_z.row(g);
        // I = z * W * z^T / S0
        double I = arma::as_scalar(z * W * z.t()) / S0;
        moran_i(g) = I;
    }

    return moran_i;
}
