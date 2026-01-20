#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Row-Normalize Sparse Matrix
//'
//' @param W Sparse matrix (n x m)
//' @return Row-normalized sparse matrix where each row sums to 1 (or 0)
// [[Rcpp::export]]
arma::sp_mat sparse_row_normalize_cpp(const arma::sp_mat& W) {
    arma::sp_mat W_norm = W;

    for (arma::uword i = 0; i < W.n_rows; i++) {
        double row_sum = arma::accu(W.row(i));
        if (row_sum > 1e-10) {
            // Iterate over non-zero elements in row
            for (arma::sp_mat::const_row_iterator it = W.begin_row(i);
                 it != W.end_row(i); ++it) {
                W_norm(i, it.col()) = (*it) / row_sum;
            }
        }
    }

    return W_norm;
}
