#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
void house_reflection_mult_inplace_cpp(arma::mat& Gamma,
                                    arma::vec v_p,
                                    bool right = true) {
    int p = Gamma.n_rows;
    int m = v_p.n_elem;

    // Adjusting v_p[1]
    v_p[0] = v_p[0] - 1;

    arma::vec x_m = - sqrt(2.0) * v_p / arma::norm(v_p, 2);

    arma::uvec ind = arma::regspace<arma::uvec>(p - m , p-1);

    if (right) {
        Gamma.rows(ind) -= x_m * (x_m.t() * Gamma.rows(ind));
    } else {
        Gamma.cols(ind) -= (Gamma.cols(ind) * x_m) * (x_m.t());
    }
}

// [[Rcpp::export]]
arma::mat house_reflection_mult_cpp(arma::mat Gamma,
                                arma::vec v_p,
                                bool right = true) {
    int p = Gamma.n_rows;
    int m = v_p.n_elem;

    // Adjusting v_p[1]
    v_p[0] = v_p[0] - 1;

    arma::vec x_m = - sqrt(2.0) * v_p / arma::norm(v_p, 2);

    arma::uvec ind = arma::regspace<arma::uvec>(p - m , p-1);

    if (right) {
        Gamma.rows(ind) -= x_m * (x_m.t() * Gamma.rows(ind));
    } else {
        Gamma.cols(ind) -= (Gamma.cols(ind) * x_m) * (x_m.t());
    }

    return Gamma;
}
/*** R
# Example usage in R
set.seed(42)
Gamma <- matrix(rnorm(12), nrow = 4)
v_p <- c(1, 2, 3)
house_reflection_mult(Gamma, v_p, right = TRUE)
*/
