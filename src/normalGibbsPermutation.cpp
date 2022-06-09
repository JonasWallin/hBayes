#include <RcppArmadillo.h>
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace Rcpp;
using namespace arma;


double normalloglik(const arma::vec& res, const double sigma)
{
    double lik = - ((double) res.size() )* log(sigma);
    lik += -0.5 *  sum(res % res)/pow(sigma,2);
    return lik;
}


// [[Rcpp::export]]
Rcpp::List permute_gibbs_beta_normal_cpp(const arma::vec& Y,
                                        arma::vec& Xbeta,
                                        const arma::mat& X,
                                        const double sigma,
                                        arma::vec& beta,
                                             arma::ivec& betaind) {

    betaind -= 1;


    const arma::uword p = X.n_cols;

    arma::vec beta_new(p, fill::zeros);
    arma::vec acc_vec(p, fill::zeros);
    arma::vec residual  = Y - Xbeta;

    double log_lik        =   normalloglik(residual, sigma);

    for (uword i = 0; i < p; ++i) {

        arma::vec  Xi = X.col(i);
        double U_  = R::runif(-0.5+1e-8, p-0.5+1e-8);

        //
        int   proposal = (int) std::round(U_);

        arma::vec  Xj = X.col(proposal);
        //log_lik_star
        arma::vec  Xbeta_star =  Xbeta + Xi * (beta(proposal) - beta(i)) - Xj * (beta(proposal) - beta(i));

        arma::vec residual_star     = Y - Xbeta_star;
        double log_lik_star     =   normalloglik(residual_star, sigma);
        double log_U = std::log(R::runif(0,1));
        double MH_ratio = log_lik_star - log_lik ;
        if(log_U  < MH_ratio){
            double beta_swap_star =beta(proposal);
            int beta_ind_prop = betaind(proposal);
            beta(proposal)  = beta(i);

            betaind(proposal) = betaind(i);
            betaind(i) = beta_ind_prop;
            beta(i)     = beta_swap_star;
            log_lik     = log_lik_star;
            residual    = residual_star;
            Xbeta = Xbeta_star;
            acc_vec(i) += 1.;
        }
    }


    betaind += 1;
    return List::create(Named("acc.vec") = wrap(acc_vec),
                        Named("beta.ind") = wrap(betaind));
}
