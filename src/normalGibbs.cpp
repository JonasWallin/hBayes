#include <RcppArmadillo.h>
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace Rcpp;
using namespace arma;


double normalloglik_(const arma::vec& res, const double sigma)
{
    double lik = - ((double) res.size() )* log(sigma);
    lik += -0.5 *  sum(res % res)/pow(sigma,2);
    return lik;
}


// [[Rcpp::export]]
Rcpp::List sample_gibbs_beta_normal_cpp(const arma::vec& Y,
                                          arma::vec& Xbeta,
                                          const arma::mat& X,
                                          const double sigma,
                                          const arma::vec& X2diag,
                                          arma::vec& beta,
                                          arma::ivec& betaind,
                                          const arma::vec& avec,
                                          const arma::vec& pivec,
                                          int innersample = 5,
                                          double interval_sample = 2.) {


    betaind -= 1;
    const arma::uword n = X.n_rows;
    const arma::uword p = X.n_cols;

    arma::vec beta_new(p, fill::zeros);
    const arma::uword m1 = pivec.size();
    arma::vec acc_vec(p, fill::zeros);
    arma::vec residual  = Y - Xbeta;

    arma::vec H  =  X2diag;
    arma::vec sd  = pow(sigma,2)/H;
    sd = arma::sqrt(sd);
    double log_lik        =   normalloglik_(residual, sigma);

    for (uword i = 0; i < p; ++i) {

        arma::vec  Xi = X.col(i);
        for(uword ii = 0; ii < innersample; ii++){
            int curpos = betaind(i);
            int newpos = curpos;
            double U_  = R::runif(-interval_sample-0.5,interval_sample+0.5);
            newpos += (int) std::round(U_);


            if(newpos < 0 || newpos >= m1)
                continue;

            double mu = beta(i)+sum(Xi % residual)/H(i);


            double upp = avec(newpos+1);
            double low = avec(newpos);
            double Phi_upp = R::pnorm(upp, mu, sd(i), 1, 0);
            double Phi_low = R::pnorm(low, mu, sd(i), 1, 0);
            double Z =Phi_upp - Phi_low;

            double U_star;
            double beta_i_star;
            double log_qxy;
            if(Phi_low > 1-1e-10){
                beta_i_star  = low;
                log_qxy      = 0;

            }else if(Phi_upp < 1e-10){
                beta_i_star = upp;
                log_qxy    = 0;

            }else{
                U_star  = R::runif(0,1);
                beta_i_star = R::qnorm(Phi_low + U_star * Z, 0., 1., 1, 0 )*sd(i) + mu;
                log_qxy =  R::dnorm(beta_i_star, mu, sd(i), 1) - std::log(Z);
            }
            arma::vec  Xbeta_star =  Xbeta + Xi * (beta_i_star- beta(i));

            arma::vec residual_star     = Y - Xbeta_star;
            double log_lik_star     =   normalloglik_(residual_star, sigma);


            double mu_star = beta_i_star + sum(Xi % residual_star)/H(i);


            upp = avec(curpos+1);
            low = avec(curpos);
            Phi_upp = R::pnorm(upp, mu_star, sd(i), 1, 0);
            Phi_low = R::pnorm(low, mu_star, sd(i), 1, 0);
            Z       = Phi_upp - Phi_low;
            double log_qyx =  R::dnorm(beta(i), mu_star, sd(i), 1) - std::log(Z);

            double log_U = std::log(R::runif(0,1));
            double MH_ratio = ( log_lik_star + pivec(newpos) )- (log_lik + pivec(curpos))  + log_qyx - log_qxy;
            if(log_U  < MH_ratio){
                beta(i)     = beta_i_star;
                betaind(i) = newpos;
                log_lik     = log_lik_star;
                residual    = residual_star;
                Xbeta = Xbeta_star;
                acc_vec(i) += 1.;
            }


        }
    }


    acc_vec /= innersample;
    betaind += 1;
    return List::create(Named("acc.vec") = wrap(acc_vec),
                        Named("beta.ind") = wrap(betaind));
}



