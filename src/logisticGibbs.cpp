#include <RcppArmadillo.h>
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace Rcpp;
using namespace arma;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//


double binomaloglik(const arma::vec& Y, const arma::vec& Xbeta, const arma::vec& expXbeta)
{
    return sum(Y % Xbeta - log1p(expXbeta));
}

// [[Rcpp::export]]
Rcpp::List sample_gibbs_beta_logistic_cpp(const arma::vec& Y,
                                 arma::vec& Xbeta,
                                 const arma::mat& X,
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
    arma::vec expXbeta    = exp(Xbeta);
    arma::vec prob        = expXbeta/ (1 + expXbeta );
    arma::vec w           = prob % (1-prob);
    arma::vec residual    = Y - prob;
    double log_lik        =  binomaloglik(Y, Xbeta, expXbeta);

    for (uword i = 0; i < p; ++i) {

        arma::vec  Xi = X.col(i);
        arma::vec Xi2 = Xi % Xi;
        for(uword ii = 0; ii < innersample; ii++){
            int curpos = betaind(i);
            int newpos = curpos;
            double U_  = R::runif(-interval_sample-0.5,interval_sample+0.5);
            newpos += (int) std::round(U_);


            if(newpos < 0 || newpos >= m1)
                continue;

            double H = sum(w % Xi2) + 1e-8;
            double mu = beta(i)+sum(Xi % residual)/H;
            double sd = std::sqrt(1/H);
            double upp = avec(newpos+1);
            double low = avec(newpos);
            double Phi_upp = R::pnorm(upp, mu, sd, 1, 0);
            double Phi_low = R::pnorm(low, mu, sd, 1, 0);
            double Z =Phi_upp - Phi_low;
            double U_star  = R::runif(0,1);
            double beta_i_star = R::qnorm(Phi_low + U_star * Z, 0., 1., 1, 0 )*sd + mu;


            double log_qxy =  R::dnorm(beta_i_star, mu, sd, 1) - std::log(Z);
            arma::vec  Xbeta_star =  Xbeta + Xi * (beta_i_star- beta(i));

            arma::vec expXbeta_star = exp(Xbeta_star);
            double log_lik_star     =   binomaloglik(Y, Xbeta_star, expXbeta_star);


            arma::vec prob_star         = expXbeta_star/(1 + expXbeta_star );
            arma::vec    w_star         = prob_star % (1-prob_star);
            arma::vec residual_star     = Y - prob_star;
            double H_star = sum(w_star % Xi2) + 1e-8;
            double mu_star = beta_i_star + sum(Xi % residual_star)/H_star;
            double sd_star = std::sqrt(1/H_star);


            upp = avec(curpos+1);
            low = avec(curpos);
            Phi_upp = R::pnorm(upp, mu_star, sd_star, 1, 0);
            Phi_low = R::pnorm(low, mu_star, sd_star, 1, 0);
            Z       = Phi_upp - Phi_low;
            double log_qyx =  R::dnorm(beta(i), mu_star, sd_star, 1) - std::log(Z);

            double log_U = std::log(R::runif(0,1));
            double MH_ratio = ( log_lik_star + pivec(newpos) )- (log_lik + pivec(curpos))  + log_qyx - log_qxy;


            if(log_U  < MH_ratio){
                beta(i)     = beta_i_star;
                betaind(i) = newpos;
                log_lik     = log_lik_star;
                w           = w_star;
                residual    = residual_star;
                Xbeta       = Xbeta_star;
                acc_vec(i) += 1.;
            }


        }
    }


    acc_vec /= innersample;
    betaind += 1;
    return List::create(Named("acc.vec") = wrap(acc_vec),
                        Named("beta.ind") = wrap(betaind));
}

