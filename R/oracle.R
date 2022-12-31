

#' Oracle Gibbs sampler derived
#' then in the gibbs sampler we randomly permute the beta coeffients (one at the
#' time
#'
#' @paran n.gibbs - (int ) number of gibbs samples to be performed
#' @param Y       -  ( n x 1) observastions
#' @param X       -  ( n x p)  covariates vector
#' @param sigma   -  (double) standard devation
#' @param beta    -  (p x 1)  the true beta vector (or a random permutation of it)
#'
gibbs.logistic.permute.oracle          <- function(n.gibbs,
                                                   Y,
                                                   X,
                                                   beta,
                                                   quite = T){


    p		     	<- dim(X)[2]
    beta.gibbs	    <- array(dim = c(n.gibbs, p))
    beta.gibbs[1,]	<- beta
    Xbeta	<- X %*% beta.gibbs[1,]

    for (g in 2:n.gibbs)
    {
        if(quite==F){
            if (g %% 100 == 1)	print(paste(g))
        }


        ##
        # sample beta
        ##
        beta.gibbs.g     <- beta.gibbs[g-1,]
        res.beta <- permute_gibbs_beta_logistic_cpp(Y,
                                                  Xbeta,
                                                  X,
                                                  beta.gibbs.g)

        beta.gibbs[g,]   <- beta.gibbs.g



    }
    return(beta.gibbs)

}
