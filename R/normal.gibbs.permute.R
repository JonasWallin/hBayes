


#' Gibbs sampler derived from Yekuteili code
#' then in the gibbs sampler we randomly permute the beta coeffients (one at the
#' time)
#' where \pi_q is a density estimated from data
#'
#' @paran n.gibbs - (int ) number of gibbs samples to be performed
#' @param Y       -  ( n x 1) observastions
#' @param X       -  ( n x p)  covariates vector
#' @param sigma   -  (double) standard devation
#' @param a.in    -  (2 x 1) lower and upper point of the density
#' @param L       -  (int)   2^L is number of bins in the tree
#' @param beta    -  (p x 1)  inital guess of beta
#' @param cpp     -  (bool) use C++
#'
gibbs.normal.permute.beta.fixed.sigma          <- function(n.gibbs,
                                                               Y,
                                                               X,
                                                               sigma,
                                                               a.int,
                                                               L,
                                                               beta,
                                                               quite = T){



    p			<- dim(X)[2]
    n			<- dim(X)[1]
    tree.data <- build.tree(L)

    I = tree.data$I


    a.vec	<- seq(a.int[1],
                 a.int[2],
                 length = I + 1)

    #truncate beta
    beta[beta< a.int[1]] = a.int[1] + 1e-8
    beta[beta> a.int[2]] = a.int[2] - 1e-8


    pi.gibbs		<- array(dim = c(n.gibbs, I))
    beta.gibbs	    <- array(dim = c(n.gibbs, p))
    delta.gibbs    <-  array(dim = c(n.gibbs, p))



    beta.gibbs[1,]	<- beta
    beta.ind <- findInterval(beta, a.vec, left.open = T)

    pi.gibbs[1,]	<- log( (c(sum(beta.gibbs[1,] < a.vec[2]),diff(rank(c(a.vec[2:I],100,beta.gibbs[1,]))[1:I] -
                                                                    (1:I))) + 100/I) / 900)
    delta.gibbs[1, ] <- beta.ind
    Xbeta	<- X %*% beta.gibbs[1,]




    for (g in 2:n.gibbs)
    {
        if(quite==F){
            if (g %% 100 == 1)	print(paste(g))
        }


        ##
        # sample beta
        ##
        pi.gibbs_g       <-  pi.gibbs[g-1,]
        beta.gibbs.g     <- beta.gibbs[g-1,]
        res.beta <- permute_gibbs_beta_normal_cpp(Y,
                                           Xbeta,
                                           X,
                                         sigma,
                                          beta.gibbs.g,
                                         beta.ind)

        beta.ind <-         res.beta$beta.ind
        beta.gibbs[g,]   <- beta.gibbs.g

        delta.gibbs[g, ] <- beta.ind


        ##
        # sample pi
        ##
        nvec <- table(factor(beta.ind, levels=1:I))
        pi.gibbs[g,] <- polya.gibs(nvec, tree.data)
    }
    return(list(a.vec=  a.vec,
                pi.gibbs = pi.gibbs,
                beta.gibbs = beta.gibbs,
                delta.gibbs = delta.gibbs))

}



