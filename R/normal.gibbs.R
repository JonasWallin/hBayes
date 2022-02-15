


 normalloglik<- function( res, sigma)
{
    lik = - length(res) * log(sigma) - 0.5 * sum(res * res)/sigma^2;
    return(lik)
}

#' single loop through the gibbs sampler
#'
#' @param Y        (n x 1) observation
#' @param Xbeta    (n x 1) X%*%beta
#' @param X        (n x p) covariets
 #' @param X2diag        (p x 1) inner product of each covariate
#' @param beta     (p x 1) current sample
#' @param beta.ind (p x 1) current location of beta in a.vec
#' @param a.vec    (m x 1) domain of pi
#' @param pi.vec   (m-1 x 1) density value of a.vec (log)
#'
sample.gibbs.normal.beta <- function(Y,
                                     Xbeta,
                                     X,
                                     sigma,
                                     X2diag,
                                     beta,
                                     beta.ind,
                                     a.vec,
                                     pi.vec,
                                     inner.sample = 5){

    d <- length(beta)
    beta.new <- beta
    m.1 <- length(pi.vec)
    acc.vec <- rep(0, d)

    residual    = Y - Xbeta
    lik.log <-  normalloglik(residual, sigma)
    sd = sigma*sqrt(1/X2diag)
    for(i in 1:d){
        X.i  = X[,i]
        for(ii in 1:inner.sample){
            #sample position move
            cur.pos <- beta.ind[i]
            new.pos <- beta.ind[i] + sample(-5:5,1)#round(3.*runif(1)-1.5)
            if(new.pos < 1 || new.pos > m.1 )
                next


            mu = beta[i]+sum(X.i * residual)/X2diag[i]

            # sample constraint
            upp = a.vec[new.pos+1]
            low = a.vec[new.pos]
            Phi.upp = pnorm((upp - mu)/sd[i])
            Phi.low = pnorm((low - mu)/sd[i])
            Z =Phi.upp - Phi.low
            U  = runif(1)
            if(Phi.low > 1-1e-10){
                beta.i.star <- low
                qxy.log      <- 0

            }else if(Phi.upp<1e-10){
                beta.i.star <- upp
                qxy.log      <- 0

            }else{
                beta.i.star <- qnorm(Phi.low + U * Z)*sd[i] + mu
                qxy.log      <-  dnorm(beta.i.star, mu, sd[i], log = T) - log(Z)
            }
            Xbeta.star =  Xbeta + X.i*(beta.i.star- beta[i])
            residual.star        = Y - Xbeta.star

            lik.star.log <-   normalloglik(residual.star, sigma)

            # compute the opposite proposal
            mu.star = beta.i.star + sum(X.i * residual.star)/X2diag[i]


            upp = a.vec[cur.pos+1]
            low = a.vec[cur.pos]
            Phi.upp = pnorm((upp - mu.star)/sd[i])
            Phi.low = pnorm((low - mu.star)/sd[i])
            Z       = Phi.upp - Phi.low
            qyx.log      <-  dnorm(beta[i], mu.star, sd[i], log = T) - log(Z)

            if(log(runif(1)) < lik.star.log + pi.vec[new.pos ] - pi.vec[cur.pos ]-  lik.log + qyx.log - qxy.log){
                beta[i]     = beta.i.star
                beta.ind[i] = new.pos
                lik.log     = lik.star.log
                residual    = residual.star
                Xbeta       = Xbeta.star
                acc.vec[i] <- acc.vec[i]+  1
            }
        }
    }
    acc.vec <- acc.vec/inner.sample
    return(list( Xbeta = Xbeta,
                 beta = beta,
                 beta.ind = beta.ind,
                 acc.vec = acc.vec))

}



# gibbs sampler derived from Yekuteili code
#'
#' @paran n.gibbs - (int ) number of gibbs samples to be performed
#' @param Y       -  ( n x 1) observastions
#' @param X       -  ( n x p)  covariates vector
#' @param a.in    -  (2 x 1) lower and upper point of the density
#' @param L       -  (int)   2^L is number of bins in the tree
#' @param beta   -  (p x 1)  inital guess of beta
#'
gibbs.normal           <- function(n.gibbs,
                                     Y,
                                     X,
                                     a.int,
                                     L,
                                     beta,
                                     cpp = T){



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


    X2diag <- rep(0, p)
    for(i in 1:p)
        X2diag[i] <- sum(X[,i]^2)

    ##
    # inital sample simga
    ##
    beta.sigma = 0.5 * sum((Y - Xbeta)^2)
    alpha.sigma = 0.5 * n + 1
    sigma2 = 1/rgamma(1, shape = alpha.sigma, rate= beta.sigma)
    sigma = sqrt(sigma2)
    sigma.gibbs <- rep(0, n.gibbs)
    sigma.gibbs[1] <- sigma



    for (g in 2:n.gibbs)
    {
        if (g %% 100 == 1)	print(paste(g))


        ##
        # sample sigma (Jeffers prior)
        ##
        beta.sigma = 0.5 * sum((Y - Xbeta)^2)
        alpha.sigma = 0.5 * n + 1
        sigma2 = 1/rgamma(1, shape = alpha.sigma, rate= beta.sigma)
        sigma = sqrt(sigma2)

        ##
        # sample beta
        ##
        pi.gibbs_g       <-  pi.gibbs[g-1,]
        beta.gibbs.g     <- beta.gibbs[g-1,]
        if(cpp){
            res.beta    <- sample_gibbs_beta_normal_cpp(Y,
                                                         Xbeta,
                                                         X,
                                                         sigma,
                                                         X2diag,
                                                         beta.gibbs.g,
                                                         beta.ind,
                                                         a.vec,
                                                         pi.gibbs_g,
                                                         innersample= 10,
                                                         interval_sample=3)
            beta.ind <-         res.beta$beta.ind
            beta.gibbs[g,]   <- beta.gibbs.g
        }else{
            res.beta    <- sample.gibbs.normal.beta(Y,
                                                         Xbeta,
                                                         X,
                                                         sigma,
                                                         X2diag,
                                                         beta.gibbs.g,
                                                         beta.ind,
                                                         a.vec,
                                                         pi.gibbs_g)
            Xbeta = res.beta$Xbeta
            beta.gibbs[g,]   <- res.beta$beta
            beta.ind <-         res.beta$beta.ind


        }
        delta.gibbs[g, ] <- beta.ind


        ##
        # sample sigma (Jeffers prior)
        ##
        beta.sigma = 0.5 * sum((Y - Xbeta)^2)
        alpha.sigma = 0.5 * n + 1
        sigma2 = 1/rgamma(1, shape = alpha.sigma, rate= beta.sigma)
        sigma = sqrt(sigma2)
        sigma.gibbs[g] <- sigma
        ##
        # sample pi
        ##
        nvec <- table(factor(beta.ind, levels=1:I))
        pi.gibbs[g,] <- polya.gibs(nvec, tree.data)
    }
    return(list(a.vec=  a.vec,
                pi.gibbs = pi.gibbs,
                beta.gibbs = beta.gibbs,
                delta.gibbs = delta.gibbs,
                sigma.gibbs    = sigma.gibbs))

}



#' Gibbs sampler derived from Yekuteili code
#' Here we sample
#' \beta from the density
#'  \pi(\beta|X, y, \sigma) \sim N(y; X\beta, \sigma^2) \pi_q(\beta)
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
gibbs.normal.fixed.sigma           <- function(n.gibbs,
                                               Y,
                                               X,
                                               sigma,
                                               a.int,
                                               L,
                                               beta,
                                               cpp = T){



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


    X2diag <- rep(0, p)
    for(i in 1:p)
        X2diag[i] <- sum(X[,i]^2)



    for (g in 2:n.gibbs)
    {
        if (g %% 100 == 1)	print(paste(g))


        ##
        # sample beta
        ##
        pi.gibbs_g       <-  pi.gibbs[g-1,]
        beta.gibbs.g     <- beta.gibbs[g-1,]
        if(cpp){
            res.beta    <- sample_gibbs_beta_normal_cpp(Y,
                                                        Xbeta,
                                                        X,
                                                        sigma,
                                                        X2diag,
                                                        beta.gibbs.g,
                                                        beta.ind,
                                                        a.vec,
                                                        pi.gibbs_g,
                                                        innersample= 10,
                                                        interval_sample=3)
            beta.ind <-         res.beta$beta.ind
            beta.gibbs[g,]   <- beta.gibbs.g
        }else{
            res.beta    <- sample.gibbs.normal.beta(Y,
                                                    Xbeta,
                                                    X,
                                                    sigma,
                                                    X2diag,
                                                    beta.gibbs.g,
                                                    beta.ind,
                                                    a.vec,
                                                    pi.gibbs_g)
            Xbeta = res.beta$Xbeta
            beta.gibbs[g,]   <- res.beta$beta
            beta.ind <-         res.beta$beta.ind


        }
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

