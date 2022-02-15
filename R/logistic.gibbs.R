###
# simple single gibbs sampler
#
#
###



#' single loop through the gibbs sampler
#'
#' @param Y        (n x 1) observation
#' @param Xbeta    (n x 1) X%*%beta
#' @param X        (n x p) covariets
#' @param beta     (p x 1) current sample
#' @param beta.ind (p x 1) current location of beta in a.vec
#' @param a.vec    (m x 1) domain of pi
#' @param pi.vec   (m-1 x 1) density value of a.vec (log)
#'
sample.gibbs.beta <- function(Y, Xbeta, X, beta, beta.ind, a.vec, pi.vec, inner.sample = 5){

    d <- length(beta)
    beta.new <- beta
    m.1 <- length(pi.vec)
    acc.vec <- rep(0, d)
    expXbeta    = exp(Xbeta)
    prob        = expXbeta/(1 + expXbeta )
    w           = prob * (1-prob)
    residual    = Y - prob
    lik.log <-  biomial.loglik(Y, Xbeta, expXbeta)
    for(i in 1:d){
        X.i  = X[,i]
        X.i2 = X.i^2
        for(ii in 1:inner.sample){
            #sample position move
            cur.pos <- beta.ind[i]
            new.pos <- beta.ind[i] + sample(-5:5,1)#round(3.*runif(1)-1.5)
            if(new.pos < 1 || new.pos > m.1 )
                next


            H = sum(w*X.i2) + 10^-8
            mu = beta[i]+sum(X.i * residual)/H
            sd = sqrt(1/H)
            # sample constraint
            upp = a.vec[new.pos+1]
            low = a.vec[new.pos]
            Phi.upp = pnorm((upp - mu)/sd)
            Phi.low = pnorm((low - mu)/sd)
            Z =Phi.upp - Phi.low
            U  = runif(1)
            beta.i.star <- qnorm(Phi.low + U * Z)*sd + mu

            #
            qxy.log      <-  dnorm(beta.i.star, mu, sd, log = T) - log(Z)
            Xbeta.star =  Xbeta + X.i*(beta.i.star- beta[i])
            expXbeta.star = exp(Xbeta.star)
            lik.star.log <-   biomial.loglik(Y, Xbeta.star, expXbeta.star)

            # compute the opposite proposal
            prob.star            = expXbeta.star/(1 + expXbeta.star )
            w.star               = prob.star* (1-prob.star)
            residual.star        = Y - prob.star
            H.star = sum(w.star*X.i2) + 10^-8
            mu.star = beta.i.star + sum(X.i * residual.star)/H.star
            sd.star = sqrt(1/H.star)

            upp = a.vec[cur.pos+1]
            low = a.vec[cur.pos]
            Phi.upp = pnorm((upp - mu.star)/sd.star)
            Phi.low = pnorm((low - mu.star)/sd.star)
            Z       = Phi.upp - Phi.low
            qyx.log      <-  dnorm(beta[i], mu.star, sd.star, log = T) - log(Z)
            if(log(runif(1)) < lik.star.log + pi.vec[new.pos] - pi.vec[cur.pos]-  lik.log + qyx.log - qxy.log){
                beta[i]     = beta.i.star
                beta.ind[i] = new.pos
                lik.log     = lik.star.log
                w           = w.star
                residual    = residual.star
                Xbeta       = Xbeta.star
                acc.vec[i] <- acc.vec[i]+  1
            }
            }
    }
    acc.vec <- acc.vec/inner.sample
    return(list( Xbeta = Xbeta, beta = beta, beta.ind = beta.ind,acc.vec = acc.vec))

}

#' single loop through the gibbs sampler (Yekutieli version)
#' using a linear approximation of the density
#'
#'
sample.gibbs.beta.Yekutieli <- function(num.y.1,
                                        ind.y.1,
                                        beta.vals,
                                  I,
                                 Xbeta,
                                 X,
                                 beta,
                                 a.vec,
                                 pi.log,
                                 beta.vals.a.ind){


    n <- dim(X)[1]
    p <- dim(X)[2]
    nvec <- rep(0,I)
    Vec1 <- rbind(rep(1, num.y.1))
    Vec2 <- rbind(rep(1, n))
    delta.gbbs <- rep(0,length(beta))
    for (k in 1:p)
    {
        Xbeta			<- Xbeta - X[,k]*beta[k]
        mult.beta.k.lin.pred	<- rep(Xbeta, I + 1) + rep(X[,k], I + 1)*rep(a.vec, each = n)

        mult.beta.k.loglik.num	 <- Vec1 %*% array(mult.beta.k.lin.pred[ind.y.1], dim = c(num.y.1,I + 1))
        mult.beta.k.loglik.denom <- Vec2 %*% array(log(1 + exp(mult.beta.k.lin.pred)), dim = c(n, I + 1))

        # standradize density
        mult.beta.k.loglik       <- mult.beta.k.loglik.num - mult.beta.k.loglik.denom
        mult.beta.k.loglik       <- mult.beta.k.loglik - max(mult.beta.k.loglik)
        # linear approximation
        beta.vals.loglik <- approx(a.vec,mult.beta.k.loglik,xout = beta.vals)$y
        # add density
        beta.pst <- exp(beta.vals.loglik + pi.log)
        beta.ind <- sample(1:length(beta.pst),1,prob = beta.pst / sum(beta.pst))

        delta.gbbs[ k] <- beta.vals.a.ind[beta.ind]
        nvec[delta.gbbs[k]]	<- nvec[delta.gbbs[k]] + 1
        beta[k] <- beta.vals[beta.ind]
        Xbeta <- Xbeta  + X[, k] * beta[k]
    }

    return(list(nvec       = nvec,
                beta       = beta,
                delta.gbbs = delta.gbbs,
                Xbeta      = Xbeta))
}

biomial.loglik <- function(Y, Xbeta, expXbeta){

    loglik = sum( Y * Xbeta -log1p(expXbeta))
    return(loglik)
}

# gibbs sampler derived from Yekuteili code
#'
#' @paran n.gibbs - (int ) number of gibbs samples to be performed
#' @param Y       -  ( n x 1) observastions
#' @param X       -  ( n x p)  covariates vector
#' @param a.in    -  (2 x 1) lower and upper point of the density
#' @param L       -  (int)   2^L is number of bins in the tree
#' @param beta0   -  (p x 1) inital guess of beta
#'
gibbs.logistic.Yekutieli <- function(n.gibbs,
                                 Y,
                                 X,
                                 a.int,
                                 L,
                                 beta){

    p			<- dim(X)[2]
    n			<- dim(X)[1]
    tree.data <- build.tree(L)

    I = tree.data$I


    a.vec	<- seq(a.int[1],
                   a.int[2],
                   length = I + 1)


    beta[beta< a.int[1]] = a.int[1] + 1e-8
    beta[beta> a.int[2]] = a.int[2] - 1e-8

    beta.vals	<- seq(a.int[1],
                       a.int[2],
                       length = I*20 + 1)
    num.vals <- length(beta.vals)
    beta.vals.a.ind	<- c(rep(1:I, each = 20), I)

    pi.gibbs		<- array(dim = c(n.gibbs, I))
    beta.gibbs	    <- array(dim = c(n.gibbs, p))
    delta.gibbs	    <- array(dim = c(n.gibbs, p))


    num.y.1		    <- sum(Y == 1)
    ind.y.1		    <- which(rep(Y, I + 1) == 1)
    beta.gibbs[1,]	<- beta
    pi.gibbs[1,]	<- log( (c(sum(beta.gibbs[1,] < a.vec[2]),diff(rank(c(a.vec[2:I],100,beta.gibbs[1,]))[1:I] -
                                                                  (1:I))) + 100/I) / 900)
    Xbeta	<- X %*% beta.gibbs[1,]

    for (g in 2:n.gibbs)
    {
        if (g %% 100 == 1)	print(paste(g))


        res.beta <- sample.gibbs.beta.Yekutieli(num.y.1,
                                                ind.y.1,
                                                beta.vals,
                                                I,
                                                Xbeta,
                                                X,
                                                beta.gibbs[g-1,],
                                                a.vec,
                                                pi.gibbs[g-1,beta.vals.a.ind],
                                                beta.vals.a.ind)



        beta.gibbs[g,] <- res.beta$beta
        delta.gibbs[g,] <- res.beta$delta.gbbs
        Xbeta <- res.beta$Xbeta

        pi.gibbs[g,] <- polya.gibs(res.beta$nvec, tree.data)
    }
    return(list(a.vec=  a.vec,
                pi.gibbs = pi.gibbs,
                beta.gibbs = beta.gibbs,
                delta.gibbs = delta.gibbs))

}



# gibbs sampler derived from Yekuteili code
#'
#' @paran n.gibbs - (int ) number of gibbs samples to be performed
#' @param Y       -  ( n x 1) observastions
#' @param X       -  ( n x p)  covariates vector
#' @param a.in    -  (2 x 1) lower and upper point of the density
#' @param L       -  (int)   2^L is number of bins in the tree
#' @param beta0   -  (p x 1) inital guess of beta
#'
gibbs.logistic           <- function(n.gibbs,
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
    for (g in 2:n.gibbs)
    {
        if (g %% 100 == 1)	print(paste(g))

        if(cpp ==F){
            res.beta <- sample.gibbs.beta (Y, Xbeta, X, beta.gibbs[g-1,], beta.ind, a.vec, pi.gibbs[g-1,])
            Xbeta <- res.beta$Xbeta
            beta.ind <- res.beta$beta.ind
            beta.gibbs[g,] <- res.beta$beta
            delta.gibbs[g, ] <- beta.ind
        }else{
            pi.gibbs_g       <-  pi.gibbs[g-1,]
            beta.gibbs.g     <- beta.gibbs[g-1,]
            res.bets.cpp     <- sample_gibbs_beta_logistic_cpp(Y, Xbeta, X, beta.gibbs.g, beta.ind, a.vec ,pi.gibbs_g, innersample= 10, interval_sample=3)
            beta.ind <-         res.bets.cpp$beta.ind
            beta.gibbs[g,]   <- beta.gibbs.g
            delta.gibbs[g, ] <- beta.ind
        }



        nvec <- table(factor(beta.ind, levels=1:I))
        pi.gibbs[g,] <- polya.gibs(nvec, tree.data)
    }
    return(list(a.vec=  a.vec,
                pi.gibbs = pi.gibbs,
                beta.gibbs = beta.gibbs,
                delta.gibbs = delta.gibbs))

}
