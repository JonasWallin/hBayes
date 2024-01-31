library(NPBayes)
library(MASS)
p <- 20
n	<- 100


#' Build a orthonormal matrix
#' using Haar measures.
#' @param  n -  dimension of the orthonormal matrix
#' @return P.tmp - nxn orthonormal matrix
smp.rnd.orth.mat	<- function(n = 10)
{

    P.tmp	<- R.rho.b(runif(1,0,2*pi),1-2*rbinom(1,1,0.5))

    if(2 < n)
        for(i in 3:n)
        {
            v.i		<- smp.v(i)
            hr		<- house.reflctn(v.i)
            P.tmp	<- hr %*% cbind(c(1,rep(0,i-1)),rbind(0,P.tmp))
        }

    return(P.tmp)
}

Gamma <- smp.rnd.orth.mat(p)
D.oracle <- p:1
Sigma <- Gamma %*% diag(D.oracle) %*% t(Gamma)

cov.vec <- rep(NA,1000)

x.smp	<- rmvnorm(n,rep(0,p),Sigma)
sigma.mle <- t(x.smp) %*% x.smp / n
gamma.mle <- eigen(sigma.mle)$vectors

aa <- oracle.metrop.sampler(x.smp,n.mtrp = 10^4)
