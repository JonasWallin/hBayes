
library(mvtnorm)

p <- 2
n <- 25

set.seed(14)
Q <- smp.rnd.orth.mat(p)

lambda  <- 3*sort(rchisq(p,4),dec = T)
#lambda  <- c(10,10,9,rep(1,p-3))
#lambda  <- p:1
Lambda  <- diag(lambda)
Sigma <- (Q) %*% Lambda %*% t(Q)
x.smp <- rmvnorm(n,rep(0,p),Sigma)
sigma.hat <- t(x.smp) %*% (x.smp)/n

Eigen.Sigma_hat <- eigen(sigma.hat)
lambda.hat <- Eigen.Sigma_hat$values
gamma.hat <- Eigen.Sigma_hat$vectors

y.smp <- x.smp %*% gamma.hat
res_oracle <- Eigen.sampler.old(y.smp,n.mcmc=30, n.E.mcmc = 10,L = 6, sample.eigenvalues=F, sample.eigenvector = T,lambda = lambda,verbose = TRUE)
