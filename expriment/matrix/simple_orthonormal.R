n.mcmc <- 1000
p <- 5*3
n	<- 100

set.seed(0)

Lambda  <- diag(p:1)
Sigma <- Lambda
x.smp	<- rmvnorm(n,rep(0,p),Sigma)
Eigen.Sigma_hat <- eigen(t(x.smp) %*% (x.smp)/n)
lambda.hat <- Eigen.Sigma_hat$values

# Run sampler
#res <- Eigen.sampler.bingham(x.smp,n.mcmc=n.mcmc,n.E.mcmc = 1e1,L = 6, sample.eigenvalues=T,shuffle = T)

res <- Eigen.sampler.old(x.smp,n.mcmc=n.mcmc,n.E.mcmc = 1e2,L = 6, sample.eigenvalues=T)

# Display eigenvalue estimation

boxplot(split(res$lambda[(n.mcmc-199):n.mcmc,],rep(1:p,each = 200)),ylim = range(c(res$lambda[(n.mcmc-199):n.mcmc,],lambda.hat)))
lines(lambda.hat,type = "b",col=2,pch = 3)
lines(diag(Lambda),type = "b",col=3,pch=3)
