

test_that("test simple eigenvalue distribution check", {
    set.seed(2)
p <- 3
n <- 1000
n.mcmc <- 100
Sigma <- diag(sample(1:p))
X <- t(chol(Sigma)%*%matrix(rnorm(n*p), nrow=p,ncol=n))

res <- Eigen.sampler(X,
                     n.mcmc=n.mcmc,
                     n.E.mcmc = 1e1,
                     L = 6)
E.mean <- matrix(0,p,p)
D <- diag(p)
for(i in 1:n.mcmc){
    for(j in 1:p)
        D[j,j] <-  res$gamma[which.max(abs(res$gamma[,j,i])),j,i]

    E.mean <- E.mean + res$gamma[,,i]%*%D
}
E.mean <- E.mean/n.mcmc

expect_equal(diag(t(round(E.mean))%*%Sigma%*%round(E.mean)), sort(diag(Sigma), decreasing = T))
})
