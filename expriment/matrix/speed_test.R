
library(mvtnorm)
library(MASS)
graphics.off()
set.seed(1)
run.old <- FALSE
mcmc.samples <- 10000
p <- 40
n <- 100
df <- 5  # degrees of freedom
S <- solve(toeplitz(c(2,1,rep(0,p-2))))

# Generate a random covariance matrix from the Wishart distribution
Sigma <- rWishart(n = 1, df = p+df, Sigma = S)[,,1]
X	<- rmvnorm(n,rep(0,p),Sigma)
E <-  eigen(t(X)%*%X/n)

set.seed(1)
timing <- system.time(res <- oracle.metrop.sampler(E$vectors,X, E$values, mcmc.samples))

plot(cumsum(res$acc[,2])/(1:length(res$acc[,2])),ylim=c(0,1),type='l')
for(i in 3:p){
    lines(cumsum(res$acc[,i])/(1:length(res$acc[,i])))
}
for(i in 2:p){
    cat('sigma[',i,'] = ', round(res$MH.objs[[i]]$sigma,3),' n = ',res$MH.objs[[i]]$n.SA,'\n')
}
if(run.old){
    set.seed(1)
   timing.old <- system.time(res_old <- oracle.metrop.sampler.old(E$vectors,X, E$values, mcmc.samples))
   cat('likelihood diff (should be zero)= ',max(abs(res$loglik-res_old$loglik)),'\n')
   cat('time_old = ',round(timing.old[3]),'\n')

}
cat('time_new = ',round(timing[3]),'\n')
