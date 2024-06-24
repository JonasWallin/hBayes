##
# student exam score
#
##
library(NPBayes)
library(devtools)
set.seed(1)
graphics.off()
#devtools::load_all()
mcmc.samples = 2000
X <- as.matrix(read.table('data/student_score.txt',header=T))
X <- scale(X, center=T,scale=F)
n <- dim(X)[1]
E <-  eigen(t(X)%*%X/n)
set.seed(1)
res <- oracle.metrop.sampler(E$vectors,X, E$values, mcmc.samples)
set.seed(1)
res_old <- oracle.metrop.sampler.old(E$vectors,X, E$values, mcmc.samples)
plot(cumsum(res$acc[,2])/(1:length(res$acc[,2])),ylim=c(0,1),type='l')
for(i in 3:5){
    lines(cumsum(res$acc[,i])/(1:length(res$acc[,i])))
}
for(i in 2:5){
    cat('sigma[',i,'] = ', round(res$MH.objs[[i]]$sigma,3),' n = ',res$MH.objs[[i]]$n.SA,'\n')
}
print(max(abs(res$loglik-res_old$loglik)))


B <- -0.5*diag(1/E$values)
A <- (t(X)%*%X)
res.bin<-rbing(mcmc.samples, A, B)
Es <- res.bin$Es
#likelihood should be approx
lik_bin <- rep(0,mcmc.samples)
p <- dim(A)[1]
for(i in 1:mcmc.samples)
    lik_bin[i] <-  sum(diag(Es[,,i]%*%B%*%t(Es[,,i])%*%A)) - n*0.5*sum(log(E$values)) -  0.5*n*p*log(2*pi)


res <- Eigen.sampler(X,
                     n.mcmc=1e3,
                     n.E.mcmc = 1e1,
                     L = 6)

res.old <- Eigen.sampler.old(X,
                     n.mcmc=1e3,
                     n.E.mcmc = 1e1,
                     L = 6)
res.2 <- Eigen.sampler.bingham.v2(X,
                                  n.mcmc=1e3,
                                  n.E.mcmc = 1e1/p,
                                  L = 6)
