


library(NPBayes)
library(mvtnorm)
library(MASS)
library(stcov)
library(Matrix)

#graphics.off()

#source("aux_functions_code.R")


###################################################################################################
###################################################################################################

# A relatively easy covariance matrix estimation problem in which the sample covariance matrix
# does a very good job and the sampler miscalculates the eigenvalues

###################################################################################################
###################################################################################################


n.mcmc <- 3000
p <- 50
n	<- 60

set.seed(11)
Q <- smp.rnd.orth.mat(p)

Lambda  <- diag((p:1))
#Lambda  <- diag(c(rep(10,p/2),rep(1,p/2)))
Sigma <- (Q) %*% Lambda %*% t(Q)
x.smp	<- rmvnorm(n,rep(0,p),Sigma)
Eigen.Sigma_hat <- eigen(t(x.smp) %*% (x.smp)/n)

MLE.loglik <- oracle.gamma.R.loglik(chol(t(x.smp) %*% x.smp),Eigen.Sigma_hat$vector,Eigen.Sigma_hat$values,n)
Sigma.loglik <- oracle.gamma.R.loglik(chol(t(x.smp) %*% x.smp),Q,diag(Lambda),n)

lambda.hat <- Eigen.Sigma_hat$values

res_old <- Eigen.sampler.old(x.smp,n.mcmc=n.mcmc, n.E.mcmc = 1e1,L = 6, sample.eigenvalues=T, sample.eigenvector = T,verbose = TRUE)

boxplot(split(res_old$lambda[(n.mcmc-1000+1):n.mcmc,],rep(1:p,each = 1000)),ylim = range(lambda.hat),main  = "old sampler: Sorted eigenvalue estimation")
lines(lambda.hat,type = "b",col=2,pch = 3,lwd = 2)
lines(diag(Lambda),type = "b",col=3,pch=3,lwd = 2)


exp.a.vec <- exp(res_old$lambda.param$a_vec)
a.vec <- res_old$lambda.param$a_vec

plot(ecdf(log(lambda.hat)),col=2,main = "Eigenvalue CDF log-scale estimation")
lines(ecdf(log(diag(Lambda))),col=3)
cumsum.pi <- cbind(0,t(apply(res_old$pi[(n.mcmc-1000+1):n.mcmc,],1,cumsum)))
lines(a.vec,apply(cumsum.pi,2,quantile,prob=0.5),type = "l")
lines(a.vec,apply(cumsum.pi,2,quantile,prob=0.25),type = "l",lty=2)
lines(a.vec,apply(cumsum.pi,2,quantile,prob=0.75),type = "l",lty=2)


loglik.profile <- rep(NA,dim(res_old$lambda)[1])
for(i in 1:dim(res_old$lambda)[1]) loglik.profile[i] <- oracle.gamma.R.loglik(chol(t(x.smp) %*% x.smp),res_old$gamma[,,i],res_old$lambda[i,],n)
plot(loglik.profile,main = "Loglikelihood old",ylab = "Loglikelihood", xlab = "MCMC iteration",
     ylim = range(MLE.loglik,Sigma.loglik,loglik.profile))
abline(h = MLE.loglik,col=2)
abline(h = Sigma.loglik,col=3)
