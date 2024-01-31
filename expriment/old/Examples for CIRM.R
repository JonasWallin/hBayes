
#install.packages(c("Rcpp"))

library(NPBayes)
library(glmnet)


##################################################################################################################

#   Run the same linear models with y.sd = 4, 2, 1

##################################################################################################################


#      Setup:

n <- 4000
k <- 800
burnin  <- 1000
n.gbbs  <- 2000
set.seed(0)

X.mat	<- array(rnorm(n*k,mean = 0, sd = sqrt(1/n)),dim=c(n,k))
beta.vec   <- c(rep(-10,100),rep(0,k-200), rep(10,100))
mu.vec     <- X.mat %*% beta.vec

y1.vec       <- mu.vec + rnorm(n,mean = 0, sd = 4)
y2.vec       <- mu.vec + rnorm(n,mean = 0, sd = 2)
y3.vec       <- mu.vec + rnorm(n,mean = 0, sd = 1)

##################################################################################################################

# Simulation 1:   y.sd = 4

lm.full     <- lm(y1.vec ~ X.mat + 0)
beta1.lse    <- lm.full$coef



res.1 <- gibbs.normal(n.gbbs,
                    y1.vec,
                    X.mat,
                    c(-20,20),
                    8,
                    as.vector(beta1.lse),cpp=T)

res.1.oracle <- gibbs.normal.permute.beta.fixed.sigma(n.gbbs, y1.vec,X.mat, sigma = 4,c(-20,20),8, beta.vec[sample(k)])
pi.cdf      <- rbind(rep(0,1000),apply(exp(res.1$pi.gibbs[burnin:n.gbbs,]),1,cumsum))
pi.quant    <- apply(pi.cdf,1,quantile,c(0.025,0.50,0.975))


#   Pi estimation plot

x.lm <- range(beta1.lse)
plot(ecdf(beta.vec),col = 2,lwd = 2,xlim = x.lm, main = "Error SD = 4",ylab = "CDF",xlab = "Regression coeff.")
lines(res.1$a.vec,pi.quant[2,],col=3,lty = 2,lwd = 2)
lines(res.1$a.vec,pi.quant[1,],col=3,lty = 1,lwd = 2)
lines(res.1$a.vec,pi.quant[3,],col=3,lty = 1,lwd = 2)
lines(ecdf(beta1.lse),col = 4,lwd = 0.5)

beta1.hBeta   <- apply(res.1$beta.gibbs[burnin:n.gbbs,],2,mean)
beta1.oracle   <- apply(res.1.oracle$beta.gibbs[burnin:n.gbbs,],2,mean)

best_lambda <- cv.glmnet(X.mat, y1.vec, alpha = 0,standardize=F,intercept=F)$lambda.min
beta1.ridge <-  glmnet(X.mat, y1.vec, alpha = 0, lambda = best_lambda,standardize=F,intercept=F)$beta

best_lambda <- cv.glmnet(X.mat, y1.vec, alpha = 1,standardize=F,intercept=F)$lambda.min
beta1.glmnet <- glmnet(X.mat, y1.vec, alpha = 1,lambda = best_lambda,standardize=F,intercept=F)$beta
beta1.lasso <- (beta1.glmnet+1)@x - 1

#   Shrinkage plot

plot(beta1.lse,beta.vec,ylim = range(beta1.lse),xlim = range(beta1.lse),main = "Error SD = 4",
     xlab = "Max. Likelihood Est.",ylab = "Regression coeff.")
points(beta1.lse,beta1.hBeta,col=3)
points(beta1.lse,beta1.lasso,col=4)
points(beta1.lse,beta1.ridge,col=2)
abline(0,1)

#  Beta estimation comparison

sqrt(mean((beta.vec - beta1.lse)^2))
sqrt(mean((beta.vec - beta1.hBeta)^2))
sqrt(mean((beta.vec - beta1.oracle)^2))
sqrt(mean((beta.vec - beta1.ridge)^2))
sqrt(mean((beta.vec - beta1.lasso)^2))

##################################################################################################################

# Simulation 2:   y.sd = 2

lm.full     <- lm(y2.vec ~ X.mat + 0)
beta2.lse    <- lm.full$coef

res.2 <- gibbs.normal(n.gbbs,
                      y2.vec,
                      X.mat,
                      c(-20,20),
                      8,
                      as.vector(beta2.lse),cpp=T)

res.2.oracle <- gibbs.normal.permute.beta.fixed.sigma(n.gbbs, y2.vec,X.mat, sigma = 2,c(-20,20),8, beta.vec[sample(k)])

pi.cdf      <- rbind(rep(0,1000),apply(exp(res.2$pi.gibbs[burnin:n.gbbs,]),1,cumsum))
pi.quant    <- apply(pi.cdf,1,quantile,c(0.025,0.50,0.975))


#   Pi estimation plot

plot(ecdf(beta.vec),col = 2,lwd = 2,xlim = x.lm, main = "Error SD = 2",ylab = "CDF",xlab = "Regression coeff.")
lines(res.2$a.vec,pi.quant[2,],col=3,lty = 2,lwd = 2)
lines(res.2$a.vec,pi.quant[1,],col=3,lty = 1,lwd = 2)
lines(res.2$a.vec,pi.quant[3,],col=3,lty = 1,lwd = 2)
lines(ecdf(beta2.lse),col = 4,lwd = 0.5)

beta2.hBeta   <- apply(res.2$beta.gibbs[burnin:n.gbbs,],2,mean)
beta2.oracle   <- apply(res.2.oracle$beta.gibbs[burnin:n.gbbs,],2,mean)

best_lambda <- cv.glmnet(X.mat, y2.vec, alpha = 0,standardize=F,intercept=F)$lambda.min
beta2.ridge <-  glmnet(X.mat, y2.vec, alpha = 0, lambda = best_lambda,standardize=F,intercept=F)$beta

best_lambda <- cv.glmnet(X.mat, y2.vec, alpha = 1,standardize=F,intercept=F)$lambda.min
beta2.glmnet <- glmnet(X.mat, y2.vec, alpha = 1,lambda = best_lambda,standardize=F,intercept=F)$beta
beta2.lasso <- (beta2.glmnet+1)@x - 1

#   Shrinkage plot

plot(beta2.lse,beta.vec,xlim = range(beta2.lse),ylim = range(beta2.lse),main = "Error SD = 2",
     xlab = "Max. Likelihood Est.",ylab = "Regression coeff.")
points(beta2.lse,beta2.hBeta,col=3)
points(beta2.lse,beta2.lasso,col=4)
points(beta2.lse,beta2.ridge,col=2)
abline(0,1)

#  Beta estimation comparison

sqrt(mean((beta.vec - beta2.lse)^2))
sqrt(mean((beta.vec - beta2.hBeta)^2))
sqrt(mean((beta.vec - beta2.oracle)^2))
sqrt(mean((beta.vec - beta2.ridge)^2))
sqrt(mean((beta.vec - beta2.lasso)^2))


##################################################################################################################


# Simulation 3:   y.sd = 1

lm.full     <- lm(y3.vec ~ X.mat + 0)
beta3.lse    <- lm.full$coef

res.3 <- gibbs.normal(n.gbbs,
                      y3.vec,
                      X.mat,
                      c(-20,20),
                      8,
                      as.vector(beta3.lse),cpp=T)
res.3.oracle <- gibbs.normal.permute.beta.fixed.sigma(n.gbbs, y3.vec,X.mat, sigma = 1,c(-20,20),8, beta.vec[sample(k)])
pi.cdf      <- rbind(rep(0,1000),apply(exp(res.3$pi.gibbs[burnin:n.gbbs,]),1,cumsum))
pi.quant    <- apply(pi.cdf,1,quantile,c(0.025,0.50,0.975))


#   Pi estimation plot

plot(ecdf(beta.vec),col = 2,lwd = 2,xlim = x.lm, main = "Error SD = 1",ylab = "CDF",xlab = "Regression coeff.")
lines(res.3$a.vec,pi.quant[2,],col=3,lty = 2,lwd = 2)
lines(res.3$a.vec,pi.quant[1,],col=3,lty = 1,lwd = 2)
lines(res.3$a.vec,pi.quant[3,],col=3,lty = 1,lwd = 2)
lines(ecdf(beta3.lse),col = 4,lwd = 0.5)

beta3.hBeta   <- apply(res.3$beta.gibbs[burnin:n.gbbs,],2,mean)
beta3.oracle   <- apply(res.3.oracle$beta.gibbs[burnin:n.gbbs,],2,mean)


best_lambda <- cv.glmnet(X.mat, y3.vec, alpha = 0,standardize=F,intercept=F)$lambda.min
beta3.ridge <-  glmnet(X.mat, y3.vec, alpha = 0, lambda = best_lambda,standardize=F,intercept=F)$beta

best_lambda <- cv.glmnet(X.mat, y3.vec, alpha = 1,standardize=F,intercept=F)$lambda.min
beta3.glmnet <- glmnet(X.mat, y3.vec, alpha = 1,lambda = best_lambda,standardize=F,intercept=F)$beta
beta3.lasso <- (beta3.glmnet+1)@x - 1

#   Shrinkage plot

plot(beta3.lse,beta.vec,ylim = range(beta3.lse),main = "Error SD = 1",xlab = "Max. Likelihood Est.",ylab = "Regression coeff.")
points(beta3.lse,beta3.hBeta,col=3)
points(beta3.lse,beta3.lasso,col=4)
points(beta3.lse,beta3.ridge,col=2)
abline(0,1)

#  Beta estimation comparison

sqrt(mean((beta.vec - beta3.lse)^2))
sqrt(mean((beta.vec - beta3.oracle)^2))
sqrt(mean((beta.vec - beta3.hBeta)^2))
sqrt(mean((beta.vec - beta3.ridge)^2))
sqrt(mean((beta.vec - beta3.lasso)^2))


