
library(NPBayes)
library(glmnet)
n.gibbs <- 102
sim <- 3
burnin  <- ceiling(0.1*n.gibbs)
n	<- 4000
p	<- 800

set.seed(1)

sim.mse.mat <- array(dim=c(5, sim, 3))
X.mat	<- array(rnorm(n*p,mean = 0, sd = sqrt(1/n)),dim=c(n,p))

betas <- array(dim=c(p,3))
betas[,1] <- c(rep(-10,100),rep(10,100),rep(0,600))
betas[,2] <- rnorm(K,mean = 3, sd = 4)
betas[,3] <- rbinom(K,1,0.5)*rnorm(K,mean = 7, sd = 1)

inv.log.odds <- function(Xbeta){
    return(1/(1+ exp(-Xbeta)))
}

for(i in 1:sim){
    for(j in 1:3){
        cat('(i,j) = (',i,',',j,')\n')
        y		  <- rbinom(n,1,inv.log.odds(X.mat %*% betas[,j]))
        beta.oracle <- gibbs.logistic.permute.oracle(n.gibbs, y, X.mat, betas[,j],F)
        beta.oracle <- colMeans(beta.oracle[burnin:n.gibbs,])

        fit	    <- glm(y ~ X.mat + 0,family = binomial)
        beta.lse <- fit$coeff


        fit.glmnet	<- cv.glmnet(X.mat,y,alpha = 1,intercept = FALSE, family = "binomial",type.measure = "class")
        beta.lasso	<- as.numeric(coef(fit.glmnet, s = "lambda.min"))[-1]

        fit.glmnet	<- cv.glmnet(X.mat,y, alpha = 0,intercept = FALSE, family = "binomial",type.measure = "class")
        beta.ridge	<- as.numeric(coef(fit.glmnet, s = "lambda.min"))[-1]


        a.min <- min(c(-24,beta.lse - 0.5))
        a.max <- max(c(24,beta.lse + 0.5))

        hBeta <- gibbs.logistic(n.gibbs,y,X.mat,c(a.min,a.max),6,as.vector(beta.ridge),cpp=T)
        beta.hBeta <- colMeans(hBeta$beta.gibbs[burnin:n.gibbs,])
        sim.mse.mat[1, i, j] <- mean((beta.lse    - betas[,j])^2)
        sim.mse.mat[2, i, j] <- mean((beta.oracle - betas[,j])^2)
        sim.mse.mat[3, i, j] <- mean((beta.hBeta  - betas[,j])^2)
        sim.mse.mat[4, i, j] <- mean((beta.lasso  - betas[,j])^2)
        sim.mse.mat[5, i, j] <- mean((beta.ridge  - betas[,j])^2)
    }
}
