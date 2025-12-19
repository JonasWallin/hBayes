
library(NPBayes)
library(glmnet)
library(brms)
n.gibbs <- 110
sim <- 1
burnin  <- ceiling(0.1*n.gibbs)
n	<- 4000
p	<- 800

prior.dense <- prior(
    R2D2(
        mean_R2  = 0.9,
        prec_R2  = 1
    ),
    class = "b"
)
prior.sparse <- prior(
    R2D2(
        mean_R2  = 0.1,
        prec_R2  = 1
    ),
    class = "b"
)



set.seed(1)
sim.mse.mat <- array(dim=c(6, sim, 3))
X.mat	<- array(rnorm(n*p,mean = 0, sd = sqrt(1/n)),dim=c(n,p))
set.seed(1)
betas <- array(dim=c(p,3))
betas[,1] <- c(rep(-10,100),rep(10,100),rep(0,600))
betas[,2] <- rnorm(p,mean = 3, sd = 4)
betas[,3] <- rbinom(p,1,0.5)*rnorm(p,mean = 7, sd = 1)

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
        beta.mle <- fit$coeff
        beta.mle.adj <- beta.mle/1.499

        fit.glmnet	<- cv.glmnet(X.mat,y,alpha = 1,intercept = FALSE, family = "binomial",type.measure = "class")
        beta.lasso	<- as.numeric(coef(fit.glmnet, s = "lambda.min"))[-1]

        fit.glmnet	<- cv.glmnet(X.mat,y, alpha = 0,intercept = FALSE, family = "binomial",type.measure = "class")
        beta.ridge	<- as.numeric(coef(fit.glmnet, s = "lambda.min"))[-1]


        a.min <- min(c(-24,beta.mle - 0.5))
        a.max <- max(c(24,beta.mle + 0.5))


        hBeta <- gibbs.logistic(n.gibbs,
                                y,
                                X.mat,
                                c(a.min,a.max),
                                6,
                                beta = as.vector(beta.ridge),
                                cpp=T)
        beta.hBeta <- colMeans(hBeta$beta.gibbs[burnin:n.gibbs,])


        dat1 <- data.frame(y = y,X = X.mat)
        fit.R2D2.sparse <- brm(
            formula = y ~ -1 + .,
            data    = dat1,
            family  = bernoulli(link = "logit"),
            prior   = prior.sparse,
            chains  = 4,
            cores   = 4,
            iter    = 2000,
            seed    = 1
        )
        beta.R2D2.sparse <- fixef(fit.R2D2.sparse)[,1]
        fit.R2D2.dense <- brm(
            formula = y ~ -1 + .,
            data    = dat1,
            family  = bernoulli(link = "logit"),
            prior   = prior.dense,
            chains  = 4,
            cores   = 4,
            iter    = 2000,
            seed    = 1
        )
        beta.R2D2.dense <- fixef(fit.R2D2.dense)[,1]
        sim.mse.mat[1, i, j] <- mean((beta.mle    - betas[,j])^2)
        sim.mse.mat[2, i, j] <- mean((beta.oracle - betas[,j])^2)
        sim.mse.mat[3, i, j] <- mean((beta.hBeta  - betas[,j])^2)
        sim.mse.mat[4, i, j] <- mean((beta.lasso  - betas[,j])^2)
        sim.mse.mat[5, i, j] <- mean((beta.ridge  - betas[,j])^2)
        sim.mse.mat[6, i, j] <- mean((beta.mle.adj  - betas[,j])^2)
    }
}
