##
# power simulation for mixed effect model
# to create Table 3: set.seed(2) q=1000, table =2 ,gives sparse mixed effect
# to create Table 4: set.seed(2) q=1000, table =1 ,gives laplace distribution
# D: 2022-02-15
##

rm(list=ls())
save.file=F
graphics.off()
library(bigstep)
library(PolyMixed)
library(RcppEigen)
library(ggplot2)
library(SSLASSO)
library(glmnet)
library(bigstep)
set.seed(2)
distribution = c(1)  # 1 - gaussian random effect, 2- Laplace randeom effect, 3 -sparse
q     <- 10    # number of simulations
#thin =  how often to observe X (i.e. cm distnace)
#type =  1 - pure measurement error (\sigma^2I )
#        2 -  mixed effect     (\sigma^2I + XX^T \tau )


thin    = 1
signal  = 1 #signal = 1 - strong signal , 0 weak signal, -1 no signal
active_propotion = 0.2 # if distribution 3- then how many non-zero signals
# parameters----------------------------------------------------------------
chr <- 10                # chromosomes
m <- 150                 # markers on chromosome
M <- chr * m             # all markers
markers <- rep(m, chr)   # markers on each chromosome
ns <- c(200)             # observations
d <- 1                   # distance between markers (cM)
sigma <- 1               # sd for error term
tau <- 0.01              # sd for polygenic effect
mu <- rep(0.004, q)      # mean for polygenic effect


nc = NULL
qtl.pos <- c(151, 824, 1274)
mag = 0.35
if(signal==1){
    mag = 0.5
}else if(signal==-1){
    mag = 0
}
qtl  <- c(mag, qtl.pos[1], 1) # qtl (beta; position: marker, trait or traits)
qtl2 <- c(mag, qtl.pos[2], 1)
qtl3 <- c(-mag, qtl.pos[3], 1)
beta.true <- c(qtl[1], qtl2[1], qtl3[1])

# --------------------------------------------------------------------------

t.mix.forward.known <- matrix(0, q, M)
t.mix.forward <- matrix(0, q, M)
sigma.est.forward <- numeric(q)
tau.est.forward <- numeric(q)
crit.mix.forward <- numeric(q)
ignoreTau <- numeric(q) # LRT test



for(n in ns){
    string_out <- paste('\n\n',n)
    Xs <- list()
    ys <- list()
    betas <- list()
    betas.lasso <- betas.SSlasso <-  list()
    dupls <- list()
    SVDXs <- list()
    for(i in 1:q) {
        # Should we allow for duplitacte?
        X <- simulateDesignMatrix(chr, m, n, d)

        Z <- rnorm(dim(X)[2])
        beta_          <- rep(0, dim(X)[2])
        if(distribution==1){
            beta_    = tau * Z
            beta_[qtl.pos] <- beta.true
        }else if(distribution==2){
            V <- rgamma(dim(X)[2],1)
            beta_ <- mu[1] + tau * sqrt(V)*Z
            beta_[qtl.pos] <- beta.true

        }else if(distribution == 3){
            index    = sample(1:dim(X)[2],ceiling(dim(X)[2]*active_propotion))
            V        = rep(0, dim(X)[2])
            V[index] = 1
            beta_    = V*mu[1]/active_propotion + tau/active_propotion * sqrt(V)*Z
            beta_[qtl.pos] <- beta.true
        }

        y <- X%*%(beta_) + sigma * rnorm(dim(X)[1])
        #y <- simulateTraits(X, q = 1, mu, tau, qtl = 0, qtl2 = 0, qtl3 = 0)

        thin_index <- seq(1,dim(X)[2],by = thin)
        X          <- X[,thin_index]
        Xs[[i]]    <- X
        ys[[i]]    <- y
        betas[[i]] <- beta_
        dupls[[i]] <- findDuplicate(X)
        SVDXs[[i]] <- svd(X, nu = dim(X)[1])


    }
    for(i in 1:q){
        # Lasso
        cvfit <- cv.glmnet(Xs[[i]], ys[[i]])
        coef.lasso <- coef(cvfit, s = "lambda.1se")
        index.lasso <- which(coef.lasso!=0)
        index.lasso <- index.lasso[-1]
        betas.lasso[[i]] <-cbind(index.lasso-1, coef.lasso[index.lasso])

        #Spike and slab lasso
        sslasso.fit <- SSLASSO(Xs[[i]], ys[[i]], variance = "unknown")
        index <- sslasso.fit$model
        betas.SSlasso[[i]] <- cbind(index, sslasso.fit$beta[sslasso.fit$model,ncol(sslasso.fit$beta)])
    }
}
