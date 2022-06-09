##
# power simulation for mixed effect model
# to create Table 3: set.seed(2) q=1000, table =2 ,gives sparse mixed effect
# to create Table 4: set.seed(2) q=1000, table =1 ,gives laplace distribution
# D: 2022-02-15
##

rm(list=ls())
save.fig=F
graphics.off()
library(NPBayes)
library(bigstep)
library(PolyMixed)
library(RcppEigen)
library(ggplot2)
library(SSLASSO)
library(glmnet)
library(bigstep)
set.seed(10)
distribution = c(1)  # 1 - gaussian random effect, 2- Laplace randeom effect, 3 -sparse
shift <- -0.01 #-0.02 # shift in the laplace distribution
q     <- 10    # number of simulations
#thin =  how often to observe X (i.e. cm distnace)
#type =  1 - pure measurement error (\sigma^2I )
#        2 -  mixed effect     (\sigma^2I + XX^T \tau )


thin    = 1
signal  = 0 #signal = 1 - strong signal , 0 weak signal, -1 no signal
active_propotion = 0.2 # if distribution 3- then how many non-zero signals
# parameters----------------------------------------------------------------
chr <- 10                 # chromosomes
m <-  150                 # markers on chromosome
M <- chr * m             # all markers
markers <- rep(m, chr)   # markers on each chromosome
ns <-  400            # observations
d <- 1                   # distance between markers (cM)
sigma <- 0.1               # sd for error term
tau <- 0.01              # sd for polygenic effect
mu <- rep(0.01, q)      # mean for polygenic effect


n.gibbs <- 1000 #samples for MCMC

nc = NULL
qtl.pos <- ceiling(M * c(0.2,0.5,0.8))
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
    Xs <- Xs_new <- list()
    ys <- ys_new <- ys_new2 <- list()
    betas <- list()
    betas.lasso <- betas.SSlasso <- betas.ridge <- beta.npbayes <- list()
    intercept.ridge <- intercept.lasso <- intercept.sslasso <- vector(length=q, mode="numeric")
    dupls <- list()
    SVDXs <- list()
    for(i in 1:q) {
        # Should we allow for duplitacte?
        X <- simulateDesignMatrix(chr, m, n, d)
        X_new <- simulateDesignMatrix(chr, m, n, d)

        Z <- rnorm(dim(X)[2])
        beta_          <- rep(0, dim(X)[2])
        if(distribution==1){
            beta_    = tau * Z
            beta_[qtl.pos] <- beta.true
        }else if(distribution==2){
            V <- rgamma(dim(X)[2],1)
            beta_ <- mu[1] + shift * (V-1) + tau * sqrt(V)*Z
            beta_[qtl.pos] <- beta.true

        }else if(distribution == 3){
            index    = sample(1:dim(X)[2],ceiling(dim(X)[2]*active_propotion))
            V        = rep(0, dim(X)[2])
            V[index] = 1
            beta_    = V*mu[1]/active_propotion + tau/active_propotion * sqrt(V)*Z
            beta_[qtl.pos] <- beta.true
        }

        y <- X%*%(beta_) + sigma * rnorm(dim(X)[1])
        y_new <- X%*%(beta_) + sigma * rnorm(dim(X)[1])
        y_new_X_new <- X_new%*%(beta_) + sigma * rnorm(dim(X)[1])
        #y <- simulateTraits(X, q = 1, mu, tau, qtl = 0, qtl2 = 0, qtl3 = 0)

        thin_index   <- seq(1,dim(X)[2],by = thin)
        X            <- X[,thin_index]
        X_new        <- X_new[,thin_index]
        Xs[[i]]      <- X
        Xs_new[[i]]  <- X_new
        ys[[i]]      <- y
        ys_new[[i]]  <- y_new
        ys_new2[[i]] <- y_new_X_new
        betas[[i]]   <- beta_


    }
    for(i in 1:q){
        # Lasso
        cvfit <- cv.glmnet(Xs[[i]], ys[[i]])
        coef.lasso <- coef(cvfit, s = "lambda.1se")
        index.lasso <- which(coef.lasso!=0)
        index.lasso <- index.lasso[-1]
        betas.lasso[[i]] <-cbind(index.lasso-1, coef.lasso[index.lasso])
        intercept.lasso[i] <-coef.lasso[1]

        #Spike and slab lasso
        sslasso.fit <- SSLASSO(Xs[[i]], ys[[i]], variance = "unknown")
        index <- sslasso.fit$model
        betas.SSlasso[[i]] <- cbind(index, sslasso.fit$beta[sslasso.fit$model,ncol(sslasso.fit$beta)])
        intercept.sslasso[i] <- sslasso.fit$intercept[1,dim(sslasso.fit$intercept)[2]]
        #ridge regression
        cvfit <- cv.glmnet(Xs[[i]], ys[[i]], alpha=0)
        coef.ridge <- coef(cvfit, s = "lambda.1se")
        betas.ridge[[i]] <- coef.ridge[-1]
        intercept.ridge[i] <- coef.ridge[1]
        #hBayes
        #res <- gibbs.normal(n.gibbs,
        #                    ys[[i]],
        #                    Xs[[i]],
        #                    c(-0.5,0.5),
        #                    7,
        #                    sslasso.fit$beta[,ncol(sslasso.fit$beta)])
        res <- gibbs.normal.fixed.sigma(n.gibbs,
                            ys[[i]],
                            Xs[[i]],
                            sigma,
                            c(-0.5,0.5),
                            8,
                            sslasso.fit$beta[,ncol(sslasso.fit$beta)])

        beta.npbayes[[i]] <-  colMeans(res$beta.gibbs[floor(n.gibbs/3):n.gibbs,])
    }
}
library(zoo)
library(ggplot2)
window.size= 10
#' @param x.sim simulations to get posterior conf int
plot.graph <- function(x,
                       window.size,
                       markers,
                       beta_true= NULL,
                       x.sim = NULL){
    x.window <- rep(0,length(x))
    beta_true.window = rep(0, length(x))
    marker.tot <- c(0, cumsum(markers))
    CI <- matrix(0, ncol=2, nrow=length(x))
    for(i in 2:length(marker.tot)){
        x.window[(1+marker.tot[i-1]):marker.tot[i]] <- zoo::rollsum(x[(1+marker.tot[i-1]):marker.tot[i]],
                                                                window.size,
                                                                na.pad=T,
                                                                fill="extend")

        if(!is.null(beta_true))
            beta_true.window[(1+marker.tot[i-1]):marker.tot[i]] <- zoo::rollsum(beta_true[(1+marker.tot[i-1]):marker.tot[i]],
                                                                    window.size,
                                                                    na.pad=T,
                                                                    fill="extend")

        if(!is.null(x.sim)){
            CI.window = t(apply(x.sim[,(1+marker.tot[i-1]):marker.tot[i]], 1,zoo::rollsum,k= window.size,
                              na.pad=T,
                              fill="extend"))
            CI[(1+marker.tot[i-1]):marker.tot[i],] = t(apply(CI.window, 2, quantile, probs=c(0.05,0.95)))
        }
    }
    df <- data.frame(locus = 1:length(x.window), beta = x.window)
    fig <- ggplot(df)
    if(!is.null(beta_true))
        fig <- fig + geom_line(data=data.frame(locus = 1:length(x.window), beta = beta_true.window),
                               aes(x=locus, y= beta),
                               colour ="red",
                               alpha = 0.5)
    fig <- fig + geom_line(aes(x=locus, y=beta))
    fig <- fig +  geom_vline(xintercept = marker.tot,
                             colour ="blue",
                             linetype = "dashed",
                             alpha = 0.5)
    if(!is.null(x.sim)){
        fig <- fig  + geom_ribbon(data = data.frame(locus=1:length(x),
                                                          low = CI[,1]
                                                          ,upp = CI[,2]),
                                        aes(x=locus,ymin = low, ymax=upp), alpha=0.2, fill='blue')
    }

    return(list(fig=fig, res = beta_true.window-x.window ) )
}
library(latex2exp)
res.lasso <- plot.graph(coef.lasso[-1],window.size, markers, beta_true = beta_)
fig.lasso <- res.lasso$fig
fig.lasso <- fig.lasso +  labs(title=paste('MSE = ',round(sqrt(mean(res.lasso$res^2)),3) , " (window size = ", window.size,")",sep=""))
fig.lasso <- fig.lasso+ theme(plot.title = element_text(hjust = 0.5))

print(fig.lasso)
if(save.fig)
    ggsave('lasso.jpeg',fig.lasso)

res.sslasso <- plot.graph(sslasso.fit$beta[,ncol(sslasso.fit$beta)],window.size, markers, beta_true = beta_)
fig.sslasso <- res.sslasso$fig
fig.sslasso <- fig.sslasso +  labs(title=paste('MSE = ',round(sqrt(mean(res.sslasso$res^2)),3) , " (window size = ", window.size,")",sep=""))
fig.sslasso <- fig.sslasso+ theme(plot.title = element_text(hjust = 0.5))

print(fig.sslasso)
if(save.fig)
    ggsave('sslasso.jpeg',fig.sslasso)


Betas.np <- colMeans(res$beta.gibbs[floor(n.gibbs/3):n.gibbs,])
res.np <- plot.graph(Betas.np,window.size, markers, beta_true = beta_, x.sim =res$beta.gibbs[floor(n.gibbs/3):n.gibbs,] )
fig.np <- res.np$fig
fig.np <- fig.np +  labs(title=paste('MSE = ',round(sqrt(mean(res.np$res^2)),3) , " (window size = ", window.size,")",sep=""))
fig.np <- fig.np+ theme(plot.title = element_text(hjust = 0.5))

print(fig.np)
if(save.fig)
    ggsave('bayes.jpeg',fig.np)
#geom_ribbon(aes(ymin = Anomaly10y - Unc10y, ymax = Anomaly10y + Unc10y), alpha = 0.2)
cat('lasso var:',round(var(zoo::rollsum(coef.lasso[-1],window.size )-zoo::rollsum(beta_,window.size)),4),'\n')
cat('nbayes var:',round(var( zoo::rollsum(colMeans(res$beta.gibbs[floor(n.gibbs/3):n.gibbs,]),window.size)-zoo::rollsum(beta_,window.size)),4),'\n')



if(M < n){
    LS <- solve(t(X)%*%X,t(X)%*%y)
    plot(zoo::rollsum( LS,window.size), xlab='locus', ylab='beta',type='l',main='LS')
    lines(zoo::rollsum(beta_,window.size),col='red')
    cat('LS var:',round(var( zoo::rollsum(LS,window.size)-zoo::rollsum(beta_,window.size)),4),'\n')
}

###
# computing predictive power
#
###
MSE_XB_lasso <- MSE_XB_SS <- MSE_XB_ridge <- MSE_XB_npbayes <- rep(0,q)
MSE_beta_lasso <- MSE_beta_SS <- MSE_beta_ridge <- MSE_beta_npbayes <- rep(0,q)
for(i in 1:q){
    beta_lasso <- rep(0,length(beta_))
    beta_lasso[betas.lasso[[i]][,1]] <- betas.lasso[[i]][,2]
    MSE_XB_lasso[i] <- sqrt(mean((Xs[[i]]%*%betas[[i]]- Xs[[i]]%*%beta_lasso - intercept.lasso[i])^2))
    MSE_beta_lasso[i] <- sqrt(mean((betas[[i]]  -beta_lasso )^2))

    beta_SS <- rep(0,length(beta_))
    beta_SS[betas.SSlasso[[i]][,1]] <- betas.SSlasso[[i]][,2]
    MSE_XB_SS[i] <- sqrt(mean((Xs[[i]]%*%betas[[i]]- Xs[[i]]%*%beta_SS - intercept.sslasso[i])^2))
    MSE_beta_SS[i] <- sqrt(mean((betas[[i]]  -beta_SS )^2))


    MSE_XB_ridge[i] <- mean((Xs[[i]]%*%betas[[i]]- Xs[[i]]%*%betas.ridge[[i]] - intercept.ridge[i])^2)
    MSE_beta_ridge[i] <- sqrt(mean((betas[[i]]  -betas.ridge[[i]] )^2))

    MSE_XB_npbayes[i] <- sqrt(mean((Xs[[i]]%*%betas[[i]]- Xs[[i]]%*%  beta.npbayes[[i]] )^2))
    MSE_beta_npbayes[i] <- sqrt(mean((betas[[i]]  -beta.npbayes[[i]] )^2))

}
