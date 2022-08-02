##
# power simulation for mixed effect model
# to create Table 3: set.seed(2) q=1000, table =2 ,gives sparse mixed effect
# to create Table 4: set.seed(2) q=1000, table =1 ,gives laplace distribution
# D: 2022-02-15
##

rm(list=ls())
save.fig=F
graphics.off()
library(doParallel)
library(purrr)
library(bigstep)
library(PolyMixed)
library(RcppEigen)
library(ggplot2)
library(bigstep)
set.seed(12)

n.gibbs <- 10000 #samples for MCMC
n_cores = 8
distribution = c(2)  # 1 - gaussian random effect, 2- Laplace randeom effect, 3 -sparse
shift <- -2. #-0.02 # shift in the laplace distribution
nu <- 0.75 #assymetry in Laplace distribution
q     <- 100   # number of simulations
n_cores <- min(q, n_cores)
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
tau <- 0.005              # sd for polygenic effect
mu <- -rep(0.01, q)      # mean for polygenic effect



nc = NULL
qtl.pos <- ceiling(M * c(0.2,0.5,0.8))
mag.pos = 0.1
mag.neg = 0.2
if(signal==1){
    mag = 0.5
}else if(signal==-1){
    mag = 0
}
qtl  <- c(mag.pos, qtl.pos[1], 1) # qtl (beta; position: marker, trait or traits)
qtl2 <- c(mag.pos, qtl.pos[2], 1)
qtl3 <- c(-mag.neg, qtl.pos[3], 1)
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
            V <- rgamma(dim(X)[2],nu,nu)
            beta_ <- mu[1] + tau * (shift * (V-1) +  sqrt(V)*Z)
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

    cl <- makeCluster(n_cores, outfile="")
    registerDoParallel(cl)
    i0 = 1
    par.out <- foreach(i = i0:q)%dopar%{

        library(glmnet)
        library(SSLASSO)
        library(NPBayes)
        library(PolyMixed)
        cat('i=',i,'\n')
        set.seed(i)
        #unadjusted lasso
        cvfit <- cv.glmnet(Xs[[i]], ys[[i]])
        coef.lasso <- coef(cvfit, s = "lambda.1se")
        index.lasso <- which(coef.lasso!=0)
        index.lasso <- index.lasso[-1]
        lasso.beta <- rep(0,dim(X)[2])
        lasso.beta[index.lasso-1] <- coef.lasso[index.lasso]
        betas.lasso_unadj <- lasso.beta
        intercept.lasso_unadj <-coef.lasso[1]


        # Lasso
        cvfit <- cv.glmnet(cbind(rowSums(Xs[[i]]),Xs[[i]]), ys[[i]])
        coef.lasso <- coef(cvfit, s = "lambda.1se")
        index.lasso <- which(coef.lasso!=0)
        index.lasso <- index.lasso[-(1:2)]
        lasso.beta <- rep( coef.lasso[2],dim(X)[2])
        lasso.beta[index.lasso-2] <- lasso.beta[index.lasso-2]  + coef.lasso[index.lasso]
        betas.lasso <- lasso.beta
        intercept.lasso <-coef.lasso[1]

        #Spike and slab lasso
        sslasso.fit <- SSLASSO(cbind(rowSums(Xs[[i]]),Xs[[i]]), ys[[i]], variance = "unknown")
        index <- setdiff(sslasso.fit$model,1)
        sslasso.beta <- rep( sslasso.fit$beta[1, ncol(sslasso.fit$beta)], dim(X)[2])
        sslasso.beta[index-1] <- sslasso.beta[index-1] + sslasso.fit$beta[setdiff(sslasso.fit$model,1),ncol(sslasso.fit$beta)]
        betas.SSlasso <- sslasso.beta
        intercept.sslasso <- sslasso.fit$intercept[1,dim(sslasso.fit$intercept)[2]]

        sslasso.fit <- SSLASSO(Xs[[i]], ys[[i]], variance = "unknown")
        index <- sslasso.fit$model
        sslasso.beta <- rep( 0, dim(X)[2])
        sslasso.beta[index-1] <- sslasso.fit$beta[setdiff(sslasso.fit$model,1),ncol(sslasso.fit$beta)]
        betas.SSlasso_unadj <- sslasso.beta
        intercept.sslasso_unadj <- sslasso.fit$intercept[1,dim(sslasso.fit$intercept)[2]]


        #ridge regression
        cvfit <- cv.glmnet(cbind(rowSums(Xs[[i]]),Xs[[i]]), ys[[i]], alpha=0)
        coef.ridge <- coef(cvfit, s = "lambda.1se")
        betas.ridge <- coef.ridge[-(1:2)] + coef.ridge[2]
        intercept.ridge <- coef.ridge[1]


        cvfit <- cv.glmnet(Xs[[i]], ys[[i]], alpha=0)
        coef.ridge <- coef(cvfit, s = "lambda.1se")
        betas.ridge_unadj <- coef.ridge[-1]
        intercept.ridge_unadj <- coef.ridge[1]
        #hBayes
        res <- gibbs.normal(n.gibbs,
                            ys[[i]],
                            Xs[[i]],
                            c(-0.5,0.5),
                            8,
                            betas.ridge_unadj)
        #res <- gibbs.normal.fixed.sigma(n.gibbs,
        #                    ys[[i]],
        #                    Xs[[i]],
        #                    sigma,
        #                    c(-0.5,0.5),
        #                    8,
        #                    sslasso.fit$beta[,ncol(sslasso.fit$beta)])

        beta.npbayes <-  base::colMeans(res$beta.gibbs[floor(n.gibbs/3):n.gibbs,])



        #forward backward
        MixGeneObj <- SetupGeneMix('Y ~ 1',
                                   data = data.frame(Y=ys[[i]]),
                                   X=Xs[[i]],
                                   meanMuOff = F,
                                   tauOff = F)

        MixGeneObj <- mixedModelForwardBackward(MixGeneObj,
                                                markers,
                                                dupl  = findDuplicate(Xs[[i]]))
        #compute posterior estimate of beta
        Q_2 = (t(Xs[[i]])%*%Xs[[i]])/MixGeneObj$sigma^2
        diag(Q_2) <- diag(Q_2)  + 1/MixGeneObj$tau^2
        Ebeta_2 = solve(Q_2, (t(Xs[[i]])%*%(ys[[i]] - MixGeneObj$beta[1] - Xs[[i]][,MixGeneObj$find,drop=F]%*%MixGeneObj$beta[-c(1:2)]))/MixGeneObj$sigma^2 + MixGeneObj$beta[2]/ MixGeneObj$tau^2 )
        Ebeta_2[MixGeneObj$find] =Ebeta_2[MixGeneObj$find] +MixGeneObj$beta[-c(1:2)]

        beta.forwardback <- Ebeta_2
        intercept.forwardback <-  MixGeneObj$beta[1]


        list( betas.lasso = betas.lasso,
              intercept.lasso =intercept.lasso,

              betas.lasso_unadj = betas.lasso_unadj,
              intercept.lasso_unadj = intercept.lasso_unadj,

              betas.SSlasso = betas.SSlasso,
              intercept.sslasso = intercept.sslasso,

              betas.SSlasso_unadj = betas.SSlasso_unadj,
              intercept.sslasso_unadj = intercept.sslasso_unadj,

              betas.ridge = betas.ridge,
              intercept.ridge = intercept.ridge,


              betas.ridge_unadj = betas.ridge_unadj,
              intercept.ridge_unadj = intercept.ridge_unadj,

              beta.npbayes = beta.npbayes,

              beta.forwardback=  beta.forwardback,
              intercept.forwardback = intercept.forwardback
              )
    }
    list2env(purrr::transpose(par.out),globalenv())
    stopCluster(cl)
    save( betas.lasso,intercept.lasso, betas.SSlasso, intercept.sslasso, betas.ridge, intercept.ridge, beta.npbayes,beta.forwardback,intercept.forwardback, file = "simulation2.RData")
}



library(zoo)
library(ggplot2)
window.size= 10
#' @param x.sim simulations to get posterior conf int
plot.graph <- function(x,
                       window.size,
                       markers,
                       beta_true= NULL,
                       name='',
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
    fig <- fig + ylab(TeX("$\\beta$"))+ theme(axis.title.y = element_text(angle=0,vjust=0.5,size=20,face="bold"),
                                              axis.title.x = element_text(size=16,face="bold"))
    if(!is.null(x.sim)){
        fig <- fig  + geom_ribbon(data = data.frame(locus=1:length(x),
                                                          low = CI[,1]
                                                          ,upp = CI[,2]),
                                        aes(x=locus,ymin = low, ymax=upp), alpha=0.2, fill='blue')
    }
    fig <- fig +  labs(title=paste(name))
    #fig <- fig +  labs(title=paste(name,', RMSE = ',round(sqrt(mean((beta_true.window-x.window)^2))/sd(beta_true.window),3) , ", window size = ", window.size,sep=""))
    fig <- fig+ theme(plot.title = element_text(hjust = 0.5, size=20))

    return(list(fig=fig, res = beta_true.window-x.window, beta_true = beta_true.window ) )
}
i=1
library(latex2exp)
res.lasso <- plot.graph(betas.lasso[[q]],window.size, markers, beta_true = beta_,name="lasso")
fig.lasso <- res.lasso$fig

print(fig.lasso)
if(save.fig)
    ggsave('lasso.jpeg',fig.lasso)

res.lasso_un <- plot.graph(betas.lasso_unadj[[q]],window.size, markers, beta_true = beta_,name="lasso unadj")
fig.lasso_un <- res.lasso_un$fig

print(fig.lasso_un)
if(save.fig)
    ggsave('lasso_un.jpeg',fig.lasso_un)


res.ridge <- plot.graph(betas.ridge[[q]],window.size, markers, beta_true = beta_,name="ridge adj")
fig.ridge <- res.ridge$fig

print(fig.ridge)
if(save.fig)
    ggsave('ridge.jpeg',fig.ridge)

res.ridge <- plot.graph(betas.ridge_unadj[[q]],window.size, markers, beta_true = beta_,name="ridge")
fig.ridge_unadj <- res.ridge$fig

print(fig.ridge_unadj)
if(save.fig)
    ggsave('ridge_un.jpeg',fig.ridge_unadj)

res.sslasso <- plot.graph(betas.SSlasso[[q]],window.size, markers, beta_true = beta_,name="SSlasso")
fig.sslasso <- res.sslasso$fig

print(fig.sslasso)
if(save.fig)
    ggsave('sslasso.jpeg',fig.sslasso)


res.sslasso <- plot.graph(betas.SSlasso_unadj[[q]],window.size, markers, beta_true = beta_,name="SSlasso unadj")
fig.sslasso_unadj <- res.sslasso$fig

print(fig.sslasso_unadj)
if(save.fig)
    ggsave('sslasso_unadj.jpeg',fig.sslasso_unadj)

res.np <- plot.graph(beta.npbayes[[q]],window.size, markers, beta_true = beta_ ,name="Polya prior")
fig.np <- res.np$fig

print(fig.np)
if(save.fig)
    ggsave('bayes.jpeg',fig.np)

library(latex2exp)
res.fb <- plot.graph(beta.forwardback[[q]],window.size, markers, beta_true = beta_,name="mixed effect")
fig.fb <- res.fb$fig

print(fig.fb)
if(save.fig)
    ggsave('fb.jpeg',fig.fb)




#geom_ribbon(aes(ymin = Anomaly10y - Unc10y, ymax = Anomaly10y + Unc10y), alpha = 0.2)
#cat('lasso var:',round(var(zoo::rollsum(coef.lasso[-1],window.size )-zoo::rollsum(beta_,window.size)),4),'\n')
#cat('nbayes var:',round(var( zoo::rollsum(colMeans(res$beta.gibbs[floor(n.gibbs/3):n.gibbs,]),window.size)-zoo::rollsum(beta_,window.size)),4),'\n')





###
# computing predictive power
#
###
MSE_XB_lasso <- MSE_XB_SS <- MSE_XB_ridge <- MSE_XB_fb <- MSE_XB_npbayes <- rep(0,q)
MSE_beta_lasso <- MSE_beta_SS <- MSE_beta_ridge <-MSE_beta_fb <- MSE_beta_npbayes <- rep(0,q)
for(i in 1:q){
    XB <- Xs[[i]]%*%betas[[i]]
    MSE_XB_lasso[i] <- sqrt(mean((XB- Xs[[i]]%*%betas.lasso[[i]] - intercept.lasso[[i]])^2))/sd(XB)
    MSE_beta_lasso[i] <- sqrt(mean((betas[[i]]  -betas.lasso[[i]] )^2))/sd(betas[[i]])

    MSE_XB_SS[i] <- sqrt(mean((XB- Xs[[i]]%*%betas.SSlasso[[i]] - intercept.sslasso[[i]])^2))/sd(XB)
    MSE_beta_SS[i] <- sqrt(mean((betas[[i]]  -betas.SSlasso[[i]]  )^2))/sd(betas[[i]])


    MSE_XB_ridge[i] <- mean((XB- Xs[[i]]%*%betas.ridge_unadj[[i]] - intercept.ridge_unadj[[i]])^2)/sd(XB)
    MSE_beta_ridge[i] <- sqrt(mean((betas[[i]]  -betas.ridge_unadj[[i]] )^2))/sd(betas[[i]])

    MSE_XB_npbayes[i] <- sqrt(mean((XB- Xs[[i]]%*%  beta.npbayes[[i]] )^2))/sd(XB)
    MSE_beta_npbayes[i] <- sqrt(mean((betas[[i]]  -beta.npbayes[[i]] )^2))/sd(betas[[i]])


    MSE_XB_fb[i] <- sqrt(mean((XB - Xs[[i]]%*%beta.forwardback[[i]] - intercept.forwardback[[i]])^2))/sd(XB)
    MSE_beta_fb[i] <- sqrt(mean((betas[[i]]  -beta.forwardback[[i]] )^2))/sd(betas[[i]])


}
Table <- rbind(c(mean(MSE_XB_lasso),mean(MSE_XB_SS),mean(MSE_XB_ridge),mean(MSE_XB_npbayes),mean(MSE_XB_fb)),
               c(mean(MSE_beta_lasso),mean(MSE_beta_SS),mean(MSE_beta_ridge),mean(MSE_beta_npbayes),mean(MSE_beta_fb)))
rownames(Table) <- c("MSE of $ {\\bf X} \\beta $","MSE of $ \\beta $")
colnames(Table) <- c("lasso", "sslasso", "ridge", "non parameteric","forward backward")
library(xtable)
print.xtable(xtable(Table, digits=3),type="latex",sanitize.text.function = function(x) x)
#plot true density
#n=10^4
#V <- rgamma(n=n,0.75,0.75)
#Z <- rnorm(n=n)
#B <- mu[1] + tau * (shift * (V-1) +  sqrt(V)*Z)

