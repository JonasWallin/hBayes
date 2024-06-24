##
# power simulation for mixed effect model
# D: 2022-02-15
# D: 2023-11-10 adding methods from https://arxiv.org/pdf/2208.10910.pdf
#
##
#devtools::install_github("stephenslab/susieR")
#devtools::install_github("stephenslab/mr.ash.alpha")
#rm(list=ls())
plot.fig = F
save.fig = T
save.data = T
graphics.off()
source("qtl_util.R")
library(SSLASSO)
library(glmnet)
library(latex2exp)
library(zoo)
library(doParallel)
library(purrr)
library(bigstep)
library(PolyMixed) #devtools::install_github('JonasWallin/PolyMixed')
library(RcppEigen)
library(ggplot2)
library(bigstep)
library(mr.ash.alpha) #devtools::install_github("stephenslab/mr.ash.alpha")
library(susieR) #devtools::install_github("stephenslab/susieR")
library(EMVS) #devtools::install_version("EMVS", repos = "http://cran.us.r-project.org")
library(varbvs)
library(horseshoe) #for horseshoe
set.seed(12)

if(plot.fig){
    q     <- 1  #we only generate a single sample to plot data.
    n_cores = 1
}else{
    q     <- 200
    n_cores = detectCores() - 6
}
distribution = c(2)  # 1 - gaussian random effect, 2- Laplace randeom effect, 3 -sparse
shift <- -2. #-0.02 # shift in the laplace distribution
nu <- 0.75 #assymetry in Laplace distribution

n_cores <- min(q, n_cores)
#thin =  how often to observe X (i.e. cm distnace)
#type =  1 - pure measurement error (\sigma^2I )
#        2 -  mixed effect     (\sigma^2I + XX^T \tau )


thin    = 1
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
mag.pos = 0.2
mag.neg = 0.2
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
        library(mr.ash.alpha)
        library(glmnet)
        library(susieR)
        library(varbvs)
        library(EMVS)
        library(SSLASSO)
        library(horseshoe)
        # Lasso
        cvfit <- cv.glmnet(cbind(rowSums(Xs[[i]]),Xs[[i]]), ys[[i]])
        coef.lasso <- coef(cvfit, s = "lambda.1se")
        index.lasso <- which(coef.lasso!=0)
        index.lasso <- index.lasso[-(1:2)]
        lasso.beta <- rep( coef.lasso[2],dim(X)[2])
        lasso.beta[index.lasso-2] <- lasso.beta[index.lasso-2]  + coef.lasso[index.lasso]
        betas.lasso <- lasso.beta
        intercept.lasso <-coef.lasso[1]

        #Mr.ASH from https://arxiv.org/pdf/2208.10910.pdf
        fit.mr.ash  = mr.ash(Xs[[i]],ys[[i]], method_q = "sigma_dep_q",max.iter=2000,standardize=T)
        betahat.ash = predict(fit.mr.ash, type = "coefficients")


        #taking options from https://github.com/stephenslab/mr-ash-workflow/blob/master/code/method_wrapper.R
        fit.susie        <- susie(X = Xs[[i]], y = ys[[i]], standardize = T, L = 20)
        #sanity: predict(fit.susie, Xs[[i]])-Xs[[i]]%*%coef(fit.susie)[-1]-fit.susie$intercept
        fit.varbvs <- varbvs(Xs[[i]], Z = NULL, ys[[i]], verbose = FALSE, logodds = seq(-log10(dim(Xs[[i]])[2]),1,length.out = 40))
        beta.varbvs  = c(rowMeans(fit.varbvs$alpha * fit.varbvs$mu))

        fit.varbvsmix <- varbvsmix(Xs[[i]], Z = NULL, ys[[i]])
        beta.varbvsmix  = c(rowSums(fit.varbvsmix$alpha * fit.varbvsmix$mu))
        #EMVS

        v0 = seq(0.1, 20, length.out = 50)
        v1 = 1000
        beta_init = rep(1, dim(Xs[[i]])[2])
        sigma_init = 1
        a = b = 1
        epsilon = 10^{-5}

        result.EMVS = EMVS(ys[[i]], Xs[[i]], v0 = v0, v1 = v1, type = "betabinomial",
                      independent = FALSE, beta_init = beta_init, sigma_init = sigma_init,
                      epsilon = epsilon, a = a, b = b)

        #
        sslasso.fit <- SSLASSO(Xs[[i]], ys[[i]], variance = "unknown")
        index <- sslasso.fit$model
        sslasso.beta <- rep( 0, dim(X)[2])
        sslasso.beta[index-1] <- sslasso.fit$beta[setdiff(sslasso.fit$model,1),ncol(sslasso.fit$beta)]
        betas.SSlasso <- sslasso.beta
        intercept.sslasso <- sslasso.fit$intercept[1,dim(sslasso.fit$intercept)[2]]

        hs.object <- horseshoe(ys[[i]],
                               cbind(1,Xs[[i]]),
                               method.tau = "truncatedCauchy",
                               method.sigma ="Jeffreys",
                               nmc = 6000)

        list(
             betas.lasso = betas.lasso,
             intercept.lasso =intercept.lasso,
              beta.ash = betahat.ash[-1],
              intercept.ash =  betahat.ash[1],
             intercept.susie = fit.susie$intercept,
             beta.susie = coef(fit.susie)[-1],
             beta.varbvs = beta.varbvs,
             pred.varbvs = predict(fit.varbvs,Xs[[i]]),
             beta.varbvsmix = beta.varbvsmix,
             pred.varbvsmix = predict(fit.varbvsmix,Xs[[i]]),
             beta.EMVS = result.EMVS$betas[which(result.EMVS$log_g_function==EMVSbest(result.EMVS)$log_g_function)[1],],
             beta.EMVS2 = result.EMVS$betas[1,],
             betas.SSlasso = betas.SSlasso,
             intercept.SSlasso = intercept.sslasso,
             intercept.horseshoe  = hs.object$BetaHat[1],
             betas.horseshoe  = hs.object$BetaHat[-1]
              )
    }
    list2env(purrr::transpose(par.out),globalenv())
    stopCluster(cl)
}


    if(plot.fig){
        window.size= 10
        res.ash_un <- plot.graph(beta.ash[[q]], window.size, markers, beta_true = beta_,name="Mr.ash")
        res.susie <- plot.graph(beta.susie[[q]], window.size, markers, beta_true = beta_,name="Susie")
        res.varbvs <- plot.graph(beta.varbvs[[q]], window.size, markers, beta_true = beta_,name="Varbv")
        res.varbvsmix <- plot.graph(beta.varbvsmix[[q]], window.size, markers, beta_true = beta_,name="Varbvmix")
        res.EMVS <- plot.graph(beta.EMVS[[q]], window.size, markers, beta_true = beta_,name="EMVS")
        res.EMVS2 <- plot.graph(beta.EMVS2[[q]], window.size, markers, beta_true = beta_,name="EMVS 2")
        #res.sslasso <- plot.graph(betas.SSlasso[[q]], window.size, markers, beta_true = beta_,name="SSLASSO")

        if(save.fig){
            ggsave('ash.jpeg',res.ash_un$fig,width = 8, height = 6)
            ggsave('susie.jpeg',res.susie$fig,width = 8, height = 6)
            ggsave('varbvs.jpeg',res.varbvs$fig,width = 8, height = 6)
            ggsave('varbvsmix.jpeg',res.varbvsmix$fig,width = 8, height = 6)
            #ggsave('sslasso.jpeg',res.sslasso$fig)
            ggsave('EMVS.jpeg',res.EMVS$fig,width = 8, height = 6)
            ggsave('EMVS2.jpeg',res.EMVS2$fig,width = 8, height = 6)
        }
    }


    ###
    # computing predictive power
    #
    ###
    MSE_XB_HS<- MSE_XB_EMVS <- MSE_XB_lasso <- MSE_XB_SSlasso <- MSE_XB_ash  <- MSE_XB_susie<- MSE_XB_varbvs<- MSE_XB_varbvsmix<- rep(0,q)
    MSE_beta_HS <- MSE_beta_EMVS <- MSE_beta_lasso <- MSE_beta_SSlasso<- MSE_beta_ash <- MSE_beta_susie <- MSE_beta_varbvs <- MSE_beta_varbvsmix <- rep(0,q)
    beta.qtl.HS <- beta.qtl.EMVS <- beta.qtl.lasso <- beta.qtl.SSlasso <- beta.qtl.ash <- beta.qtl.susie <- beta.qtl.varbvs <- beta.qtl.varbvsmix <- matrix(0,nrow=q, ncol = length(qtl.pos))
    beta.smooth.qtl.HS <- beta.smooth.qtl.EMVS <- beta.smooth.qtl.true <-  beta.smooth.qtl.lasso <- beta.smooth.qtl.SSlasso <- beta.smooth.qtl.ash <- beta.smooth.qtl.varbvs <- beta.smooth.qtl.varbvsmix<- beta.smooth.qtl.susie <- matrix(0,nrow=q, ncol = length(qtl.pos))
    MSE_beta_smooth.HS <-MSE_beta_smooth.EMVS <- MSE_beta_smooth.lasso <- MSE_beta_smooth.SSlasso <- MSE_beta_smooth.ash <- MSE_beta_smooth.susie <- MSE_beta_smooth.varbvs<- MSE_beta_smooth.varbvsmix<- rep(0,q)

    for(i in 1:q){
        XB <- Xs[[i]]%*%betas[[i]]
        res.true <- smooth.beta(betas[[i]], betas[[i]], markers, window.size=10)
        beta.smooth.qtl.true[i, ] <- res.true$x[qtl.pos]


        Hs_res <- scores(XB, Xs[[i]], betas.horseshoe[[i]], intercept.horseshoe[[i]], betas[[i]])
        MSE_XB_HS[i] <- Hs_res[1]
        MSE_beta_HS[i] <- Hs_res[2]
        MSE_beta_smooth.HS[i] <- Hs_res[3]
        beta.qtl.HS[i, ] <- Hs_res[4:(3+length(qtl.pos))]
        beta.smooth.qtl.HS[i,] <- Hs_res[(4+length(qtl.pos)):(3+2*length(qtl.pos))]

        EMVS_res <- scores(XB, Xs[[i]], beta.EMVS[[i]], 0, betas[[i]])
        MSE_XB_EMVS[i] <- EMVS_res[1]
        MSE_beta_EMVS[i] <- EMVS_res[2]
        MSE_beta_smooth.EMVS[i] <- EMVS_res[3]
        beta.qtl.EMVS[i, ] <- EMVS_res[4:(3+length(qtl.pos))]
        beta.smooth.qtl.EMVS[i,] <- EMVS_res[(4+length(qtl.pos)):(3+2*length(qtl.pos))]


        lasso_res <- scores(XB, Xs[[i]], betas.lasso[[i]], intercept.lasso[[i]], betas[[i]])

        MSE_XB_lasso[i] <- lasso_res[1]
        MSE_beta_lasso[i] <- lasso_res[2]
        MSE_beta_smooth.lasso[i] <- lasso_res[3]
        beta.qtl.lasso[i, ] <- lasso_res[4:(3+length(qtl.pos))]
        beta.smooth.qtl.lasso[i,] <- lasso_res[(4+length(qtl.pos)):(3+2*length(qtl.pos))]


        lassoSS_res <- scores(XB, Xs[[i]], betas.SSlasso[[i]], intercept.SSlasso[[i]], betas[[i]])

        MSE_XB_SSlasso[i] <- lassoSS_res[1]
        MSE_beta_SSlasso[i] <- lassoSS_res[2]
        MSE_beta_smooth.SSlasso[i] <- lassoSS_res[3]
        beta.qtl.SSlasso[i, ] <- lassoSS_res[4:(3+length(qtl.pos))]
        beta.smooth.qtl.SSlasso[i,] <- lassoSS_res[(4+length(qtl.pos)):(3+2*length(qtl.pos))]

        ash_res <- scores(XB, Xs[[i]], beta.ash[[i]], intercept.ash[[i]], betas[[i]])

        MSE_XB_ash[i] <- ash_res[1]
        MSE_beta_ash[i] <- ash_res[2]
        MSE_beta_smooth.ash[i] <- ash_res[3]
        beta.qtl.ash[i, ] <- ash_res[4:(3+length(qtl.pos))]
        beta.smooth.qtl.ash[i,] <- ash_res[(4+length(qtl.pos)):(3+2*length(qtl.pos))]


        susie_res <- scores(XB, Xs[[i]], beta.susie[[i]], intercept.susie[[i]], betas[[i]])

        MSE_XB_susie[i] <- susie_res[1]
        MSE_beta_susie[i] <- susie_res[2]
        MSE_beta_smooth.susie[i] <- susie_res[3]
        beta.qtl.susie[i, ] <- susie_res[4:(3+length(qtl.pos))]
        beta.smooth.qtl.susie[i,] <- susie_res[(4+length(qtl.pos)):(3+2*length(qtl.pos))]

        varbvs_res <- scores(XB, Xs[[i]], beta.varbvs[[i]], NULL, betas[[i]], pred.varbvs[[i]])

        MSE_XB_varbvs[i] <- varbvs_res[1]
        MSE_beta_varbvs[i] <- varbvs_res[2]
        MSE_beta_smooth.varbvs[i] <- varbvs_res[3]
        beta.qtl.varbvs[i, ] <- varbvs_res[4:(3+length(qtl.pos))]
        beta.smooth.qtl.varbvs[i,] <- varbvs_res[(4+length(qtl.pos)):(3+2*length(qtl.pos))]


        varbvsmix_res <- scores(XB, Xs[[i]], beta.varbvsmix[[i]], NULL, betas[[i]], pred.varbvsmix[[i]])

        MSE_XB_varbvsmix[i] <- varbvsmix_res[1]
        MSE_beta_varbvsmix[i] <- varbvsmix_res[2]
        MSE_beta_smooth.varbvsmix[i] <- varbvsmix_res[3]
        beta.qtl.varbvsmix[i, ] <- varbvsmix_res[4:(3+length(qtl.pos))]
        beta.smooth.qtl.varbvsmix[i,] <- varbvsmix_res[(4+length(qtl.pos)):(3+2*length(qtl.pos))]



    }

    Table <- cbind(c(mean(MSE_XB_lasso),mean(MSE_beta_lasso),mean(MSE_beta_smooth.lasso)),
                   c(mean(MSE_XB_ash),mean(MSE_beta_ash),mean(MSE_beta_smooth.ash)),
                   c(mean(MSE_XB_susie),mean(MSE_beta_susie),mean(MSE_beta_smooth.susie)),
                   c(mean(MSE_XB_varbvs),mean(MSE_beta_varbvs),mean(MSE_beta_smooth.varbvs)),
                   c(mean(MSE_XB_varbvsmix),mean(MSE_beta_varbvsmix),mean(MSE_beta_smooth.varbvsmix)),
                   c(mean(MSE_XB_EMVS),mean(MSE_beta_EMVS),mean(MSE_beta_smooth.EMVS)),
                   c(mean(MSE_XB_SSlasso),mean(MSE_beta_SSlasso),mean(MSE_beta_smooth.SSlasso)),
                   c(mean(MSE_XB_HS),mean(MSE_beta_HS),mean(MSE_beta_smooth.HS))
                   )


                   rownames(Table) <- c("RMSE of $ {\\bf X} \\beta $","RMSE of $ \\beta $","RMSE of $ \\tilde{\\beta} $")
    colnames(Table) <- c("lasso","mr.ash","susie","varbvs","varbvsmix","EMVS","SSlasso","horseshoe")

    library(xtable)
    print.xtable(xtable(Table, digits=3),type="latex",sanitize.text.function = function(x) x)

    sqrt(rowMeans(apply(beta.qtl.lasso, 1,function(x){(x-betas[[1]][qtl.pos])^2})))
    sqrt(rowMeans(apply(beta.qtl.ash, 1,function(x){(x-betas[[1]][qtl.pos])^2})))

    sqrt(rowMeans(apply(beta.qtl.EMVS, 1,function(x){(x-betas[[1]][qtl.pos])^2})))
    sqrt(rowMeans(apply(beta.smooth.qtl.EMVS, 1,function(x){(x-betas[[1]][qtl.pos])^2})))
    sqrt(rowMeans(apply(beta.qtl.SSlasso, 1,function(x){(x-betas[[1]][qtl.pos])^2})))
    if(save.data){
    saveRDS(list(
                 MSE_beta_HS               = MSE_beta_HS,
                 MSE_XB_HS                 = MSE_XB_HS,
                 MSE_beta_smooth_HS        = MSE_beta_smooth.HS,
                 MSE_beta_ash              = MSE_beta_ash,
                 MSE_XB_ash                = MSE_XB_ash,
                 MSE_beta_smooth_ash       = MSE_beta_smooth.ash,
                 MSE_beta_susie            = MSE_beta_susie,
                 MSE_XB_susie.             = MSE_XB_susie,
                 MSE_beta_smooth_susie     = MSE_beta_smooth.susie,
                 MSE_beta_varbvs           = MSE_beta_varbvs,
                 MSE_XB_varbvs             = MSE_XB_varbvs,
                 MSE_beta_smooth_varbvs    = MSE_beta_smooth.varbvs,
                 MSE_beta_varbvsmix        = MSE_beta_varbvsmix,
                 MSE_XB_varbvsmix          = MSE_XB_varbvsmix,
                 MSE_beta_smooth_varbvsmix = MSE_beta_smooth.varbvsmix,
                 MSE_beta_EMVS             = MSE_beta_EMVS,
                 MSE_XB_EMVS               = MSE_XB_EMVS,
                 MSE_beta_smooth_EMVS      = MSE_beta_smooth.EMVS)
            ,"MSE_sim2.rds")}
