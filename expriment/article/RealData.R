###
#  real data analysis
#
#
# D: 2018-01-08
###
rm(list=ls())
set.seed(1)
save.fig=F
graphics.off()
library(ggplot2)
library(NPBayes)
library(glmnet)
library(latex2exp)
L	<- 8
I	<- 2^L
burnin <-  50000
n.gibbs <-  500000
thin <- 3 # take ever x cm
source("misc/3-exp_matrix.R")

##
# setting up data
##
# data

d1 <- read.table("../../data/Zeng/bm4zb.out", skip = 9, row.names = 1)
d2 <- read.table("../../data/Zeng/bm6zb.out", row.names = 1)
d3 <- read.table("../../data/Zeng/bs4zb.out", skip = 9, row.names = 1)
d4 <- read.table("../../data/Zeng/bs6zb.out", row.names = 1)

# X
X <- as.matrix(rbind(d1[, 1:45], d2[, 1:45]))
#X <- as.matrix(rbind(d1[, 1:45], d2[, 1:45], d3[, 1:45], d4[, 1:45]))
X <- X[, -(1:6)] # odrzucam 1. chromosom X
X[X == 9] <- NA
X <- replaceNA(X) - 1
n <- nrow(X)
M <- ncol(X)
# Y
y <- as.matrix(c(d1[, 46], d2[, 46]), ncol = 1)
#y <- as.matrix(c(d1[, 46], d2[, 46], d3[, 46], d4[, 46]), ncol = 1)
y <- scale(y)
# map
markers <- c(6, 16, 23)[-1] # bez 1. chromosomu
chr <- length(markers)
len <- c(0, 3.60, 10.60, 9.20, 17.20, 18.70, 0, 6.98, 10.10, 4.94, 6.51, 6.19,
         20.46, 12.78, 3.90, 4.55, 7.49, 30.02, 16.85, 4.34, 3.71, 7.03, 0, 4.99,
         9.34, 6.97, 7.44, 14.46, 6.79, 3.55, 6.32, 11.86, 4.58, 6.85, 6.35,
         11.79, 12.88, 9.15, 3.30, 7.98, 13.09, 10.04, 3.70, 9.79, 3.43)[-(1:6)]
len_cum <- unlist(tapply(len, rep(1:chr, markers), cumsum), use.names = FALSE)
by <- 2


res <- pseudomarkers(X, markers, len_cum, by)
X <- res$X
M <- ncol(X)
markers <- res$markers
mark_1 <- seq(1,73,by=thin)
mark_2 <- seq(74,73+88,by=thin)
X <- X[,c(mark_1, mark_2)]
markers <- c(length(mark_1), length(mark_2))

q <- 1
findDuplicate2 <- function(X) {
  X[X < -0.7] <- -1
  X[X > 0.7] <- 1
  X[X > -0.3 & X < 0.3] <- 0
  dupl.ind <- which(duplicated(X, MARGIN = 2)) # ktore kolumny sie powtarzaja
  dupl.repl <- which(duplicated(X, MARGIN = 2, fromLast = TRUE)) # z tymi sie powtarzaja
  dupl <- cbind(dupl.ind, dupl.repl)
  return(dupl)
}
dupl <- findDuplicate2(X)


##
# estimating simple mode
##

X <- scale(X)
beta0 <- glmnet(X, y, alpha = 0, lambda = 0,standardize=F,intercept=F)$beta
res <- gibbs.normal(n.gibbs,
                    y,
                    X,
                    c(-0.4,0.4),
                    L,
                    as.vector(beta0),cpp=T)
par(mfrow=c(3,2))
plot(res$beta.gibbs[,which.max(as.vector(abs(beta0)))])
hist(res$beta.gibbs[,which.max(as.vector(abs(beta0)))])
plot(res$beta.gibbs[,4])
hist(res$beta.gibbs[,4])
plot(res$sigma.gibbs)
hist(res$sigma.gibbs)
par(mfrow=c(1,1))
plot(res$a.vec[-1],colMeans(exp(res$pi.gibbs[burnin:n.gibbs,])),type='l',xlab='beta',xlim=c(-0.2,0.2),ylab='pi(beta)')
EX = sum(res$a.vec[-1]*colMeans(exp(res$pi.gibbs[burnin:n.gibbs,])))
EX2 = sum(res$a.vec[-1]^2*colMeans(exp(res$pi.gibbs[burnin:n.gibbs,])))
EX3 = sum(res$a.vec[-1]^3*colMeans(exp(res$pi.gibbs[burnin:n.gibbs,])))
VX  = EX2 - EX^2
skewness = (EX3  - 3 * EX * VX - EX^3)/(VX^(3/2))
Kurt = sum((res$a.vec[-1]-EX)^4*colMeans(exp(res$pi.gibbs[burnin:n.gibbs,])))/VX^2
cat('E[beta] = ',round(EX,2),'\n')
cat('V[beta] = ',round(EX2-EX^2,2),'\n')


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
            CI[(1+marker.tot[i-1]):marker.tot[i],] = t(apply(CI.window, 2, quantile, probs=c(0.025,0.975)))
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
    res <- list(fig=fig )
    if(!is.null(x.sim)){
        fig <- fig  + geom_ribbon(data = data.frame(locus=1:length(x),
                                                    low = CI[,1]
                                                    ,upp = CI[,2]),
                                  aes(x=locus,ymin = low, ymax=upp), alpha=0.2, fill='blue')
        res <- list(fig=fig, res = beta_true.window-x.window )
    }

    return(res )
}



Betas.np <-  c(apply(res$beta.gibbs[burnin:n.gibbs,], 2, quantile, probs=0.5))


plot.graph2 <- function(x,
                       markers,
                       beta_true= NULL,
                       x.sim = NULL,
                       probs = c(0.025,0.975)){
    x.window <- rep(0,length(x))
    beta_true.window = rep(0, length(x))
    marker.tot <- c(1, cumsum(markers))
    CI <- matrix(0, ncol=2, nrow=length(x))

        if(!is.null(x.sim)){
            CI = t(apply(x.sim, 2, quantile, probs=c(probs[1],probs[2])))
        }

    df <- data.frame(locus = 1:length(x), beta = x)
    fig <- ggplot(df)
    fig <- fig + geom_line(aes(x=locus, y=beta),
                               colour ="blue",alpha=0.5)
    fig <- fig +  geom_vline(xintercept = marker.tot,
                             colour ="black",
                             linetype = "dashed",
                             alpha = 0.5)
    res <- list(fig=fig )
    if(!is.null(x.sim)){
        fig <- fig  + geom_ribbon(data = data.frame(locus=1:length(x),
                                                    low = CI[,1]
                                                    ,upp = CI[,2]),
                                  aes(x=locus,ymin = low, ymax=upp), alpha=0.2, fill='blue')
        res <- list(fig=fig, res = beta_true.window-x.window )
    }

    return(res )
}
CI_int = list()
CI_int[[1]] = c(0.25,0.75)
CI_int[[2]] = c(0.05,0.95)
fig.np <- list()
for(i in 1:2){
    res.np_CI  <- plot.graph2(Betas.np, markers, beta_true = NULL, x.sim =res$beta.gibbs[burnin:n.gibbs,] ,probs=CI_int[[i]])
    fig.np[[i]] <- res.np_CI$fig
    fig.np[[i]] <- fig.np[[i]] + ylab(TeX("$\\beta$"))
    fig.np[[i]] <- fig.np[[i]] + theme_bw()+ theme(axis.line = element_line(colour = "black"),
                         panel.border = element_blank(),
                         panel.background = element_blank(),
                         axis.text.x  = element_text(size=20),
                         axis.text.y  = element_text(size=20),
                         axis.title.y = element_text(angle=0,vjust=0.5,size=22,face="bold"),
                         axis.title.x = element_text(size=20,face="bold")) + ylim(c(-0.3,0.2))
}


###
#
# Run previous analysis
#
##

library(PolyMixed)
Forward_mix <- function (GeneMixObj, markers, dCM = 1, dupl = NULL, liberal = TRUE,
                         qval = NULL, alpha = 0.05, estParam = T)
{
    M <- ncol(GeneMixObj$LikObj$UX)
    find <- c()
    GeneMixObj <- mixedModel(GeneMixObj, estPar = estParam, dupl = dupl)
    if (GeneMixObj$tauOff == T)
        GeneMixObj$tau = 0
    t <- GeneMixObj$t
    maxv <- max(abs(t))
    if (is.null(qval))
        qval <- qnorm(1 - 0.05/(M * 2) * 5)
    while (maxv > qval) {
        find <- c(GeneMixObj$find, setdiff(which(abs(t) == maxv),
                                           dupl[, 1]))
        GeneMixObj <- mixedModel(GeneMixObj, find = find, dupl = dupl)
        t <- GeneMixObj$t
        if (!is.null(dupl)) {
            while (sum(abs(t[dupl[, 1]] - t[dupl[, 2]])) > 0) t[dupl[,
                                                                     1]] <- t[dupl[, 2]]
        }
        maxv <- max(abs(t[setdiff(1:M, c(find, dupl[, 1]))]),
                    na.rm = TRUE)
    }
    beta_critval = 0.02
    if (GeneMixObj$tauOff == FALSE) {
        beta_critval = NULL
    }
    crit <- calculateCrit(t, markers, d = dCM, beta = beta_critval,
                          alpha = alpha)
    find <- GeneMixObj$find
    return(GeneMixObj)
}
MixGeneObj <- SetupGeneMix('Y ~ 1',
                           data = data.frame(Y=y),
                           X=X,
                           meanMuOff = F,
                           tauOff = F)
MixGeneObj_forward <-Forward_mix(MixGeneObj,
                                               markers,
                                               dupl = dupl)
MixGeneObj <- mixedModelForwardBackward(MixGeneObj,
                                        markers,
                                        dupl = dupl)

pi.cdf      <- rbind(0,apply(exp(res$pi.gibbs[burnin:n.gibbs,]),1,cumsum))
pi.quant    <- apply(pi.cdf,1,quantile,c(0.025,0.50,0.975))


#   Pi estimation plot
mu <- MixGeneObj$beta[2]
tau <- MixGeneObj$tau
#phi.x <-
# if(save.fig)
#     pdf('CDF.pdf')
# par(mfrow=c(1,1))
# plot(res$a.vec,pi.quant[2,],col = 1,lwd = 2,xlim=c(-0.3,0.2), main = "",ylab = "CDF",xlab = expression(beta), type='l',
#      cex.lab = 1.5)
# lines(res$a.vec, pnorm(res$a.vec,mean=mu, sd = tau), col=2, lwd=2)
# abline(v=0,lw=2,col = adjustcolor("black", alpha = 0.3),lty=2)
# n_betas_forward <- length(MixGeneObj_forward$beta) - 2
#
# forward_q <- ((dim(X)[2]-n_betas_forward)/dim(X)[2])*pnorm(res$a.vec,mean=MixGeneObj_forward$beta[2], sd = MixGeneObj_forward$tau)
# for(i in 3:length(MixGeneObj_forward$beta) ){
#     forward_q <- forward_q + (1/dim(X)[2]) * (MixGeneObj_forward$beta[i] <res$a.vec)
# }
#
# lines(res$a.vec,forward_q, col=3, lwd=2)
# lines(res$a.vec,pi.quant[1,],col=1,lty = 2,lwd = 2)
# lines(res$a.vec,pi.quant[3,],col=1,lty = 2,lwd = 2)
# if(save.fig)
#     dev.off()

EX_f  = sum(res$a.vec[-1]*diff(forward_q))
EX2_f = sum(res$a.vec[-1]^2*diff(forward_q))
EX3_f = sum(res$a.vec[-1]^3*diff(forward_q))
VX_f  = EX2_f - EX_f^2
skewness_f = (EX3_f  - 3 * EX_f * VX_f - EX_f^3)/(VX_f^(3/2))
Kurt_f = sum((res$a.vec[-1]-EX)^4*colMeans(exp(res$pi.gibbs[burnin:n.gibbs,])))/VX^2

beta_forward <- rep(0, dim(X)[2])+MixGeneObj_forward$beta[2]
beta_forward[ MixGeneObj_forward$find] <-beta_forward[ MixGeneObj_forward$find] + MixGeneObj_forward$beta[-(1:2)]



#  (y- X\beta)(y-X\beta)/sigma^2 \beta^T\beta
# Q = X'X/sigma^2 + I \tau^2
# E[\beta|X] = Q * X/sigma^2 \tilde y
Q = (t(X)%*%X)/MixGeneObj$sigma^2
diag(Q) <- diag(Q)  + 1/MixGeneObj$tau^2
Ebeta = solve(Q, t(X)%*%(y - MixGeneObj$beta[1])/MixGeneObj$sigma^2 + MixGeneObj$beta[2]/ MixGeneObj$tau^2 )
Sigma <- solve(Q)
for(i in 1:2){
    fig.np[[i]] <- fig.np[[i]] +   geom_ribbon(data = data.frame(locus=1:dim(X)[2],
                                                      low = qnorm(CI_int[[i]][1], mean=Ebeta, sd=sqrt(diag(Sigma))),
                                                      upp = qnorm(CI_int[[i]][2], mean=Ebeta, sd=sqrt(diag(Sigma)))),
                                     aes(x=locus,ymin = low, ymax=upp), alpha=0.1, fill='red')
fig.np[[i]] <- fig.np[[i]] +   geom_line(aes(x=1:dim(X)[2], y=Ebeta),
                                          colour ="red",alpha=0.5)
}

if(save.fig){
    for(i in 1:2)
        ggsave(paste('CV_',i,'.pdf',sep=""), fig.np[[i]], height=6, width=10)
}

# Create a data frame for plotting
df <- data.frame(
    x = res$a.vec,
    y1 = pi.quant[2,],
    y2 = pnorm(res$a.vec, mean = mu , sd =tau),
    y3 = pi.quant[1,],
    y4 = pi.quant[3,]
)

# Plot using ggplot2
fig_cdf <- ggplot(df, aes(x = x)) +
    geom_line(aes(y = y1), col = "blue", lwd = 1) +
    geom_line(aes(y = y2), col = "red", lwd = 1) +
    geom_line(aes(y = y3), col = "blue", lwd = 1, linetype=2) +
    geom_line(aes(y = y4), col = "blue", lwd = 1, linetype=2) +
    labs(
        title = "",
        y = "CDF",
        x = expression(beta)
    ) +
    xlim(c(-0.3, 0.2)) +
    theme_minimal() +  # Adjust the theme as needed
     theme(axis.text.x  = element_text(size=20),
           axis.title.x = element_text(size=22,face="bold"),
            axis.text.y  = element_text(size=20),
            axis.title.y = element_text(size=22,face="bold")) +
    geom_vline(xintercept = 0, col = "black", linetype = 2, size = 1, alpha = 0.3)
print(fig_cdf)
if(save.fig){
        ggsave("CDF.pdf",fig_cdf, height=6, width=10)
}

