###
#  real data analysis
#
#
# D: 2018-01-08
###
rm(list=ls())
graphics.off()
library(ggplot2)
library(NPBayes)
library(glmnet)
L	<- 8
I	<- 2^L
burnin <- 2000
n.gibbs <- 4000
source("3-exp_matrix.R")

##
# setting up data
##
# data

d1 <- read.table("../data/Zeng/bm4zb.out", skip = 9, row.names = 1)
d2 <- read.table("../data/Zeng/bm6zb.out", row.names = 1)
d3 <- read.table("../data/Zeng/bs4zb.out", skip = 9, row.names = 1)
d4 <- read.table("../data/Zeng/bs6zb.out", row.names = 1)

# X
#X <- as.matrix(rbind(d1[, 1:45], d2[, 1:45]))
X <- as.matrix(rbind(d1[, 1:45], d2[, 1:45], d3[, 1:45], d4[, 1:45]))
X <- X[, -(1:6)] # odrzucam 1. chromosom X
X[X == 9] <- NA
X <- replaceNA(X) - 1
n <- nrow(X)
M <- ncol(X)
# Y
#y <- as.matrix(c(d1[, 46], d2[, 46]), ncol = 1)
y <- as.matrix(c(d1[, 46], d2[, 46], d3[, 46], d4[, 46]), ncol = 1)
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
                    c(-0.2,0.2),
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
plot(res$a.vec[-1],colMeans(exp(res$pi.gibbs[burnin:n.gibbs,])),type='l',xlab='beta',xlim=c(-0.1,0.1),ylab='pi(beta)')
EX = sum(res$a.vec[-1]*colMeans(exp(res$pi.gibbs[burnin:n.gibbs,])))
EX2 = sum(res$a.vec[-1]^2*colMeans(exp(res$pi.gibbs[burnin:n.gibbs,])))
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
window.size=5
Betas.np <- colMeans(res$beta.gibbs[burnin:n.gibbs,])
res.np <- plot.graph(Betas.np,window.size, markers, beta_true = NULL, x.sim =res$beta.gibbs[burnin:n.gibbs,] )
fig.np <- res.np$fig
fig.np <- fig.np +  labs(title=paste('MSE = ',round(sqrt(mean(res.np$res^2)),3) , " (window size = ", window.size,")",sep=""))
fig.np <- fig.np+ theme(plot.title = element_text(hjust = 0.5))

print(fig.np)

#ggsave('np_real_w5.jpeg',fig.np)
window.size=1
Betas.np <- colMeans(res$beta.gibbs[burnin:n.gibbs,])
res.np <- plot.graph(Betas.np,window.size, markers, beta_true = NULL, x.sim =res$beta.gibbs[burnin:n.gibbs,] )
fig.np <- res.np$fig
fig.np <- fig.np +  labs(title=paste('MSE = ',round(sqrt(mean(res.np$res^2)),3) , " (window size = ", window.size,")",sep=""))
fig.np <- fig.np+ theme(plot.title = element_text(hjust = 0.5))
#ggsave('np_real_w1.jpeg',fig.np)


CI <- t(apply(exp(res$pi.gibbs[burnin:n.gibbs,]), 2, quantile, probs=c(0.05,0.5,0.95)))
df <- data.frame(x = res$a.vec[-1], pi = CI[,2])
fig <- ggplot(df) + geom_line(aes(x=x, y=pi))
fig <- fig  + geom_ribbon(data = data.frame(x = res$a.vec[-1],
                                            low = CI[,1]
                                            ,upp = CI[,3]),
                          aes(x=x,ymin = low, ymax=upp), alpha=0.2, fill='blue')
