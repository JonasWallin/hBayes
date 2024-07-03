#' Comparing samplers ordering
#'
#'
library(NPBayes)
graphics.off()
save.fig <- T
values <- c(100,10,5,1,rep(0.1,20))
n <- 100
eigA <- list(values = 1/values,
     vectors = diag(nrow=length(values),ncol=length(values)))
eigB <-  list(values = -values*n/2,
             vectors = diag(nrow=length(values),ncol=length(values)))

p <- length(values)
n.mcmc <- 20000 * p
burnin <- 1000*p
res.iter <- rbing.iter(n.mcmc/p, NULL, NULL,eigA = eigA, eigB = eigB,ret = 2)
res      <- rbing(n.mcmc, NULL, NULL,eigA = eigA, eigB = eigB,ret = 2)
lag <- 50
r <- acf(abs(res.iter$Es[1,1,seq(burnin,n.mcmc,by=p)]),lag)
if(save.fig)
    pdf('looping1.pdf')
plot(r,main=paste('ACF for looping Gibbs ', expression(X[11]),sep=''))
if(save.fig)
    dev.off()

r <-   acf(abs(res$Es[1,1,seq(burnin,n.mcmc,by=p)]),lag)
if(save.fig)
    pdf('random1.pdf')
plot(r,main=paste('ACF for random Gibbs ', expression(X[11]),sep=''))
if(save.fig)
    dev.off()

r <- acf(abs(res.iter$Es[3,3,seq(burnin,n.mcmc,by=p)]),lag)
if(save.fig)
    pdf('looping3.pdf')
plot(r,main=paste('ACF for looping Gibbs ', expression(X[33]),sep=''))
if(save.fig)
    dev.off()

r <-   acf(abs(res$Es[3,3,seq(burnin,n.mcmc,by=p)]),lag)
if(save.fig)
    pdf('random3.pdf')
plot(r,main=paste('ACF for random Gibbs ', expression(X[33]),sep=''))
if(save.fig)
    dev.off()

r <- acf(abs(res.iter$Es[15,15,seq(burnin,n.mcmc,by=p)]),lag)
if(save.fig)
    pdf('looping15.pdf')
plot(r,main=paste('ACF for looping Gibbs ', expression(X[1515]),sep=''))
if(save.fig)
    dev.off()

r <-   acf(abs(res$Es[15,15,seq(burnin,n.mcmc,by=p)]),lag)
if(save.fig)
    pdf('random15.pdf')
plot(r,main=paste('ACF for random Gibbs ', expression(X[1515]),sep=''))
if(save.fig)
    dev.off()


vec <- rep(0, length(values))
for(j in 2:length(values)){
        ind.1 <- (res$inds[1,]==1 & res$inds[2, ]==j) |(res$inds[1,]==j & res$inds[2,]==1)
        vec[j]<-  mean(1-abs(res$Zs[1,ind.1]))
}
if(save.fig)
    pdf('Z1.pdf')
plot(2:length(values),vec[-1], xlab='i',ylab='1-|Z|')
if(save.fig)
    dev.off()
vec <- rep(0, length(values))
for(j in 1:length(values)){
    if(j==15)
        next
    ind.1 <- (res$inds[1,]==15 & res$inds[2, ]==j) |(res$inds[1,]==j & res$inds[2,]==15)
    vec[j]<-  mean(1-abs(res$Zs[1,ind.1]))
}
if(save.fig)
    pdf('Z2.pdf')
plot(setdiff(1:length(values),15),vec[-15], xlab='i',ylab='1-|Z|')
if(save.fig)
    dev.off()
