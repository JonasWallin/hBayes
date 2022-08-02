
rm(list=ls())
save.fig=F
graphics.off()
library(NPBayes)
n.gibbs=100000
n=100
X = matrix(rnorm(n*n),nrow=n,ncol=n)
y = rnorm(n)

res <- gibbs.normal(n.gibbs,
                    y,
                    X,
                    c(-4,4),
                    8,
                    rep(0,n))

pi.cdf      <- rbind(0,apply(exp(res$pi.gibbs),1,cumsum))
pi.quant    <- apply(pi.cdf,1,quantile,c(0.025,0.50,0.975))


plot(res$a.vec,pi.quant[2,],col = 1,lwd = 2, main = "",ylab = "CDF",xlab = expression(beta), type='l')
lines(res$a.vec,pi.quant[1,],col=1,lty = 2,lwd = 2)
lines(res$a.vec,pi.quant[3,],col=1,lty = 2,lwd = 2)
