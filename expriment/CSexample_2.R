rm(list=ls())
library(NPBayes)
set.seed(1)
L	<- 6
I	<- 2^L

log.odds <- function(p) {
    log(p / (1 - p))
}
inv.log.odds <- function(theta)	{
    exp(theta) / (1 + exp(theta))
}

n.gibbs <- 200
n <- 400
kappa <- 0.1 #p/n
gamma2 <- 5 #gamma^2
eps <- 1/2 #czesc zmiennych niezerowych

p <- round(kappa*n) #liczba zmiennych
k <- round(eps*p) #liczba zmiennych z niezerowymi wspolczynnikami
gamma <- sqrt(gamma2) #gamma
amplitude <- sqrt(gamma2)/(sqrt(kappa*eps)) #wartosci niezerowych wspolczynnikow

beta.1 <- c(rep(amplitude,times = k), rep(0,times = p - k)) #wektor wspolczynnikow, p-k z nich to 0, a k to 'amplitude'
X <- matrix(rnorm(n*p)/sqrt(n), n, p)
Y		<- rbinom(n,1,inv.log.odds(X %*% beta.1))
fit.1	<- glm(Y ~ X + 0,family = binomial) #change later
res <- gibbs.logistic(n.gibbs,
                      Y,
                      X,
                      c(-24,24),
                      L,
                      fit.1$coefficients)
res.Yekutieli <- gibbs.logistic.Yekutieli(n.gibbs,
                                          Y,
                                          X,
                                          c(-24,24),
                                          L,
                                          fit.1$coefficients)


par(mfrow=c(2,1))
plot(res$beta.gibbs[,1],ylim=c(-24,24))
points(res.Yekutieli$beta.gibbs[,1], col='red')
plot(res$delta.gibbs[,1],ylim=c(1,55))
points(res.Yekutieli$delta.gibbs[,1], col='red')
