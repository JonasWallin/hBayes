##
# power simulation for mixed effect model
# to create Table 3: set.seed(2) q=1000, table =2 ,gives sparse mixed effect
# to create Table 4: set.seed(2) q=1000, table =1 ,gives laplace distribution
# D: 2020-02-10
##



#install.packages(c("bigstep","RcppEigen","ggplot2"))
#install.packages(c("bigstep","glmnet"))
#install.packages(c("devtools"))

#library(devtools)

#install_github("JonasWallin/PolyMixed")


library(bigstep)
library(PolyMixed)
library(RcppEigen)
library(ggplot2)
library(bigstep)
library(glmnet)
library(NPBayes)

set.seed(2)

# parameters----------------------------------------------------------------

n.gibbs <- 500
burnin  <- 200
L <- 7
chr <- 10                # chromosomes
m <- 150                 # markers on chromosome
M <- chr * m             # all markers
markers <- rep(m, chr)   # markers on each chromosome
n         <- 400         # observations
d <- 1                   # distance between markers (cM)
sigma <- 1               # sd for error term
tau <- 0.01              # sd for polygenic effect
mu <-  0.004             # mean for polygenic effect


qtl.pos <- c(151, 824, 1274)

beta.true <- 0.5


X <- simulateDesignMatrix(chr, m, n, d)
X <- scale(X)
Z <- rnorm(dim(X)[2])
beta_          <- rep(0, dim(X)[2])
beta_[qtl.pos] <- beta.true
y <- X%*%(beta_ + mu[1] + tau *Z) + sigma * rnorm(dim(X)[1])

beta0 <- glmnet(X, y, alpha = 0, lambda = 0,standardize=F,intercept=F)$beta
res <- gibbs.normal(n.gibbs,
                    y,
                    X,
                    c(-4,4),
                    L,
                    as.vector(beta0),cpp=T)
par(mfrow=c(3,2))
plot(res$a.vec[-1],colMeans(exp(res$pi.gibbs)),type='l',xlab='beta',ylab='pi(beta)')
EX = sum(res$a.vec[-1]*colMeans(exp(res$pi.gibbs[burnin:n.gibbs,])))
EX2 = sum(res$a.vec[-1]^2*colMeans(exp(res$pi.gibbs[burnin:n.gibbs,])))
cat('E[beta] = ',round(EX,2),'\n')
cat('V[beta] = ',round(EX2-EX^2,2),'\n')
hist(res$beta.gibbs[burnin:n.gibbs,qtl.pos[1]])
hist(res$beta.gibbs[burnin:n.gibbs,qtl.pos[2]])
hist(res$beta.gibbs[burnin:n.gibbs,qtl.pos[3]])
hist(res$beta.gibbs[burnin:n.gibbs,1400])
plot(colMeans(res$beta.gibbs[burnin:n.gibbs,]))
