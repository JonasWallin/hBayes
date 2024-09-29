set.seed(1)
library(NPBayes)
library(glmnet)
L	<- 10
a.dom <- c(-5,5)
n.gibbs <- 500
p <- 50
n <- 100
X <- matrix(rnorm(n*p), n, p)
U <- rep(0, p)
U[sample(p,5)] <- 1
beta <- rnorm(p)*(U-1) + U * abs(rnorm(p,mean = 4,sd =1))

X <- scale(X)
y <- X%*%beta +1*rnorm(n)
beta0 <- glmnet(X, y, alpha = 0, lambda = 0,standardize=F,intercept=F)$beta
res <- gibbs.normal(n.gibbs,
                    y,
                    X,
                    a.dom,
                    L,
                    as.vector(beta0),cpp=F)
