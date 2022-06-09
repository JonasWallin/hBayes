
library(NPBayes)

n <- 100
p <- 20
sigma <- 1
X <- matrix(rnorm(n*p), n, p)
beta <- 1*rnorm(p)
beta_true <- sample(beta)
Y <- X%*%beta_true + sigma*rnorm(n)

res <- gibbs.normal.permute.beta.fixed.sigma(100,
                                 Y,
                                X,
                                sigma,
                                c(min(beta)-0.1,max(beta)+0.1),
                                8,
                                beta,
                                F)


