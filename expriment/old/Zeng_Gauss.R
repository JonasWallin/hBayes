library(NPBayes)
library(glmnet)
#set.seed(1)

burnin <- 2000
n.gibbs <- 4000
L	<- 8
I	<- 2^L
real.Y=T
#set to soruce location
replaceNA <- function(X) {
    ## replace missing values with the mean
    ## (function uses replaceOneNA)
    na_pos <- apply(X, 1, function(x) which(is.na(x)))
    for (i in seq_along(na_pos)) {
        if (length(na_pos[[i]]) == 0) next
        X[i, na_pos[[i]]] <- sapply(na_pos[[i]], replaceOneNA, x = X[i, ])
    }
    return(X)
}

replaceOneNA <- function(na, x) {
    ## replace one missing value
    if (!is.na(x[na])) stop("Value is not NA.")
    n <- length(x)
    left <- which(!is.na(x[na:1]))[1] - 1
    right <- which(!is.na(x[na:n]))[1] - 1
    if (is.na(left) & is.na(right)) stop("All values are NA.")
    if (is.na(left)) left <- Inf else if (is.na(right)) right <- Inf
    if (left < right) {
        replace <- x[na - left]
    } else if (left > right) {
        replace <- x[na + right]
    } else {
        replace <- (x[na - left] + x[na + right])/2
    }
    return(replace)
}

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

n <- dim(X)[1]
beta.1 <- rep(0, dim(X)[2])
#beta.1[sample(1:length(beta.1),5)] <- 10
beta.1[1] <- 1
if(real.Y==F){
    Y		<- X%*%beta.1  + rnorm(n);
}else{
    Y <- as.matrix(c(d1[, 46], d2[, 46], d3[, 46], d4[, 46]), ncol = 1)
    #y <- as.matrix(d1[, 46], ncol = 1)
    #Y <- as.matrix(c(d1[, 46], d2[, 46]), ncol = 1)
    Y <- scale(Y)
}
X <- scale(X)
beta0 <- glmnet(X, Y, alpha = 0, lambda = 0,standardize=F,intercept=F)$beta
res <- gibbs.normal(n.gibbs,
                      Y,
                      X,
                      c(-1,1),
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
plot(res$a.vec[-1],colMeans(exp(res$pi.gibbs[burnin:n.gibbs,])),type='l',xlab='beta',xlim=c(-0.25,0.25),ylab='pi(beta)')
EX = sum(res$a.vec[-1]*colMeans(exp(res$pi.gibbs[burnin:n.gibbs,])))
EX2 = sum(res$a.vec[-1]^2*colMeans(exp(res$pi.gibbs[burnin:n.gibbs,])))
cat('E[beta] = ',round(EX,2),'\n')
cat('V[beta] = ',round(EX2-EX^2,2),'\n')
