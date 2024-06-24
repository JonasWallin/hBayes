#'
#' generating figure exploring acceptance rate for 2x2 Bingham,
#' and compare to Hoff (2009)
#'
library(NPBayes)
graphics.off()
sim <- 10^3
a.grid <- c(0.00001,10) #seq(2,5*10^4,length.out = 10)
rbing.2diagmatrix.count <- function( a, b){
    a_ <- (a[1]*b[1]+a[2]*b[2])
    b_ <- (a[1]*b[2]+a[2]*b[1])
    a_b <- a_ - b_
    theta <- rangle(-a_b, count.sample=T)

    return(theta[2])
}
# taken from rstiefel
rbing.Op <-function (A, B)
{
    b <- diag(B)
    bmx <- max(b)
    bmn <- min(b)
    if (bmx > bmn) {
        A <- A * (bmx - bmn)
        b <- (b - bmn)/(bmx - bmn)
        vlA <- eigen(A)$val
        diag(A) <- diag(A) - vlA[1]
        vlA <- eigen(A)$val
        nu <- max(dim(A)[1] + 1, round(-vlA[length(vlA)]))
        del <- nu/2
        M <- solve(diag(del, nrow = dim(A)[1]) - A)/2
        rej <- TRUE
        cholM <- chol(M)
        nrej <- 0
        while (rej) {
            Z <- matrix(rnorm(nu * dim(M)[1]), nrow = nu, ncol = dim(M)[1])
            Y <- Z %*% cholM
            tmp <- eigen(t(Y) %*% Y)
            U <- tmp$vec %*% diag((-1)^rbinom(dim(A)[1], 1, 0.5))
            L <- diag(tmp$val)
            D <- diag(b) - L
            lrr <- sum(diag((D %*% t(U) %*% A %*% U))) - sum(-sort(diag(-D)) *
                                                                 vlA)
            rej <- (log(runif(1)) > lrr)
            nrej <- nrej + 1
        }
    }
    return(nrej)
}
b <- c(2,0.001)
rej.new <- rej.rstiefel <- rep(0,length(a.grid))
for(i in 1:length(a.grid)){

    a <- c(a.grid[i],10000)

    E <- rbing.2diagmatrix(a, b)
    for(j in 1:sim){
        rej.new[i] <- rej.new[i] + rbing.2diagmatrix.count(-a,b)
        rej.rstiefel[i] <- rej.rstiefel[i] + rbing.Op(-E%*%diag(a)%*%t(E),diag(b))
    }
    rej.new[i] <-rej.new[i]/sim
    rej.rstiefel[i] <- rej.rstiefel[i]/sim
}

plot(a.grid, rej.rstiefel, type='l')
lines(a.grid, rej.new, col='red')
