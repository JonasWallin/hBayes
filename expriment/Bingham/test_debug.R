lambda <- c(3,1)
a <- 1/lambda
b <- -4*0.5*lambda
a_ <- (a[1]*b[1]+a[2]*b[2])
b_ <- (a[1]*b[2]+a[2]*b[1])
c <- a_ - b_
theta.seq <- seq(0,2*pi, length.out=200)
n <- 100000
theta.vec <- rep(0,n)
for(i in 1:n){
    theta.vec[i] <- rangle(-c )

    U <- runif(1)
    if(U < 0.25)
    {
        theta.vec[i] <- pi-theta.vec[i]
    }else if(U < 0.5){
        theta.vec[i] <-  2*pi - theta.vec[i]
    }else if(U < 0.75){
        theta.vec[i] <- theta.vec[i] + 1.*pi
    }

}

par(mfrow=c(1,1))
res.h<- hist(theta.vec,freq = F,breaks = theta.seq)
dens.e <- exp(c * (1-sin(theta.seq)^2))
C <- sum(dens.e)*(theta.seq[2]-theta.seq[1])

lines(theta.seq, dens.e/C,col='red')
