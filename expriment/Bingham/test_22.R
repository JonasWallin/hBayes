graphics.off()
sim <- 50000
theta.seq <- seq(0,2*pi, length.out=200)
E.seq <- seq(-1,1,length.out=200)
A <- matrix(c(0.9402562,  21.9719868,  21.9719868, 999.4983204),2,2)
A <- A + t(A)
A <- A/2
EA <- eigen(A)
EA$values[1] <- 10
EA$values[2] <- 0.1
A <- (EA$vectors)%*%diag(EA$values)%*%t(EA$vectors)
#A <- diag(EA$values)
B <- diag(c(5e-04, 5e-01))




dens.theta <- rep(0, length(theta.seq))
dens.e <- dens.theta
for(i in 1:length(theta.seq)){
    dens.theta[i] <- dbing_2x2(theta.seq[i],A,B)
    dens.e[i]     <- dbing_2x2(E.seq[i],A,B,Ebasis=TRUE)
}
E.sim <- rep(0, sim)
for(i in 1:sim){
    E <-rbing.2(A,diag(B))
    E.sim[i] <- E[1,1]
}
theta.sim <- acos(E.sim) + pi*(runif(sim)<0.5)

C <- sum(dens.theta)*(theta.seq[2]-theta.seq[1])
par(mfrow=c(1,2))
res.h<- hist((theta.sim),freq = F,breaks = theta.seq)
lines(theta.seq,dens.theta/C,col='red')

C.e <- sum(dens.e)*(E.seq[2]-E.seq[1])
res.h<- hist((E.sim),freq = F,breaks = E.seq)
lines(E.seq,dens.e/C.e,col='red')
