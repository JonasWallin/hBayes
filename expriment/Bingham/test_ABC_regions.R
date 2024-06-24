set.seed(1)
graphics.off()
S <- 1
a <- -1*c(20, 0.1)
b <- c(0.2,0.1)
C <- matrix(c(-1,-1,50,-1.1),2,2)
#interval 0,pi/2

ab.above.cos <- cbind(c(pi/2,1,-3,-3*pi/2) ,c(-1,-2/pi,2/pi,1))
ab.below.cos <- cbind(c(1,pi/2,-3*pi/2,-3),c(-2/pi,-1,1,2/pi))
ab.above.sin <- cbind(c(0,pi,2,-4),c(1,-1,-2/pi,2/pi))
ab.below.sin <- cbind(c(0,2,pi,-2*pi),c(2/pi,-2/pi,-1,1))

pi.interval <- c(0,pi/2,pi,3/2*pi,2*pi)

theta.vec <- seq(0,2*pi,length.out=200)#seq(pi.interval[interval], pi.interval[interval+1], length.out=100)
#'
#' density of general ..
#'
#'
dens <- function(theta,S, a, b, C){
    Z_1 <- cos(theta)
    Z_2 <- sin(theta)
    E_1 <- c(Z_1,Z_2)
    E_2 <- c(S*E_1[2],-S*E_1[1])
    E <- cbind(E_1,E_2)
    B <- diag(b,nrow = 2)
    A <- diag(a,nrow = 2)
    a_ <- (a[1]*b[1]+a[2]*b[2])
    b_ <- (a[1]*b[2]+a[2]*b[1])
    c_1 <- a_ - b_

    return( exp((sum(diag(t(C)%*%E)+diag(B%*%t(E)%*%A%*%E)))))
}
dens2 <- function(theta,S, a, b, C){
    a_ <- (a[1]*b[1]+a[2]*b[2])
    b_ <- (a[1]*b[2]+a[2]*b[1])
    c_1 <- a_ - b_
    c_2s <- C[1,1] - C[2,2]*S
    c_3s <- C[2,1] + C[1,2]*S
    if( c_1 < 0){
        return(exp(c_1*cos(theta)^2 + c_2s*cos(theta) +c_3s*sin(theta)))
    }else{
        return(exp(-c_1*sin(theta)^2 + c_2s*cos(theta) +c_3s*sin(theta) + c_1) )
    }
}

dens <- Vectorize(dens, "theta")
dens2 <- Vectorize(dens2, "theta")
#plot(theta.vec,dens(theta.vec,-1, a, b, C)+dens(theta.vec,1, a, b, C),type='l',ylim=c(0,2))
d1 <- dens(theta.vec,S, a, b, C)
d2 <- dens2(theta.vec,S, a, b, C)
plot(theta.vec,d1,type='l',col='red',ylim=c(0,max(d1,d2)))
#plot(theta.vec,dens(theta.vec,1, a, b, C),type='l',col='red')
lines(theta.vec,d2,type='l',col='green')
a_ <- (a[1]*b[1]+a[2]*b[2])
b_ <- (a[1]*b[2]+a[2]*b[1])
c_1 <- a_ - b_
c_2s <- C[1,1] - C[2,2]*S
c_3s <- C[2,1] + C[1,2]*S
Norm_C <- c(0,0,0,0)
dens_ <- rep(0,length(theta.vec))
dens_2 <- rep(0,length(theta.vec))
sigma2 <- mu <- Norm_C
for(interval in c(1:4)){
    ind <- theta.vec >= pi.interval[interval] & theta.vec <= pi.interval[interval+1]
    if(c_1 < 0){
        if(interval%in%c(2,3)){
            a_1 <- ab.above.cos[interval, 1]
            b_1 <- ab.above.cos[interval, 2]
        }else{
            a_1 <- ab.below.cos[interval, 1]
            b_1 <- ab.below.cos[interval, 2]
        }
    }else{
        if(interval%in%c(3,4)){
            a_1 <- ab.above.sin[interval, 1]
            b_1 <- ab.above.sin[interval, 2]
        }else{
            a_1 <- ab.below.sin[interval, 1]
            b_1 <- ab.below.sin[interval, 2]
        }
    }
    if(c_2s < 0){
        a_2 <- ab.below.cos[interval, 1]
        b_2 <- ab.below.cos[interval, 2]
    }else{
        a_2 <- ab.above.cos[interval, 1]
        b_2 <- ab.above.cos[interval, 2]
    }
    if(c_3s < 0){
        a_3 <- ab.below.sin[interval, 1]
        b_3 <- ab.below.sin[interval, 2]
    }else{
        a_3 <- ab.above.sin[interval, 1]
        b_3 <- ab.above.sin[interval, 2]
    }

    sigma2[interval] <- 1/(2*abs(c_1)*b_1^2)
    mu[interval] <- (-2*abs(c_1)*b_1*a_1 + c_2s*b_2+c_3s*b_3) * sigma2[interval]
    #truncated normal
    ProbI <- (pnorm(pi.interval[interval+1],mu[interval],sqrt(sigma2[interval]),log=F) - pnorm(pi.interval[interval],mu[interval],sqrt(sigma2[interval])) )
    d_ <- dnorm(theta.vec[ind],mean = mu[interval], sd = sqrt(sigma2[interval]))/ProbI

    Norm_C[interval] <- ProbI*exp(c_2s*a_2+c_3s*a_3 - abs(c_1)*a_1^2 + 0.5*(c_2s*b_2+c_3s*b_3-2*abs(c_1)*a_1*b_1)^2*sigma2[interval])*sqrt(2*pi*sigma2[interval])
    d_2 <- exp(-abs(c_1) * (a_1 + b_1*theta.vec[ind])^2 + c_2s*(a_2+b_2*theta.vec[ind])+ c_3s*(a_3+b_3*theta.vec[ind]))
    dens_[ind] <- d_*Norm_C[interval]
    dens_2[ind] <- d_2


}
lines(theta.vec, dens_,col='yellow')
lines(theta.vec, dens_2,col='blue',lty = 3)
abline(v=pi.interval)
