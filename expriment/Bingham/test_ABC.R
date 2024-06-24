set.seed(1)
graphics.off()
a <- -c(2, 0.1)
b <- c(0.2,0.1)
C <- 0.1*matrix(c(0.5,-1,2,-1),2,2)

dens <- function(theta,S, a, b, C){
    Z_1 <- cos(theta)
    Z_2 <- sin(theta)
    E_1 <- c(Z_1,Z_2)
    E_2 <- c(S*E_1[2],-S*E_1[1])
    E <- cbind(cbind(E_1,E_2))
    B <- diag(b,nrow = 2)
    A <- diag(a,nrow = 2)
    return( exp((sum(diag(C%*%E)+diag(B%*%t(E)%*%A%*%E)))))
}


#theta.vec <- seq(pi,3*pi/2, length.out=100)
theta.vec <- seq(0,2*pi, length.out=1000)

lik <- rep(0, length(theta.vec))
for(i in 1:length(lik))
    lik[i] <- dens(theta.vec[i],1,a,b,C)

#' approximate cos by linear function (above)
#' returns just the a+b\theta for each interval [0,\pi/2),[pi/2,pi),[pi,3/2*pi),[3/2*pi, 2*pi]
#' @return  ab - (4 x 2)
ab.above.cos <- cbind(c(pi/2,1,-3,-3*pi/2) ,c(-1,-2/pi,2/pi,1))

c_1 <- (2/pi)^2
abc.above.cos <- cbind(c(1,-(pi/2)^2,0,0),c(((pi/2)^2*c_1-1)*2/pi,pi,0,0),c(-c_1,-1,0,0))
abc.below.cos <- cbind(c(1,0,0,0),c(0,0,0,0),c(-1/2,0,0,0))

#' approximate cos by linear function (below)
#' returns just the a+b\theta for each interval [0,\pi/2),[pi/2,pi),[pi,3/2*pi),[3/2*pi, 2*pi]
ab.below.cos <- cbind(c(1,pi/2,-3*pi/2,-3),c(-2/pi,-1,1,2/pi))


#' approximate sin by linear function (above)
#' returns just the a+b\theta for each interval [0,\pi/2),[pi/2,pi),[pi,3/2*pi),[3/2*pi, 2*pi]
ab.above.sin <- cbind(c(0,pi,2,-4),c(1,-1,-2/pi,2/pi))


#' approximate sin by linear function (below)
#' returns just the a+b\theta for each interval [0,\pi/2),[pi/2,pi),[pi,3/2*pi),[3/2*pi, 2*pi]
ab.below.sin <- cbind(c(0,2,pi,-2*pi),c(2/pi,-2/pi,-1,1))

ab.cos.approx <- function(theta, ab){

    interval <- c(0,pi/2,pi,3/2*pi,2*pi,Inf)
    res <- rep(0, length(theta))
    for(i in 1:4){
        ind <- interval[i] <= theta & theta <= interval[i+1]
        res[ind] <- ab[i,1] + ab[i,2] * theta[ind]
    }
    return(res)
}

abc.cos.approx <- function(theta, abc){

    interval <- c(0,pi/2,pi,3/2*pi,2*pi,Inf)
    res <- rep(0, length(theta))
    for(i in 1:4){
        ind <- interval[i] <= theta & theta <= interval[i+1]
        res[ind] <- abc[i,1] +  abc[i,2] * theta[ind] + abc[i,3] * theta[ind]^2
    }
    return(res)
}

par(mfrow=c(1,2))
#plot(theta.vec,lik,type='l',xlab='theta',main='cos')
b <- 0
plot(theta.vec,cos(theta.vec)+b,type='l',ylim=c(-1.5,1.5),xlab='theta',main='cos')
lines(theta.vec,ab.cos.approx(theta.vec, ab.above.cos),col='red')
lines(theta.vec,ab.cos.approx(theta.vec, ab.below.cos),col='blue')
lines(theta.vec,abc.cos.approx(theta.vec, abc.above.cos),col='red',lty=2)
lines(theta.vec,abc.cos.approx(theta.vec, abc.below.cos),col='blue',lty=2)

b <- 0
plot(theta.vec,sin(theta.vec)+b,type='l',ylim=c(-1.5,1.5),xlab='theta',main='sin')
lines(theta.vec,ab.cos.approx(theta.vec, ab.above.sin),col='red')
lines(theta.vec,ab.cos.approx(theta.vec, ab.below.sin),col='blue')

#lines(theta.vec,(pi/2)*(1-(2/pi)*theta.vec)+b,col='blue')


#lines(theta.vec,(-1+(2/pi)*(theta.vec-pi))+b,col='red')


#plot(theta.vec,sin(theta.vec),type='l',ylim=c(0,2))
#lines(theta.vec,((2/pi)*theta.vec),col='red')
#lines(theta.vec,(pi/2)*((2/pi)*theta.vec),col='blue')

