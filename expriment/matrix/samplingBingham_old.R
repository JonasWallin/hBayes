library(rstiefel)
library(matrixcalc)
library(truncnorm)
graphics.off()
set.seed(4)

sim <- 100
sim2 <- 200
p <- 2
A <- diag((p:1))#diag((p:1)^10)
E.A <- eigen(A)
B<-diag(sort(rexp(p),decreasing=TRUE))
E<-rbing.Op(-0.5* A,B)
Es <- rstiefel::rbing.matrix.gibbs(-0.5* A,B,E)


loglik_ <- -0.5* (sum(diag(B%*%t(E)%*%A%*%E)))
loglik_2 <- -0.5* ( t(vec(E))%*%kronecker(B,A)%*%vec(E))
loglik_3 <- -0.5* ( t(vec(t(E.A$vectors)%*%E))%*%kronecker(B,diag(E.A$values))%*%vec(t(E.A$vectors)%*%E))
print(loglik_ - loglik_2)
print(loglik_3 - loglik_2)

dbing <- function(E,A,B){
    return(-0.5* (sum(diag(B%*%t(E)%*%A%*%E))))
}



#'
#' The exp(tr( diag(b)%*%t(E)%*%A%*%E) )
#' We paramatrise the switch as
#' E = E%*%Z
#' Z = [cos(theta), s sin(theta);
#'      sin(theta), -s cos(theta) ];
#' gives
#' p(Z) = exp(-0.5*tr( diag(b)%*%t(Z)%*% (t(E)%*%A%*%E ) %*%Z) )
#'      = exp(-0.5* t(vec(Z))%*%kron(diag(b), t(E)%*%A%*%E) %*%vec(Z))
#'      make the transformation
#'      Z' = (t(E)%*%A%*%E ) =  t(E_A') D_A E_A
#'      sample ZE_A' which is unitary matrix with density
#'      = exp(-0.5* t(vec(Z'))%*%kron(diag(b), D_A) %*%vec(Z'))
#'      = exp(-0.5 (a*cos(theta)^2 + b sin(theta)^2))
#'      = exp(-0.5 * (a+b cos(theta)^2) )
#'      i.e.
#' Taken from D.Hoff 2009
#' @param E - (p x 2) swap
#' @param A - (p x 1) paramaeter
#' @param b - (2 x 1) B diagonal values of b
rbing.shuffle <- function(E , A, b){

    AE1 <- A*E[,1]
    AE2 <- A*E[,2]
    EtAE <- matrix(0,2,2)
    EtAE[1,1] <- sum(AE1*E[,1] )
    EtAE[2,2] <- sum(AE2*E[,2] )
    EtAE[1,2] <- sum(AE1*E[,2] )
    EtAE[2,1] <- EtAE[1,2]
    #theta
}

dbing.theta <- function(theta,A,B){
    E <- cbind(c(cos(theta),sin(theta)),c(sin(theta),-cos(theta)))
    return(-0.5* (sum(diag(B%*%t(E)%*%A%*%E))))
}
dbing.theta2 <- function(theta,A,B){
    E <- cbind(c(cos(theta),sin(theta)),c(sin(theta),-cos(theta)))
    return(exp(-0.5* t(vec(E))%*%kronecker(B,A)%*%vec(E) ))
}

#' very special two dim Bingham
#' @param theta
dtheta.bing <- function(theta, a_b){
    return(exp(-0.5*a_b * cos(theta)^2))
}
rbing2 <-function (A, B, a = NULL, E = NULL)
{
    if (is.null(a)) {
        trA <- A[1, 1] + A[2, 2]
        lA <- 2 * sqrt(trA^2/4 - A[1, 1] * A[2, 2] + A[1, 2]^2)
        a <- lA * (B[1, 1] - B[2, 2])
        E <- diag(2)
        if (A[1, 2] != 0) {
            E <- cbind(c(0.5 * (trA + lA) - A[2, 2], A[1, 2]),
                       c(0.5 * (trA - lA) - A[2, 2], A[1, 2]))
            E[, 1] <- E[, 1]/sqrt(sum(E[, 1]^2))
            E[, 2] <- E[, 2]/sqrt(sum(E[, 2]^2))
        }
    }
    b <- min(1/a^2, 0.5)
    beta <- 0.5 - b
    lrmx <- a
    if (beta > 0) {
        lrmx <- lrmx + beta * (log(beta/a) - 1)
    }
    lr <- -Inf
    count <- 0
    while (lr < log(runif(1))) {
        w <- rbeta(1, 0.5, b)
        lr <- a * w + beta * log(1 - w) - lrmx
        count <- count + 1
    }
    cat('samples=',count,'\n')
    u <- c(sqrt(w), sqrt(1 - w)) * (-1)^rbinom(2, 1, 0.5)
    x1 <- E %*% u
    x2 <- (x1[2:1] * c(-1, 1) * (-1)^rbinom(1, 1, 0.5))
    return(cbind(x1, x2))
}


#'sample theta
#' TODO: sampling negative a_b
# improve using cos(2\theta)/2 \propto  cos(\theta)^2 on [3/2*pi, 2pi]
rtheta.bing <- function(a_b){
    # sample f(\theta). =exp(-0.5*a_b * cos(theta)^2) on pi/2 -> pi
    # f(pi/2) = exp(0) = 1,
    # f(acos(sqrt(-log(0.5)/(0.5*a_b)))) = 0.5
    # f(pi) = 0
    x_1 <- pi/2
    c_ <- -2*log(0.5)/a_b
    if(c_ < 1){
        x_2 <- acos(-sqrt(c_)) #position f(x)=0.5
        if(x_2-x_1 <0.5){
            x_3 <- min(acos(-sqrt(-log(x_2-x_1)/(0.5*a_b))),3*pi/4)
        }else{
            x_3 <- 3*pi/4
        }
        g_2 <- 0.5
    }else{
        x_2 <- pi/2+0.1
        g_2 <- dtheta.bing(x_2,a_b)
        x_3 <- 3*pi/4
    }


    #ensure that the minimum probability of p_1>p_3
    x_4 <- pi
    p_1 =  (x_2-x_1)
    p_2 = g_2 * (x_3-x_2)
    g_3 <-  dtheta.bing(x_3,a_b)
    p_3 =  g_3 * (x_4-x_3)
    p_sum <- (p_1+p_2 + p_3)
    p_2 <- p_2/p_sum
    p_1 <- p_1/p_sum
    while(TRUE){
        U <- runif(1)
        V <- runif(1)
        if(U < p_1){
            theta<- V*(x_2 - x_1) + x_1
            g <- 1
        }else if(U<p_1+ p_2){
            theta <- V*(x_3 - x_2) + x_2
            g <- g_2
        }else{
            theta <- V*(x_4 - x_3) + x_3
            g <- g_3
        }
        if(g*runif(1) <  dtheta.bing(theta,a_b))
            break
    }
    U <- runif(1)
    if(U < 0.25)
    {
        theta <- pi-theta
    }else if(U < 0.5){
        theta <-  2*pi - theta
    }else if(U < 0.75){
        theta <- theta + 1.*pi
    }


    return(theta)
}
#' compute the eigenvalues and eigenvectors of 2x2 symmetric matrix
#'
eigen2by2 <- function(A){

}

#sample the Bingham distribution for neg def diagonal (B or A) pos def diagonal(B or A)
#' @param E - (p x p) inital guess
#' @param A - (p x 1) eigenvalues of A
#' @param B - (p x 1) eigenvalues of B
rbing.diagmatrix <- function( A, B){
    a <- (A[1]*B[1]+A[2]*B[2])
    b <- (A[1]*B[2]+A[2]*B[1])
    a_b <- a - b
    theta <- rtheta.bing(a_b)
    S <- 2 * (runif(1)< 0.5) - 1
    # for possible values are theta that is uniform
    # sample s
    # build E
    Z_1 <- cos(theta)
    Z_2 <- sin(theta)
    E_1 <- c(Z_1,Z_2)
    E_2 <- c(S*E_1[2],-S*E_1[1])
    return(cbind(E_1,E_2))
}


# compute eigenvalue eigenector of 2d matrix
# todo create sampler for negative and positive eigenvalues
# (should reduce down to negative and positive a_b) should be easy to create by shifting
#
Ev <- array(dim=c(dim = c(dim(Es),sim)))
Ev[,,1] <-  Es
for(i in 2:sim)
    Ev[,,i] <- rstiefel::rbing.matrix.gibbs(-0.5* A,B,Ev[,,i-1])
cat('sd_emp[1,] = ',round(apply(Ev[1,,], 1, var),2),'\n')
cat('sd_emp[2,] = ',round(apply(Ev[2,,], 1, var),2),'\n')

#par(mfrow=c(2,1))
hist(Ev[1,1,],50,freq=F)

Es_2 <- rbing.matrix(diag(A),diag(B))
Ev_2 <- array(dim=c(dim = c(dim(Es_2),sim)))
Ev_2[,,1] <-  Es_2
for(i in 2:sim)
    Ev_2[,,i] <- rbing.matrix(diag(A),diag(B))


hist(Ev_2[1,2,],50,freq=F)
cat('sd_2_emp[1,] = ',round(apply(Ev_2[1,,], 1, var),2),'\n')
cat('sd_2_emp[2,] = ',round(apply(Ev_2[2,,], 1, var),2),'\n')
lik_1 <- rep(0,sim)
lik_2 <- rep(0,sim)
for(i in 1:sim){
    lik_1[i] <- dbing(Ev[,,i],A,B)
    lik_2[i] <- dbing(Ev_2[,,i],A,B)
}
plot(lik_1)
plot(lik_2)
theta_vec <- seq(0,2, length.out=sim2)
lik <- lik2 <- rep(0,sim2)
a <- (A[1,1]*B[1,1]+A[2,2]*B[2,2])
b <- (A[1,1]*B[2,2]+A[2,2]*B[1,1])
a_b <- a - b
for(i in 1:sim2){
    lik[i]  <- dtheta.bing(pi*theta_vec[i],-a_b)
    lik2[i] <- dbing.theta2(pi*theta_vec[i], A, B)
}

theta_s <- rep(0,sim2)
for(i in 1:sim2)
    theta_s[i] <- rtheta.bing(a_b)
res <- hist(theta_s,100, freq = F)
lines(pi*theta_vec,lik/(sum(lik)*diff(pi*theta_vec)[1]))
lines(pi*theta_vec,lik2/(sum(lik2)*diff(pi*theta_vec)[1]),col='red')

x_1 <- pi/2
c_ <- -2*log(0.5)/a_b
if(c_ < 1){
    x_2 <- acos(-sqrt(c_)) #position f(x)=0.5
    if(x_2-x_1 <0.5){
        x_3 <- min(acos(-sqrt(-log(x_2-x_1)/(0.5*a_b))),3*pi/4)
    }else{
        x_3 <- 3*pi/4
    }
    g_2 <- 0.5
}else{
    x_2 <- pi/2+0.1
    g_2 <- dtheta.bing(x_2,a_b)
    x_3 <- 3*pi/4
}
rbing2(-0.5* A,B)


if(0){
    a_b_seq <- exp(seq(log(10^-10),log(1),length.out=1000))
    p_vec <- rep(0, length(a_b_seq))
    for(i in 1:length(a_b_seq)){
        a_b <- a_b_seq[i]
        #ensure that the minimum probability of p_1>p_3
        x_1 <- pi/2
        c_ <- -2*log(0.5)/a_b
        if(c_ < 1){
            x_2 <- acos(-sqrt(c_)) #position f(x)=0.5
            if(x_2-x_1 <0.5){
                x_3 <- min(acos(-sqrt(-log(x_2-x_1)/(0.5*a_b))),3*pi/4)
            }else{
                x_3 <- 3*pi/4
            }
            g_2 <- 0.5
            #exp(-0.5*a_b*cos(x_3)^2+(0.5)*(1-1/(a_b*(x_3-x_1)^2))*a_b*(x_3-x_1)^2)
        }else{
            x_2 <- pi/2+0.1
            g_2 <- dtheta.bing(x_2,a_b)
            x_3 <- 3*pi/4
        }
        x_4 <- pi
        p_1 =  (x_2-x_1)
        p_2 = g_2 * (x_3-x_2)
        g_3 <-  dtheta.bing(x_3,a_b)
        p_3 =  g_3 * (x_4-x_3)
        p_sum <- (p_1+p_2 + p_3)
        p_2 <- p_2/p_sum
        p_1 <- p_1/p_sum
        p_vec[i] <- p_1+p_2
    }
    plot(a_b_seq,p_vec)
    cat("p_1+p_2 = ",round(p_1+p_2,4),'\n')
}
