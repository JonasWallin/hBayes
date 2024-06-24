

################################

# 1.    Aux. functions

################################


stein.loss <- function(sigma.hat,sigma)
{
  m <- dim(sigma)[1]
  sigma.prod <-sigma.hat %*% solve(sigma)

  return( sum(diag(sigma.prod)) - log(det(sigma.prod)) - m)
}

# Quick stein loss computation for Sigma and Sigma.hat  with same D.oracle

oracle.stein.loss <- function(gamma.hat,gamma)
{
  p <- dim(gamma)[1]
  sigma.prod <- gamma.hat %*% diag(D.oracle) %*% t(gamma.hat) %*% gamma %*% diag(1/D.oracle) %*% t(gamma)

  return( sum(diag(sigma.prod)) - log(det(sigma.prod)) - p)
}

#	 Stein (1961) estimator

stein.est <- function(x.smp)
{
  n <- dim(x.smp)[1]
  p <- dim(x.smp)[2]
  
  S <- t(x.smp) %*% x.smp
  K <- chol(S)
  Delta.vec <- 1 / (1+n+p-2*(1:p))
  
  return( t(K) %*% diag(Delta.vec) %*% K )
}




#	1.1 Functions used for generating random orthogonal matrix using METHOD B in Page 22 of Diaconis and Shahshahani (1987)

acos.2pi <- function(cos.rho,sin.rho)
{
  if(0 <= sin.rho) return(acos(cos.rho))
  if(0 > sin.rho) return(2*pi - acos(cos.rho))
}

eval.orth.mat.2D	<- function(orth.mat)	c(acos.2pi(orth.mat[1,1],orth.mat[1,2]), orth.mat[2,2] / orth.mat[1,1])

R.rho.b	<- function(rho,b,m = NULL)
{
  if(is.null(m)) gamma.tmp <- cbind(c(cos(rho),-b * sin(rho)),c(sin(rho),b * cos(rho)))

  if(!is.null(m))
  {
    gamma.tmp <- diag(rep(1,m))
    gamma.tmp[(m-1):m,(m-1):m] <- cbind(c(cos(rho),-b * sin(rho)),c(sin(rho),b * cos(rho)))
  }
  return(gamma.tmp)
}


smp.v.old	<- function(p)
{
  rz	<- rnorm(p,mean = 0, sd = 1)
  return(rz / sqrt(sum(rz^2)))
}


smp.near.v <- function(vv,sd.smp = 0.01)
{
  p <- length(vv)
  rz <- rnorm(p,mean = vv, sd = sqrt(sd.smp^2/p))
  return(rz / sqrt(sum(rz^2)))
}


house.reflctn	<- function(v.p,m = NULL)
{
  p	<- length(v.p)
  if(is.null(m))
  {
    e1	<- c(1,rep(0,p-1))
    v.m <- v.p
    m   <- p
  }
  if(!is.null(m))
  {
    e1	<- c(rep(0,m-p),1,rep(0,p-1))
    v.m <- c(rep(0,m-p),v.p)
  }
  x.m	<- (e1 - v.m) / sqrt(sum((e1 - v.m)^2))
  return( diag(rep(1,m)) - 2* x.m %*% t(x.m))
}


smp.rnd.orth.mat	<- function(n = 10)
{

  P.tmp	<- R.rho.b(runif(1,0,2*pi),1-2*rbinom(1,1,0.5))

  if(2 < n)
    for(i in 3:n)
    {
      v.i		<- smp.v(i)
      hr		<- house.reflctn(v.i)
      P.tmp	<- hr %*% cbind(c(1,rep(0,i-1)),rbind(0,P.tmp))
    }

  return(P.tmp)
}


smp.rnd.orth.mat.vp 	<- function(p = 10)
{

  b.smp <- 1-2*rbinom(1,1,0.5)
  rho.smp <- runif(1,0,2*pi)

  P.tmp	<- R.rho.b(rho.smp,b.smp)
  vp.vec <- c(rho.smp,b.smp)

  if(2 < p)
    for(i in 3:p)
    {
      v.i		<- smp.v(i)
      hr		<- house.reflctn(v.i)
      P.tmp	<- hr %*% cbind(c(1,rep(0,i-1)),rbind(0,P.tmp))

      vp.vec <- c(v.i,vp.vec)
    }

  return(list(P.tmp,vp.vec))
}

gamma.vp.dist <- function(vp.vec.1,vp.vec.2)
{
  k <- length(vp.vec.1)
  m <- (-1 + sqrt(1 + 4*2*(k+1))) / 2
  v.p.mat <- array(dim = c(m,m))
  ind0.vec <- c(m:3,1,1)
  ind1.vec <- cumsum(c(1,m:4,3,1))
  ind2.vec <- cumsum(c(m:3,1,1))

  dist <- 0
  for(i in 1:(m-2)) dist <- dist + sum((vp.vec.1[ind1.vec[i]:ind2.vec[i]] - vp.vec.2[ind1.vec[i]:ind2.vec[i]])^2)*i
  dist <- dist + abs(vp.vec.1[ind1.vec[m-1]] - vp.vec.2[ind1.vec[m-1]])/(2*pi)

  return(dist)
}


trns.v.p.mat.vec <- function(v.p.mat)
{
  m <- dim(v.p.mat)[1]
  v.p.vec <- v.p.mat[1,(m-1):m]
  for(i in (m-2):1) v.p.vec <- c(v.p.mat[1:(m+1-i),i],v.p.vec)
  return(v.p.vec)
}

trns.v.p.vec.mat <- function(v.p.vec.mat)
{
  k <- length(v.p.vec)
  m <- (-1 + sqrt(1 + 4*2*(k+1))) / 2
  v.p.mat <- array(dim = c(m,m))
  ind0.vec <- c(m:3,1,1)
  ind1.vec <- cumsum(c(1,m:4,3,1))
  ind2.vec <- cumsum(c(m:3,1,1))

  for(i in 1:(m-2)) v.p.mat[1:ind0.vec[i],i] <- v.p.vec[ind1.vec[i]:ind2.vec[i]] / sqrt(sum(v.p.vec[ind1.vec[i]:ind2.vec[i]]^2))
  v.p.mat[1:ind0.vec[m-1],m-1] <- v.p.vec[ind1.vec[m-1]:ind2.vec[m-1]] %% (2*pi)
  v.p.mat[1:ind0.vec[m],m] <- sign(v.p.vec[ind1.vec[m]:ind2.vec[m]])

  return(v.p.mat)
}


#	1.2 Express rotation matrix in terms of sequence of house.matrix v.n vectors

# Express orthogonal matrix as sequence of v.p vectors


vp.repr.orth.mat <- function(orth.mat)
{
  m 		<- dim(orth.mat)[1]
  v.p.mat	<- array(dim=c(m,m))

  if(m == 2) dimnames(v.p.mat) <- list(NULL,c("rho","b"))

  if(2 < m)
  {
    dimnames(v.p.mat) <- list(NULL,c(paste("v.",m:3,sep=""),"rho","b"))
    for(i in 1:(m-2))
    {
      v.p.mat[1:(m-i+1),i]	<- 	orth.mat[,1]
      orth.mat					<- (house.reflctn(v.p.mat[1:(m-i+1),i]) %*% orth.mat)[2:(m-i+1),2:(m-i+1)]
    }
  }

  v.p.mat[1,c(m-1,m)] <- 	 eval.orth.mat.2D(orth.mat)
  return(v.p.mat)
}



#	1.3 Construct gamma matrices from sequence of v.p vectors


const.gamma	<- function(v.p.mat)
{
  m <- dim(v.p.mat)[1]

  gamma.rght <- R.rho.b(v.p.mat[1,m-1],v.p.mat[1,m],m)

  if(2 < m) for(i in (m-2):1) gamma.rght <- house.reflctn(v.p.mat[1:(m-i+1),i],m)	%*% gamma.rght

  return(gamma.rght)
}

# 1.4  Compute marginal likelihood for p-dimensional DS decomposition vector v.p for 3 <= p <= m
#       where eigenvector matrix is Gamma  and thus x.smp %*% Gamma ~ N(0, diag(D.vec))
#       Assumes data matrix is x.smp


gamma.lambda.loglik <- function(gamma.tmp,lambda.vec)
{
  x.trns	    <- x.smp %*% gamma.tmp
  return(sum(dnorm(x.trns,mean = rep(0,p*n),sd = rep(sqrt(lambda.vec),each = n),log = T)))
}


oracle.gamma.x.loglik <- function(x.tmp,gamma.tmp, D)
{
  x.trns	    <- x.tmp %*% gamma.tmp
  return(sum(dnorm(x.trns,mean = rep(0,p*n),sd = rep(sqrt(D),each = n),log = T)))
}

