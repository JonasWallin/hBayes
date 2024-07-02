library(truncnorm)

#'sample the Bingham distribution for diagonal matrices A,B
#' f(E|A,B) \propto exp(B%*%t(E)%*%A%*%E)I( t(E)%*%E=I)
#' @param  a - (2 x 1) diagonal of A
#' @param  b - (2 x 1) diagonal of B
#' @return E - (2 x 2) orthonormal matrix
rbing.2diagmatrix <- function( a, b){
    a_ <- (a[1]*b[1]+a[2]*b[2])
    b_ <- (a[1]*b[2]+a[2]*b[1])
    a_b <- a_ - b_
    theta <- rangle(-a_b)

    U <- runif(1)
    if(U < 0.25)
    {
        theta <- pi-theta
    }else if(U < 0.5){
        theta <-  2*pi - theta
    }else if(U < 0.75){
        theta <- theta + 1.*pi
    }

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

#'
#' Sampling exp(-c cos(\theta)^2) I(\theta \in [pi/2,pi])
#' @param c            - ( 1 x 1) constant can be negative
#' @param count.sample - (bool) count number of accept reject
#' @return theta       - (1 x 1 2 x 1) sample theta, if sample.count=T return number rejections
rangle <- function(c, count.sample=FALSE){

    count <- 0
    if(c==0){
        theta <- pi/2 *(1+ runif(1))
        if(count.sample)
            return(c(theta, count))
        return(theta)
    }
    neg.c <- FALSE
    if(c < 0)
    {
        neg.c <- TRUE
        c <- -c
    }
    sd.c <- pi/(2 *sqrt(2*c) )
    v.c <- sd.c^2
    while(TRUE){
        theta <- truncnorm::rtruncnorm(1,
                                       a = pi/2,
                                       b = pi,
                                       mean = pi/2,
                                       sd   = sd.c)
        U <-  runif(1)
        if(log(U) <  -c * cos(theta)^2 -( - (theta - pi/2)^2/(2*v.c))){
            if(neg.c)
                theta <- (3/2*pi-theta)
            if(count.sample)
                return(c(theta, count))
            return(theta)
        }
        count <- count + 1
    }
}

#' propto exp(-c cos(\theta)^2) I(\theta \in [pi/2,pi])
#' does not give the normalizing constant
#' @param c - ( 1 x 1) constant can be negative
dangle <- function(theta, c, log=F ){
    if(log)
        return(-c * cos(theta)^2)
    return(exp(-c * cos(theta)^2))
}


eigen_2x2 <- function(A) {

    if(A[1,2] != 0){
        # Calculate the trace and determinant of the matrix
        tr <- A[1,1] + A[2,2]

        # Compute the eigenvalues
        C <- sqrt(4*A[1,2]^2 + (A[1,1]-A[2,2])^2)
        eigenvalues <- 0.5 * c((tr + C), (tr - C))

        C1 <- eigenvalues[1] - A[1,1]
        E1 <- c(A[1,2], C1)
        E1 <- E1/sqrt(E1[1]^2 + E1[2]^2)
        E2 <- c(-E1[2], E1[1])
    }else{
        if(A[1,1]>A[2,2]){
            eigenvalues <- c(A[1,1],A[2,2])
            E1 <- c(1,0)
            E2 <- c(0,1)
        }else{
            eigenvalues <- c(A[2,2],A[1,1])
            E2 <- c(1,0)
            E1 <- c(0,1)
        }
    }
    return(list(values = eigenvalues, vectors= matrix(c(E1,E2),2,2)))
}


#'sample the Bingham distribution for a 2x2 for matrices A,and diagonal B
#' f(E|A,B) \propto exp(B%*%t(E)%*%A%*%E)I( t(E)%*%E=I)
#' @param  A - (2 x 2) A symmetric matrix
#' @param  b - (2 x 1) diagonal of B
#' @return E - (2 x 2) orthonormal matrix
rbing.2 <- function( A, b){
    Eigen.A <- eigen_2x2(A)
    E <- rbing.2diagmatrix(Eigen.A$values, b)
    return(Eigen.A$vectors%*%E)
}


#' Gibbs sampler Bingham distribution (pxp) for matrices A, and B
#' f(E|A,B) \propto exp(B%*%t(E)%*%A%*%E)I( t(E)%*%E=I)
#' @param  n   - (int) how iteration samples
#' @param  A   - (p x p) A symmetric matrix
#' @param  B   - (p x p) B symmetric matrix
#' @param eigA - (list) eigenvalues eigenvector output of eigen(A)
#' @param eigB - (list) eigenvalues eigenvector output of eigen(B)
#' @param E0   - (p x p) initial guess
#' @param EtAE - (p x p) E0^T eigen(A) E0
#' @param ret. - (1 x 1) 1 - return E, 2 - E_A E E^T_B
#' @return E   - (p x p) orthonormal matrix
rbing <- function(n, A, B, eigA=NULL, eigB=NULL, E0=NULL, EtAE=NULL, ret = 1){

    if(is.null(eigA))
        eigA = eigen(A)
    if(is.null(eigB))
        eigB = eigen(B)
    p <- dim(eigA$vectors)[1]

    EtAE_update <- matrix(0, nrow=p, ncol=p)
    if(is.null(E0)){
        E0 <- diag(x=1,nrow=p)
        EtAE = diag(eigA$values)
        EtAE_update <- EtAE_update + 1
    }else{
        if(is.null(EtAE)){
            EtAE = matrix(0, nrow=p, ncol=p)
        }else{
            EtAE_update <- EtAE_update + 1
        }
    }
    #compute E^T \Lambda_A E (which if)
    E <- E0
    Es <- array(0, dim=c(p,p,n))
    if(ret==1){
        tE_B <- t(eigB$vectors)
    }else{
    }
    Zs   <- array(0,dim=c(2,n))
    inds <- array(0,dim=c(2,n))
    for(i in 1:n){
        ind <- sample(p,2)
        if(EtAE_update[ind[1],ind[1]]==0){
            EtAE[ind[1],ind[1]] <- sum(E[,ind[1]]^2 * eigA$values)
            EtAE_update[ind[1],ind[1]] <- 1
        }
        if(EtAE_update[ind[1],ind[2]]==0){
            EtAE[ind[1],ind[2]] <- sum(E[,ind[1]] * E[,ind[2]]  * eigA$values)
            EtAE[ind[2],ind[1]] <- EtAE[ind[1],ind[2]]
            EtAE_update[ind[1],ind[2]] <- 1
            EtAE_update[ind[2],ind[1]] <- 1
        }
        if(EtAE_update[ind[2],ind[2]]==0){
            EtAE[ind[2],ind[2]] <- sum(E[,ind[2]]^2 * eigA$values)
            EtAE_update[ind[2],ind[2]] <- 1
        }

        Z <- rbing.2(EtAE[ind, ind], eigB$values[ind])
        Zs[,i] <- Z[1,]
        inds[,i] <- ind
        EtAE[ind, ind] <- t(Z)%*%EtAE[ind, ind]%*%Z
        EtAE_temp <- EtAE[-ind, ind]%*%Z
        EtAE[-ind, ind] <- EtAE_temp
        EtAE[ind, -ind] <- t(EtAE_temp)
        E[,ind] <- E[,ind]%*%Z
        if(ret==1){
            Es[,,i] <- eigA$vectors%*%E%*%tE_B
        }else{
            Es[,,i] <- E
        }
    }

return(list(Es=Es,Zs=Zs,inds=inds))
}


#' Gibbs sampler Bingham distribution (pxp) for matrices A, and B
#' f(E|A,B) \propto exp(B%*%t(E)%*%A%*%E)I( t(E)%*%E=I)
#' here which index is looped
#' @param  n.p   - (int) how iteration samples (multiple of p) (p * n.p)
#' @param  A   - (p x p) A symmetric matrix
#' @param  B   - (p x p) B symmetric matrix
#' @param eigA - (list) eigenvalues eigenvector output of eigen(A)
#' @param eigB - (list) eigenvalues eigenvector output of eigen(B)
#' @param E0   - (p x p) initial guess
#' @param EtAE - (p x p) E0^T eigen(A) E0
#' @param ret. - (1 x 1) 1 - return E, 2 - E_A E E^T_B
#' @return E   - (p x p) orthonormal matrix
rbing.iter <- function(n.p, A, B, eigA=NULL, eigB=NULL, E0=NULL, EtAE=NULL, ret = 1){

    if(is.null(eigA))
        eigA = eigen(A)
    if(is.null(eigB))
        eigB = eigen(B)
    p <- dim(eigA$vectors)[1]

    EtAE_update <- matrix(0, nrow=p, ncol=p)
    if(is.null(E0)){
        E0 <- diag(x=1,nrow=p)
        EtAE = diag(eigA$values)
        EtAE_update <- EtAE_update + 1
    }else{
        if(is.null(EtAE)){
            EtAE = matrix(0, nrow=p, ncol=p)
        }else{
            EtAE_update <- EtAE_update + 1
        }
    }
    #compute E^T \Lambda_A E (which if)
    E <- E0
    Es <- array(0, dim=c(p,p,n.p*p))
    if(ret==1){
        tE_B <- t(eigB$vectors)
    }else{
    }
    Zs   <- array(0,dim=c(2,n.p*p))
    inds <- array(0,dim=c(2,n.p*p))
    for(i in 1:(n.p*p)){
        ind <- c(1 + (i-1)%%(p-1),
                 2 + (i-1)%%(p-1) )
        if(EtAE_update[ind[1],ind[1]]==0){
            EtAE[ind[1],ind[1]] <- sum(E[,ind[1]]^2 * eigA$values)
            EtAE_update[ind[1],ind[1]] <- 1
        }
        if(EtAE_update[ind[1],ind[2]]==0){
            EtAE[ind[1],ind[2]] <- sum(E[,ind[1]] * E[,ind[2]]  * eigA$values)
            EtAE[ind[2],ind[1]] <- EtAE[ind[1],ind[2]]
            EtAE_update[ind[1],ind[2]] <- 1
            EtAE_update[ind[2],ind[1]] <- 1
        }
        if(EtAE_update[ind[2],ind[2]]==0){
            EtAE[ind[2],ind[2]] <- sum(E[,ind[2]]^2 * eigA$values)
            EtAE_update[ind[2],ind[2]] <- 1
        }

        Z <- rbing.2(EtAE[ind, ind], eigB$values[ind])
        Zs[,i] <- Z[1,]
        inds[,i] <- ind
        EtAE[ind, ind] <- t(Z)%*%EtAE[ind, ind]%*%Z
        EtAE_temp <- EtAE[-ind, ind]%*%Z
        EtAE[-ind, ind] <- EtAE_temp
        EtAE[ind, -ind] <- t(EtAE_temp)
        E[,ind] <- E[,ind]%*%Z
        if(ret==1){
            Es[,,i] <- eigA$vectors%*%E%*%tE_B
        }else{
            Es[,,i] <- E
        }
    }

    return(list(Es=Es,Zs=Zs,inds=inds))
}
