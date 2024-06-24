library(rstiefel)
set.seed(2)
sim <- 10000
if(0){
Z<-matrix(rnorm(4),2,2) ; A<-t(Z)%*%Z
A <- matrix(c(2,0,0,0.1),2,2)
#A <- diag(sort(rexp(2),decreasing=TRUE))
B<-diag(sort(rexp(2),decreasing=TRUE))
U<-rbing.Op(A,B)
Es <- rbing(sim,A,B)

Es <- array(0, dim=c(2,2,sim))
Us <- array(0, dim=c(2,2,sim))
for(i in 1:sim){
    U<-rbing.Op(A,B)
    Us[,,i] <- U
    E <- rbing.2(A, diag(B))
    Es[,,i] <- E
}
for(i in 1:2){
    for(j in 1:2){
        cat('mean(E-U[',i,',',j,'] ) = (',round(mean(Es[i,j,]-Us[i,j,])-2*sd(Es[i,j,]-Us[i,j,])/sqrt(sim),4),',',round(mean(Es[i,j,]-Us[i,j,])+2*sd(Es[i,j,]-Us[i,j,])/sqrt(sim),4),')\n')

    }
}

Es2 <- rbing(sim,A,B)
}
if(1){
    A <- diag(c(3,2,1))#toeplitz(c(3,1,1))
    B <- diag(c(3,2,1))
    U<-rbing.Op(A,B)
    Es <- rbing(sim,A,B)

    Us <- array(0, dim=c(3,3,sim))
    for(i in 1:sim){
        U<-rbing.matrix.gibbs(A,B,U)
        Us[,,i] <- U
    }
    for(i in 1:2){
        for(j in 1:2){
            cat('mean(E-U[',i,',',j,'] ) = (',round(mean(Es[i,j,]-Us[i,j,])-2*sd(Es[i,j,]-Us[i,j,])/sqrt(sim),4),',',round(mean(Es[i,j,]-Us[i,j,])+2*sd(Es[i,j,]-Us[i,j,])/sqrt(sim),4),')\n')

        }
    }

}
